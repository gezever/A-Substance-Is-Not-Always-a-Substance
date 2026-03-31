# ==============================================================================
# 13a_classyfire_coverage.R
# Analysis 13: ChemOnt class coverage via local ClassyFire cache
#
# PURPOSE
# -------
# ClassyFire (Djoumbou Feunang et al. 2016, J. Cheminform. 8:61,
# DOI: 10.1186/s13321-016-0174-y) assigns every chemical structure to a
# position in the ChemOnt ontology (ClassyFire 2.1) based on its InChIKey.
# This analysis asks which ChemOnt `direct_parent` classes contain the most
# regulatory substances — i.e., which chemical families are most heavily
# regulated across European instruments.
#
# The analysis uses a local cache of ClassyFire API responses rather than
# querying the API at runtime, making it reproducible without network access.
# The cache is populated by the companion script
# `13b_before_prioritisation_create_scheme.R`.
#
# DATA PROVENANCE
# ---------------
# Input 1: `data/processed/all_substances.rds`
#          Provides (inchikey, source) pairs for all structure-defined entries.
# Input 2: `data/cache/classyfire/*.json`
#          One JSON file per InChIKey, containing the full ClassyFire
#          classification including `direct_parent.chemont_id` and
#          `direct_parent.name`.  Files produced by the ClassyFire REST API
#          at https://cfb.fiehnlab.ucdavis.edu/entities/{inchikey}.
#
# METHODOLOGY
# -----------
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **ChemOnt 2.1 / ClassyFire** (Djoumbou Feunang et al. 2016, *ClassyFire: automated chemical classification with a comprehensive, computable taxonomy*, J. Cheminform. 8:61. DOI: 10.1186/s13321-016-0174-y) | Hierarchical chemical ontology providing a `direct_parent` class for each InChIKey via structure-based rules.  The cache contains JSON responses from the ClassyFire REST API (https://cfb.fiehnlab.ucdavis.edu). |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | `direct_parent` as the aggregation level | One level below the kingdom in ChemOnt; offers sufficient specificity to discriminate chemical families while remaining broadly interpretable.  Superclass/kingdom are too coarse (groups all organic acids together); subclass/class fragment data into hundreds of near-empty bins.  Validated by confirming that the top classes correspond to well-known regulatory families (organohalogens, fatty acids, steroids). |
# | Coverage counted on unique InChIKeys (not records) | Avoids double-counting a substance that appears in multiple regulatory sources; the metric measures chemical family breadth, not administrative repetition. |
# | Top-30 bar chart with colour = n_sources | Bar length captures the absolute size of each class in the regulatory dataset; colour adds a second dimension (regulatory breadth) without requiring a separate plot. |
#
# INTERPRETATION
# --------------
# **Analysis 13a — Top 30 ChemOnt classes by substance coverage**
# Long bars indicate chemical families that are heavily regulated in absolute
# terms.  Dark colour (high n_sources) means the family is regulated across
# many different instruments simultaneously; light colour (low n_sources)
# means regulation is concentrated in one or two instruments.
#
# A family with a long bar but light colour may be an opportunity for
# cross-instrument harmonisation: many substances in the same chemical class
# but regulated by only one instrument.  A family with a short bar but dark
# colour is a small but broadly regulated class — likely a legacy or
# high-concern family.
#
# **Analysis 13b — Distribution of substance coverage across classes**
# A highly right-skewed histogram (many classes with 1–5 substances, a few
# with hundreds) is expected and reflects the power-law structure of
# chemical space.  The long right tail consists of the dominant families
# shown in 13a.  If the histogram is approximately uniform, the ChemOnt
# cache is likely incomplete (few InChIKeys resolved) and coverage should
# be increased before drawing conclusions.
#
# OUTPUTS
# -------
# output/figures/Analysis_13a_ChemOnt_prioritisation_top30.pdf
# output/figures/Analysis_13b_ChemOnt_coverage_histogram.pdf
# output/tables/Analysis_13_class_coverage.csv
# output/tables/Analysis_13a_top30.csv
# data/processed/class_coverage.rds  (consumed by 14_prioritization.R)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(purrr)
library(jsonlite)
library(readr)
library(scales)

# Null-safe fallback operator — returns b when a is NULL or length-0
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 13: Build ChemOnt coverage table from ClassyFire cache
# ==============================================================================
# INTENT
# Determine which chemical families (direct_parent classes) are most
# represented in European regulatory lists by counting how many unique
# InChIKeys — and how many distinct sources — each class covers.
# High coverage across many sources indicates a chemical family that is a
# regulatory priority across multiple instruments simultaneously.
# ==============================================================================

classyfire_cache_dir <- here("data", "cache", "classyfire")
cf_files <- list.files(classyfire_cache_dir,
                       pattern = "\\.json$", full.names = TRUE)

# Extract direct_parent from each cached ClassyFire response
cf_lookup <- map_dfr(cf_files, function(f) {
  parsed <- tryCatch(jsonlite::read_json(f), error = function(e) NULL)
  if (is.null(parsed) ||
      is.null(parsed[["inchikey"]]) ||
      is.null(parsed[["direct_parent"]])) return(NULL)
  tibble(
    inchikey          = sub("^InChIKey=", "", parsed[["inchikey"]]),
    direct_parent_id  = parsed[["direct_parent"]][["chemont_id"]] %||%
                          NA_character_,
    direct_parent_lbl = parsed[["direct_parent"]][["name"]]       %||%
                          NA_character_
  )
})

# Join with all_substances to get (inchikey, source) context
substance_class <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(inchikey, source) |>
  inner_join(cf_lookup, by = "inchikey")

# Count substances and source breadth per ChemOnt class
class_coverage <- substance_class |>
  filter(!is.na(direct_parent_id)) |>
  group_by(direct_parent_id, direct_parent_lbl) |>
  summarise(
    n_inchikeys = n_distinct(inchikey),   # unique substances in this class
    n_sources   = n_distinct(source),     # number of regulatory sources
    .groups     = "drop"
  ) |>
  arrange(desc(n_inchikeys))

# Save coverage table (used by 14_prioritization.R and reporting)
saveRDS(class_coverage, here("data", "processed", "class_coverage.rds"))

# Columns (Analysis_13_class_coverage.csv):
#   direct_parent_id  — ChemOnt identifier for the direct_parent class
#                       (e.g., CHEMONTID:0000000); primary key
#   direct_parent_lbl — human-readable class name (e.g., "Organochlorides")
#   n_inchikeys       — number of unique InChIKeys in all_substances assigned
#                       to this class; a substance appearing in multiple sources
#                       is counted only once
#   n_sources         — number of distinct regulatory sources containing at
#                       least one InChIKey from this class; ranges 1–n_sources
write_csv(class_coverage,
          here("output", "tables", "Analysis_13_class_coverage.csv"))

# ==============================================================================
# Analysis 13a: Top 30 ChemOnt classes by regulatory substance coverage
# ==============================================================================
# INTENT
# The top-30 bar chart gives a direct overview of which chemical families
# dominate European regulatory lists.  Colour encodes the breadth of
# regulatory coverage (number of distinct sources), revealing whether a class
# is intensively regulated by a single instrument or broadly covered across
# many instruments.  Both dimensions matter for prioritisation.
# ==============================================================================

top30_classes <- class_coverage |> slice_head(n = 30)

p13a <- ggplot(
  top30_classes,
  aes(x = reorder(direct_parent_lbl, n_inchikeys),
      y = n_inchikeys,
      fill = n_sources)
) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_inchikeys), hjust = -0.1, size = 3) +
  scale_fill_gradient(low = "#b3cde3", high = "#084594",
                      name = "Number of\nsources") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title    = "Analysis 13a: Top 30 ChemOnt classes by regulatory substance coverage",
    subtitle = "Colour = number of sources; bar length = number of unique InChIKeys",
    x        = NULL,
    y        = "Number of unique substances (InChIKey)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p13a)
ggsave(p13a,
       filename = here("output", "figures",
                       "Analysis_13a_ChemOnt_prioritisation_top30.pdf"),
       device = "pdf",
       height = 8, width = 12, units = "in")

# Columns (Analysis_13a_top30.csv):
#   direct_parent_id  — ChemOnt identifier (same as in Analysis_13_class_coverage.csv)
#   direct_parent_lbl — class name; the 30 entries with the highest n_inchikeys
#   n_inchikeys       — unique substances in this class (bar length in the figure)
#   n_sources         — regulatory breadth (colour encoding in the figure)
write_csv(top30_classes,
          here("output", "tables", "Analysis_13a_top30.csv"))

# ==============================================================================
# Analysis 13b: Distribution of substance coverage across ChemOnt classes
# ==============================================================================
# INTENT
# The histogram shows whether regulatory substances are concentrated in a few
# dominant chemical families or spread across many.  A highly skewed
# distribution (most classes cover only 1–2 substances; a few cover hundreds)
# would support a targeted approach focusing on the high-coverage classes.
# ==============================================================================

p13b <- ggplot(class_coverage, aes(x = n_inchikeys)) +
  geom_histogram(binwidth = 5, fill = "#4a90d9", colour = "white") +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Analysis 13b: Distribution of substance coverage across ChemOnt classes",
    subtitle = "How many classes cover 1, 5, 10 \u2026 substances?",
    x        = "Number of unique substances per ChemOnt class",
    y        = "Number of ChemOnt classes"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p13b)
ggsave(p13b,
       filename = here("output", "figures",
                       "Analysis_13b_ChemOnt_coverage_histogram.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

message(sprintf(
  "Analysis 13: %d ChemOnt classes with coverage; top class '%s' covers %d substances",
  nrow(class_coverage),
  class_coverage$direct_parent_lbl[1],
  class_coverage$n_inchikeys[1]
))

message("13a_classyfire_coverage.R: analysis 13 completed")
