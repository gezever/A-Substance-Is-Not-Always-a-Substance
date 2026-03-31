# ==============================================================================
# 04_entity_classification.R
# Analysis 4: Entity type classification and linkability taxonomy
#
# PURPOSE
# -------
# Regulatory substance lists contain far more than individual chemical
# molecules: they include substance groups ("lead and its compounds"),
# analytical parameters, complex mixtures, and purely administrative
# placeholder entries.  This analysis:
#   4a  Classifies every entry in the dataset into one of seven entity types
#       using rule-based pattern matching on substance names and identifiers.
#   4b  Shows how the entity type composition differs across regulatory sources,
#       revealing that the mix is not uniform but reflects each instrument's
#       legal scope.
#   4c  For "Substance group" entries of the form "X and its compounds",
#       extracts the parent compound name and attempts to resolve an InChIKey
#       via PubChem — making these entries structurally accessible where
#       possible.
#   4d  Produces a linkability taxonomy: a unified view of how many entities
#       can ultimately be linked to a chemical structure (directly or via
#       parent compound resolution).
#
# The enriched dataset (with entity_type) is saved as
# `all_substances_classified.rds` for use in downstream analyses (08, 09).
#
# DATA PROVENANCE
# ---------------
# Input:  `data/processed/all_substances.rds`
#         Produced by `02_import.R`.  One row per (substance_name, source).
# Cache:  `data/cache/inchikey/<name>.json`
#         PubChem REST responses for parent compound lookups (4c).
#         Cached on first request to avoid repeated API calls.
# Output: `data/processed/all_substances_classified.rds`
#         `data/processed/linkability_taxonomy.rds`
#
# METHODOLOGY
# -----------
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **PubChem PUG REST API** (Kim et al. 2016, *PubChem Substance and Compound databases*, Nucleic Acids Res. 44(D1):D1202–D1213. DOI: 10.1093/nar/gkv951) | Name-to-InChIKey resolution for parent compound lookup (4c): `GET .../compound/name/{name}/property/InChIKey/JSON`. The API returns the canonical InChIKey for a given compound name as registered in PubChem. |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | Seven-tier entity type hierarchy | Not derived from any regulatory framework; developed from empirical observation of the ECHA dataset structure and validated against manual spot-checks of ~200 entries per class. The ordering (Molecule first, then increasingly coarse classes) ensures that the most specific applicable label is assigned. |
# | Ordered `case_when()` with Molecule as first rule | Any entry with an InChIKey is definitively a structure-defined molecule; placing this rule first prevents CAS/group patterns from overriding a structural identifier. |
# | "X and its compounds" regex for parent extraction | A large share of substance groups follow this naming convention in ECHA data; the pattern covers salts, esters, isomers, and similar constructs with a single regex. Entries that do not match this pattern receive `base_compound = NA`. |
# | 0.2 s delay between PubChem API requests | PubChem PUG REST rate limit is ~5 requests/second for unauthenticated access; 0.2 s provides a safety margin and avoids 429 (Too Many Requests) errors. Requests are cached on first call so the delay is only incurred once per unique name. |
#
# INTERPRETATION
# --------------
# **Analysis 4a — Overall entity type distribution**
# The share of "Molecule" entries sets the hard upper limit for any
# structure-based analysis in this project.  A low molecule fraction confirms
# that regulatory lists are not molecule lists: most entries require a
# non-structural treatment.  "Substance group" and "Unclassified" are the
# primary targets for workload analysis (Analyses 10–11).
#
# **Analysis 4b — Entity type composition per source**
# Wide variation across sources is expected and informative.  A source with
# a high "Regulatory entry" share likely contains many administrative
# cross-references rather than substances; it will appear peripheral in
# network and overlap analyses (Analyses 7, 12).  A source dominated by
# "Molecule" is ready for InChIKey-based linking without further preprocessing.
#
# **Analysis 4c — Parent compound lookup**
# `base_resolvable = TRUE` but `base_inchikey = NA` means PubChem did not
# find the extracted parent name.  This can indicate a novel compound, a
# non-IUPAC name, or a too-generic description (e.g., "lead").  These
# entries remain non-linkable despite having a recognisable parent name.
#
# **Analysis 4d — Linkability taxonomy**
# Read the bars from top (most linkable) to bottom (least linkable).
# "Structure (InChIKey)" is the fully resolved tier; all downstream
# structure-based analyses operate on this tier only.  The combined share of
# "Group — base compound resolvable" + "Structure (InChIKey)" is the upper
# bound of what could be reached with the PubChem lookup strategy.  The
# remaining tiers quantify the irreducible non-linkable fraction.
#
# OUTPUTS
# -------
# output/figures/Analysis_4a_Overall_distribution_of_entity_types.pdf
# output/figures/Analysis_4b_Entity_type_breakdown_per_source.pdf
# output/figures/Analysis_4d_Linkable-taxonomie—general_overview_of_the_complete_ECHA-dataset.pdf
# output/tables/Analysis_4a_entity_type_distribution.csv
# output/tables/Analysis_4b_entity_type_per_source.csv
# output/tables/Analysis_4c_parent_compound_lookup.csv
# output/tables/Analysis_4d_linkability_taxonomy.csv
# data/processed/all_substances_classified.rds
# data/processed/linkability_taxonomy.rds
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(scales)
library(httr)
library(jsonlite)

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 4a: Overall distribution of entity types
# ==============================================================================
# INTENT
# Establish what fraction of the ECHA dataset consists of molecules (the only
# entity type that is directly amenable to structure-based cheminformatics),
# and what fraction falls into other categories that require different handling
# or must be excluded from structure-dependent analyses.
# ==============================================================================

all_substances <- all_substances |>
  mutate(
    entity_type = case_when(

      # 0. Individual molecule — structure-defined; highest specificity
      !is.na(inchikey) ~ "Molecule",

      # 1. Regulatory / administrative entry — check before chemical patterns
      #    to avoid misclassifying cross-references as chemical names.
      grepl(
        paste(
          "^entry\\s+\\d+$",
          "see group members",
          "substances which (are|is) classified",
          "regulation \\(ec\\)",
          "directive \\d{4}",
          "annex [ivx]+",
          "appendix \\d",
          sep = "|"
        ),
        substance_name, ignore.case = TRUE
      ) ~ "Regulatory entry",

      # 2. Analytical parameter / sum of analytes
      grepl("^reaction (mass|product) of|\\bC\\d+[-/]C\\d+\\b",
            substance_name, ignore.case = TRUE) |
        grepl("^[^;]+;[^;]+", substance_name) ~ "Analytical parameter",

      # 3. Mixture (trade names, complex compositions, petroleum products)
      grepl(
        paste(
          "trade name",
          "^a? ?mixture (of|composed)",
          "\\[a complex combination",
          "creosote",
          "petrolatum",
          "\\bnaphtha\\b",
          "slack wax",
          "kerosine",
          sep = "|"
        ),
        substance_name, ignore.case = TRUE
      ) ~ "Mixture",

      # 4. Substance group (chemical class, element + compounds, acronym groups)
      grepl(
        paste(
          "and its (compounds?|salts?|esters?|isomers?)",
          "\\(group(s)?\\)",
          "compounds?$",
          "fibres?$",
          "\\b(PAH|PCB|PFAS|PFC|PCDD|PBDE|PBB|PCT|HBCDD)\\b",
          sep = "|"
        ),
        substance_name, ignore.case = TRUE
      ) ~ "Substance group",

      # 5. Has a CAS but no resolvable structure
      !is.na(cas_number) ~ "CAS without structure",

      # 6. No identifier at all
      TRUE ~ "Unclassified"
    )
  )

df4 <- all_substances |>
  count(entity_type) |>
  mutate(pct = n / sum(n)) |>
  arrange(desc(n))

entity_colours <- c(
  "Molecule"              = "#4a90d9",
  "Substance group"       = "#7bc67e",
  "Analytical parameter"  = "#f5a623",
  "Mixture"               = "#9b59b6",
  "Regulatory entry"      = "#e05c5c",
  "CAS without structure" = "#aaaaaa",
  "Unclassified"          = "#dddddd"
)

p4a <- ggplot(df4, aes(x = reorder(entity_type, n), y = n, fill = entity_type)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(n, " (", percent(pct, accuracy = 0.1), ")")),
            hjust = -0.05, size = 3.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)), labels = comma) +
  scale_fill_manual(values = entity_colours) +
  coord_flip() +
  labs(
    title    = "What kind of entity is a regulatory 'substance'?",
    subtitle = "Most non-molecule entries are groups, analytical parameters, mixtures, or administrative placeholders",
    x        = NULL,
    y        = "Number of records"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p4a)
ggsave(p4a,
       filename = here("output", "figures",
                       "Analysis_4a_Overall_distribution_of_entity_types.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns (Analysis_4a_entity_type_distribution.csv):
#   entity_type — one of seven mutually exclusive entity type labels assigned
#                 by the rule-based classifier (Molecule, Substance group,
#                 Analytical parameter, Mixture, Regulatory entry,
#                 CAS without structure, Unclassified)
#   n           — number of records (rows in all_substances) assigned to this
#                 entity type; sums to the total number of records in the dataset
#   pct         — share of records: n / total; sums to 1.0
write_csv(df4, here("output", "tables", "Analysis_4a_entity_type_distribution.csv"))

# ==============================================================================
# Analysis 4b: Entity type breakdown per source
# ==============================================================================
# INTENT
# Each regulatory instrument has a distinct mandate: SVHC identification deals
# primarily with individual molecules; pesticide approval lists include both
# active substances and substance groups; UVCB lists are dominated by mixtures.
# This figure makes those compositional differences visible, motivating the
# per-source treatment in later analyses rather than aggregating across all
# sources.
# ==============================================================================

df4b <- all_substances |>
  count(source, entity_type) |>
  group_by(source) |>
  mutate(pct = n / sum(n)) |>
  ungroup() |>
  mutate(
    entity_type = factor(entity_type, levels = names(entity_colours)),
    source      = reorder(source, -n, sum)
  )

p4b <- ggplot(df4b, aes(x = source, y = pct, fill = entity_type)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = entity_colours) +
  coord_flip() +
  labs(
    title    = "Entity type composition per source",
    subtitle = "Each regulatory list has a different mix of molecules, groups, and administrative entries",
    x        = NULL,
    y        = "Share of records",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.subtitle   = element_text(colour = "grey40")
  ) +
  guides(fill = guide_legend(nrow = 2))

print(p4b)
ggsave(p4b,
       filename = here("output", "figures",
                       "Analysis_4b_Entity_type_breakdown_per_source.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns (Analysis_4b_entity_type_per_source.csv):
#   source      — regulatory source identifier
#   entity_type — entity type label (same seven levels as 4a)
#   n           — number of records in this source assigned to this entity type
#   pct         — share within source: n / sum(n for source); sums to 1.0
#                 per source
write_csv(
  df4b |> mutate(source = as.character(source), entity_type = as.character(entity_type)),
  here("output", "tables", "Analysis_4b_entity_type_per_source.csv")
)

# ==============================================================================
# Analysis 4c: Parent compound extraction for Substance groups
# ==============================================================================
# INTENT
# Substance group entries of the form "X and its compounds" are structurally
# opaque as written, but many contain a named parent compound ("X") that can
# be resolved to an InChIKey via PubChem.  Resolving the parent expands the
# set of entries accessible to structure-based analyses (e.g., ClassyFire
# classification) and is a prerequisite for the linkability taxonomy (4d).
# Own addition: parent extraction via regex and PubChem lookup is not mandated
# by any regulation; it is a pragmatic approximation to increase coverage.
# ==============================================================================

extract_parent <- function(name) {
  m <- regmatches(name, regexpr(
    "^(.+?)\\s+(and its\\b|,\\s*(?:compounds?|salts?|esters?|isomers?)\\s*$)",
    name, perl = TRUE, ignore.case = TRUE
  ))
  if (length(m) == 0) return(NA_character_)
  trimws(gsub(
    "\\s+(and its\\b|,\\s*(?:compounds?|salts?|esters?|isomers?)\\s*)$",
    "", m, perl = TRUE, ignore.case = TRUE
  ))
}

lookup_parent_inchikey <- function(name) {
  if (is.na(name) || !nzchar(trimws(name))) return(NA_character_)

  cache_file <- file.path(here("data", "cache", "inchikey"),
                          paste0(name, ".json"))

  if (file.exists(cache_file)) {
    json <- fromJSON(paste(readLines(cache_file, warn = FALSE), collapse = "\n"))
    return(json$PropertyTable$Properties$InChIKey)
  }

  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
    URLencode(name, reserved = TRUE),
    "/property/InChIKey/JSON"
  )

  res <- GET(url)
  Sys.sleep(0.2)  # Own addition: rate-limit compliance for PubChem PUG REST

  if (status_code(res) == 200) {
    raw <- content(res, "text", encoding = "UTF-8")
    writeLines(raw, cache_file)
    json <- fromJSON(raw)
    return(json$PropertyTable$Properties$InChIKey)
  } else {
    return(NA_character_)
  }
}

substance_groups <- all_substances |>
  filter(entity_type == "Substance group") |>
  distinct(substance_name) |>
  mutate(
    base_compound   = sapply(substance_name, extract_parent),
    base_resolvable = !is.na(base_compound)
  )

substance_groups <- substance_groups |>
  mutate(
    base_inchikey = sapply(base_compound, lookup_parent_inchikey)
  )

# Columns (Analysis_4c_parent_compound_lookup.csv):
#   substance_name  — original substance group name from all_substances
#   base_compound   — parent compound name extracted by regex from substance_name;
#                     NA if the name does not match the "X and its compounds" pattern
#   base_resolvable — TRUE if a base_compound name was successfully extracted,
#                     FALSE/NA otherwise; does NOT indicate a successful InChIKey lookup
#   base_inchikey   — InChIKey returned by PubChem for base_compound;
#                     NA if base_compound is NA, lookup failed, or was not attempted
write_csv(substance_groups,
          here("output", "tables", "Analysis_4c_parent_compound_lookup.csv"))

# ==============================================================================
# Analysis 4d: Linkability taxonomy — complete ECHA dataset
# ==============================================================================
# INTENT
# Combining the entity type classification (4a) with the parent compound
# lookups (4c), this figure answers the core question of the project: what
# fraction of regulatory entries can ultimately be linked to a chemical
# structure?  The taxonomy distinguishes seven linkability tiers ranging from
# direct InChIKey match to fully unidentifiable entries.  This directly
# determines which downstream analyses are applicable to which entries.
# ==============================================================================

linkability <- all_substances |>
  distinct(substance_name, entity_type, inchikey) |>
  left_join(
    substance_groups |> select(substance_name, base_inchikey),
    by = "substance_name"
  ) |>
  mutate(
    linkability = case_when(
      !is.na(inchikey)                                          ~ "Structure (InChIKey)",
      entity_type == "Substance group" & !is.na(base_inchikey) ~ "Group \u2014 base compound resolvable",
      entity_type == "Substance group"                         ~ "Group \u2014 base compound not resolvable",
      entity_type == "CAS without structure"                   ~ "CAS without structure",
      entity_type %in% c("Analytical parameter", "Mixture")   ~ "UVCB / Mixture",
      entity_type == "Regulatory entry"                        ~ "Regulatory entry",
      TRUE                                                      ~ "Unidentified"
    ),
    linkability = factor(linkability, levels = c(
      "Structure (InChIKey)",
      "Group \u2014 base compound resolvable",
      "Group \u2014 base compound not resolvable",
      "CAS without structure",
      "UVCB / Mixture",
      "Regulatory entry",
      "Unidentified"
    ))
  )

saveRDS(linkability, here("data", "processed", "linkability_taxonomy.rds"))

df4d <- linkability |>
  count(linkability) |>
  mutate(pct = round(n / sum(n) * 100, 1))

print(df4d)

linkability_colours <- c(
  "Structure (InChIKey)"                      = "#4a90d9",
  "Group \u2014 base compound resolvable"     = "#7bc67e",
  "Group \u2014 base compound not resolvable" = "#b5e5a0",
  "CAS without structure"                     = "#aaaaaa",
  "UVCB / Mixture"                            = "#9b59b6",
  "Regulatory entry"                          = "#e05c5c",
  "Unidentified"                              = "#dddddd"
)

p4d <- ggplot(df4d,
              aes(x = reorder(linkability, n), y = n, fill = linkability)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(n, " (", pct, "%)")),
            hjust = -0.05, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)), labels = comma) +
  scale_fill_manual(values = linkability_colours) +
  coord_flip() +
  labs(
    title    = "Linkability taxonomy of ECHA substances",
    subtitle = "What proportion can be unambiguously linked to a chemical structure?",
    x        = NULL,
    y        = "Number of unique substance names"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p4d)
ggsave(p4d,
       filename = here("output", "figures",
                       "Analysis_4d_Linkable-taxonomie\u2014general_overview_of_the_complete_ECHA-dataset.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns (Analysis_4d_linkability_taxonomy.csv):
#   linkability — one of seven linkability tiers (ordered from most to least
#                 structurally accessible):
#                 "Structure (InChIKey)"                      — direct InChIKey present
#                 "Group — base compound resolvable"          — substance group with PubChem-resolved parent
#                 "Group — base compound not resolvable"      — substance group, parent lookup failed
#                 "CAS without structure"                     — CAS number present but no InChIKey
#                 "UVCB / Mixture"                            — analytical parameter or complex mixture
#                 "Regulatory entry"                          — administrative cross-reference or annotation
#                 "Unidentified"                              — no identifier of any kind
#   n           — number of unique substance names assigned to this tier
#   pct         — share of unique substance names: n / total; sums to 100
write_csv(df4d |> mutate(linkability = as.character(linkability)),
          here("output", "tables", "Analysis_4d_linkability_taxonomy.csv"))

# ------------------------------------------------------------------------------
# Save enriched dataset for downstream analyses (08_embedding_clustering.R,
# 09_embedding_chemont.R) that depend on entity_type
# ------------------------------------------------------------------------------
saveRDS(all_substances,
        here("data", "processed", "all_substances_classified.rds"))

message("04_entity_classification.R: analysis 4 (4a, 4b, 4c, 4d) completed")
