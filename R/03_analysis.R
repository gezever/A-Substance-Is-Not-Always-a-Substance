# ==============================================================================
# 03_analysis.R
# Analyses 1–3: Substance identity and identifier consistency
#
# PURPOSE
# -------
# This script answers three foundational questions about the ECHA substance
# dataset before any further analysis:
#   1. What share of regulatory entries are chemically identifiable
#      (structure-defined) versus non-structure entities?
#   2. Does interoperability — the ability to link entries across regulatory
#      lists via InChIKey — differ systematically between regulatory sources?
#   3. How consistent is the mapping between CAS numbers and InChIKeys across
#      the dataset?
# These questions establish the epistemic limits of all subsequent analyses:
# any scoring or linking approach is constrained by the ~21 % InChIKey coverage.
#
# DATA PROVENANCE
# ---------------
# Input: `data/processed/all_substances.rds`
# Produced by `02_import.R` from ECHA open-data downloads combined with
# CAS/InChIKey resolution.  One row per (substance_name, source) combination;
# a substance may appear across multiple regulatory sources.
# An absent `inchikey` means the entry has no resolvable chemical structure.
#
# METHODOLOGY
# -----------
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **Cleveland (1993)** *Visualizing Data* (Hobart Press, Summit NJ) | Log₁₀ y-axis transformation in analyses 3a and 3b: counts of CAS/InChIKey mappings are strongly right-skewed; log scaling makes the long tail visible without compressing the dominant mode. |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | Binary split (InChIKey present / absent) for Analysis 1 | The InChIKey is the only unambiguous structural identifier in the dataset; the split directly measures what fraction is accessible to structure-based analyses. No finer taxonomy is needed at this stage (see Analysis 4 for the full entity-type taxonomy). |
# | Per-source stacked bar sorted by `pct_with` (Analysis 2) | Sorting by descending InChIKey coverage makes the gradient from structure-rich to structure-poor instruments immediately visible; alphabetical ordering would hide this pattern. |
# | Separate histograms for CAS→InChIKey and InChIKey→CAS (Analyses 3a/3b) | The two inconsistency directions have different causes and different remediation strategies; combining them into one plot would obscure which direction drives which substances. |
#
# INTERPRETATION
# --------------
# **Analysis 1 — Structure-defined vs non-structure entities**
# The bar heights establish the baseline constraint for all downstream analyses.
# The non-structure fraction (~79 %) is not a data quality problem — it reflects
# how regulations define "substance" (groups, mixtures, parameters).  Any
# approach that claims to cover "all regulatory substances" must address this
# fraction explicitly.
#
# **Analysis 2 — Interoperability by source**
# A strong gradient (some sources near 100 %, others near 0 %) means that
# cross-list linking is feasible for some instruments but structurally
# impossible for others.  Sources with low `pct_with` cannot meaningfully
# participate in InChIKey-based overlap analyses (7, 12); they appear in
# workload analyses (10) instead.
#
# **Analysis 3a — CAS → multiple InChIKeys**
# Most CAS numbers map to exactly one InChIKey (mode at x = 1).  CAS numbers
# with n_inchikeys > 1 are ambiguous: the same registry number has been
# associated with structurally different entries across sources.  These are
# candidates for manual review before using CAS as a linking key.
#
# **Analysis 3b — InChIKey → multiple CAS numbers**
# The same structure legitimately carries multiple CAS numbers (historical
# assignments, stereoisomers grouped under one entry).  A high n_cas is not
# necessarily an error, but it means CAS-based deduplication will
# under-collapse the dataset relative to InChIKey-based deduplication.
#
# **Combined reading for 3a + 3b:** motivates using InChIKey (not CAS) as the
# primary linking key in all subsequent cross-list analyses.
#
# OUTPUTS
# -------
# output/figures/Analysis_1_What_is_a_substance.pdf
# output/figures/Analysis_2_interoperability_differs_by_regulation.pdf
# output/figures/Analysis_3a_CAS↔InChIKey_inconsistencies_multiple_INCHIKEYS.pdf
# output/figures/Analysis_3b_CAS↔InChIKey_inconsistencies_multiple_CAS_numbers.pdf
# output/tables/Analysis_1_structure_vs_nonstructure.csv
# output/tables/Analysis_2_interoperability_per_source.csv
# output/tables/Analysis_3a_cas_to_inchikey.csv
# output/tables/Analysis_3b_inchikey_to_cas.csv
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(readr)
library(scales)

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 1: What is a "substance"?
# Structure-defined vs non-structure entities
# ==============================================================================
# INTENT
# Not all regulatory entries refer to chemically identifiable substances.
# This figure quantifies how many entries carry an InChIKey (and are therefore
# linkable via structure) versus those that do not.  The split determines what
# fraction of the dataset is accessible to cheminformatics-based approaches in
# later analyses.
# ==============================================================================

df1 <- all_substances |>
  mutate(has_inchikey = !is.na(inchikey)) |>
  count(has_inchikey) |>
  mutate(
    label = if_else(has_inchikey,
                    "Structure-defined\n(InChIKey present)",
                    "Non-structure entity\n(InChIKey absent)"),
    pct   = n / sum(n)
  )

p1 <- ggplot(df1, aes(x = label, y = n, fill = has_inchikey)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(n, "\n(", percent(pct, accuracy = 0.1), ")")),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("FALSE" = "#e05c5c", "TRUE" = "#4a90d9")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Analysis 1: What is a 'substance'?",
    subtitle = "A large proportion of regulatory entries are not chemically identifiable substances",
    x        = NULL,
    y        = "Number of records"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p1)
ggsave(p1,
       filename = here("output", "figures", "Analysis_1_What_is_a_substance.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns (Analysis_1_structure_vs_nonstructure.csv):
#   has_inchikey — TRUE = record has a resolvable InChIKey (structure-defined);
#                  FALSE = record has no InChIKey (non-structure entity)
#   label        — human-readable category label used in the figure
#   n            — number of records in this category; the two rows sum to the
#                  total number of records in all_substances
#   pct          — share of total records: n / sum(n); the two rows sum to 1.0
write_csv(df1, here("output", "tables", "Analysis_1_structure_vs_nonstructure.csv"))

# ==============================================================================
# Analysis 2: Per source — interoperability differs by regulation
# ==============================================================================
# INTENT
# If InChIKey coverage were uniform across regulatory lists, the choice of
# list would not affect the feasibility of cross-list linking.  This figure
# tests that assumption: it shows whether chemical identifiability varies
# systematically with the regulatory instrument.  A strong gradient here
# implies that interoperability is not a data quality issue but a structural
# consequence of how different regulations define "substance".
# ==============================================================================

df2_wide <- all_substances |>
  group_by(source) |>
  summarise(
    with_inchikey    = sum(!is.na(inchikey)),
    without_inchikey = sum(is.na(inchikey)),
    total            = n(),
    .groups = "drop"
  ) |>
  mutate(pct_with = with_inchikey / total) |>
  arrange(desc(pct_with))

df2 <- df2_wide |>
  pivot_longer(c(with_inchikey, without_inchikey),
               names_to = "type", values_to = "n") |>
  mutate(
    type   = recode(type,
                    with_inchikey    = "InChIKey present",
                    without_inchikey = "InChIKey absent"),
    source = factor(source, levels = unique(source[type == "InChIKey present"]))
  )

p2 <- ggplot(df2, aes(x = source, y = n, fill = type)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c("InChIKey present" = "#4a90d9",
                               "InChIKey absent"  = "#e05c5c")) +
  coord_flip() +
  labs(
    title    = "Interoperability by source",
    subtitle = "Regulation determines whether substances are chemically identifiable",
    x        = NULL,
    y        = "Share of records",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    plot.subtitle    = element_text(colour = "grey40")
  )

print(p2)
ggsave(p2,
       filename = here("output", "figures",
                       "Analysis_2_interoperability_differs_by_regulation.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns (Analysis_2_interoperability_per_source.csv):
#   source           — regulatory source identifier
#   with_inchikey    — number of records in this source that have an InChIKey
#   without_inchikey — number of records in this source without an InChIKey
#   total            — with_inchikey + without_inchikey; total records in source
#   pct_with         — with_inchikey / total; share of structure-defined entries;
#                      1.0 = fully structure-defined, 0.0 = no linkable entries
write_csv(df2_wide, here("output", "tables", "Analysis_2_interoperability_per_source.csv"))

# ==============================================================================
# Analysis 3: CAS ↔ InChIKey inconsistencies
# ==============================================================================
# INTENT
# CAS numbers are widely used as substance identifiers in regulatory databases,
# but they were not designed to be structurally unique: a CAS may be assigned
# to a substance class rather than a single structure, and the same structure
# may have received multiple CAS numbers over time.  This analysis quantifies
# both directions of inconsistency (CAS → multiple InChIKeys; InChIKey →
# multiple CAS numbers) to motivate the use of InChIKey as the primary linking
# key in subsequent cross-list analyses.
# ==============================================================================

# 3a: One CAS → multiple InChIKeys
# Own addition: log₁₀ y-axis to expose the long tail of CAS numbers with many
# mapped structures, which would be invisible on a linear scale given the
# strong mode at n_inchikeys = 1.
cas_to_inchi <- all_substances |>
  filter(!is.na(cas_number), !is.na(inchikey)) |>
  group_by(cas_number) |>
  summarise(n_inchikeys = n_distinct(inchikey), .groups = "drop")

p3a <- ggplot(cas_to_inchi, aes(x = n_inchikeys)) +
  geom_histogram(binwidth = 1, fill = "#4a90d9", colour = "white") +
  scale_y_log10(labels = comma) +
  labs(
    title    = "CAS \u2192 InChIKey (one CAS, how many InChIKeys?)",
    subtitle = "A CAS number should uniquely refer to a single structure",
    x        = "Number of distinct InChIKeys per CAS",
    y        = "Number of CAS numbers (log\u2081\u2080)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p3a)
ggsave(p3a,
       filename = here("output", "figures",
                       "Analysis_3a_CAS\u2194InChIKey_inconsistencies_multiple_INCHIKEYS.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns (Analysis_3a_cas_to_inchikey.csv):
#   cas_number  — CAS registry number (only rows where both CAS and InChIKey are
#                 present in all_substances are included)
#   n_inchikeys — number of distinct InChIKeys observed for this CAS number across
#                 all regulatory sources; 1 = consistent (one structure per CAS),
#                 > 1 = inconsistency (same CAS mapped to multiple structures)
write_csv(cas_to_inchi, here("output", "tables", "Analysis_3a_cas_to_inchikey.csv"))

# 3b: One InChIKey → multiple CAS numbers
inchi_to_cas <- all_substances |>
  filter(!is.na(cas_number), !is.na(inchikey)) |>
  group_by(inchikey) |>
  summarise(n_cas = n_distinct(cas_number), .groups = "drop")

p3b <- ggplot(inchi_to_cas, aes(x = n_cas)) +
  geom_histogram(binwidth = 1, fill = "#e07b3a", colour = "white") +
  scale_y_log10(labels = comma) +
  labs(
    title    = "InChIKey \u2192 CAS (one structure, how many CAS numbers?)",
    subtitle = "The same chemical structure can carry multiple CAS numbers",
    x        = "Number of distinct CAS numbers per InChIKey",
    y        = "Number of InChIKeys (log\u2081\u2080)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p3b)
ggsave(p3b,
       filename = here("output", "figures",
                       "Analysis_3b_CAS\u2194InChIKey_inconsistencies_multiple_CAS_numbers.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns (Analysis_3b_inchikey_to_cas.csv):
#   inchikey — InChIKey of the substance (only rows where both InChIKey and CAS
#              are present in all_substances are included)
#   n_cas    — number of distinct CAS numbers observed for this InChIKey across
#              all regulatory sources; 1 = consistent (one CAS per structure),
#              > 1 = the same structure has been assigned multiple CAS numbers
write_csv(inchi_to_cas, here("output", "tables", "Analysis_3b_inchikey_to_cas.csv"))

message(sprintf(
  "Analysis 3: %d CAS numbers map to >1 InChIKey; %d InChIKeys have >1 CAS",
  sum(cas_to_inchi$n_inchikeys > 1),
  sum(inchi_to_cas$n_cas > 1)
))

message("03_analysis.R: analyses 1\u20133 completed")
