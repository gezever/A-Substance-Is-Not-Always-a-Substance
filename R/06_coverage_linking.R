# ==============================================================================
# 06_coverage_linking.R
# Analysis 6: Coverage of structure-based linking per source
#
# PURPOSE
# -------
# Even within the set of entries that have an InChIKey, the ratio of unique
# structures to total records varies across regulatory sources.  A source
# where many records map to the same InChIKey exhibits structural redundancy
# (the same substance appears multiple times under different names or
# conditions).  This analysis visualises that relationship on a log–log scale,
# placing each source relative to the 1:1 line where every record is a
# distinct structure.
#
# DATA PROVENANCE
# ---------------
# Input: `data/processed/all_substances.rds`
# Records without an InChIKey are ignored in the unique structure count (they
# cannot be linked), but contribute to the total record count.  Sources where
# all records lack InChIKeys therefore appear on the x-axis at (total, 0).
#
# METHODOLOGY
# -----------
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **Cleveland (1993)** *Visualizing Data* (Hobart Press, Summit NJ) | Principles for log-transforming skewed, right-tailed distributions so that differences at all scales are equally readable; justifies the log₁₀ axis transformation used here. |
# | **Tufte (2001)** *The Visual Display of Quantitative Information*, 2nd ed. (Graphics Press) | Reference-line principle: adding a diagonal at slope = 1 gives the reader an immediate reference frame without adding ink that carries no data, making deviations from the 1:1 ideal directly interpretable. |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | 1:1 dashed diagonal as reference | Not prescribed by any framework for this specific context; chosen because the theoretical maximum (every record is a distinct structure) defines the upper bound and makes duplication immediately visible as vertical distance below the line. |
# | Label every point directly | With ~10–15 sources, individual labels are readable; a legend would require additional eye movement without adding information. |
# | Colour and alpha fixed (not mapped to a variable) | Source identity is already encoded by position and text label; adding a colour dimension would suggest a third variable that does not exist here. |
#
# INTERPRETATION
# --------------
# The plot positions each regulatory source in a two-dimensional space defined
# by list size (x) and structural coverage (y).  Four reading patterns matter:
#
# 1. **On the diagonal (linking_ratio ≈ 1.0)**
#    Every record in this source is a chemically distinct structure.  The list
#    is essentially a set of unique molecules with minimal administrative
#    duplication.  Example pattern: a curated restriction list where each entry
#    is one specific substance.
#
# 2. **Below the diagonal (linking_ratio < 1.0)**
#    The source contains more records than unique structures.  Causes include:
#    - The same molecule is listed multiple times under different names or CAS
#      numbers (administrative duplication).
#    - Substance group entries: one InChIKey represents many nominal records.
#    - A large fraction of records lack an InChIKey entirely (non-structured
#      entries pull n_records up without raising n_unique_inchikeys).
#    A low linking_ratio does not imply poor data quality; it reflects the
#    nature of the regulatory instrument (e.g., a sum-of-members list will
#    legitimately have many records per InChIKey).
#
# 3. **Small x, high y/x ratio (top-left region)**
#    Small but structurally dense lists: few records, but nearly all are
#    structure-defined and distinct.  These are the most amenable to automated
#    cross-list matching.
#
# 4. **Large x, low y/x ratio (bottom-right region)**
#    Large lists with poor structural coverage.  A significant share of entries
#    cannot be matched to other sources via InChIKey.  Workload estimates for
#    these sources (see Analysis 10) must account for the non-linkable fraction
#    separately.
#
# **Actionable reading:** sources far below the diagonal with large x-values
# represent the highest workload for manual harmonisation; sources on or near
# the diagonal are candidates for automated cross-source deduplication.
#
# OUTPUTS
# -------
# output/figures/Analysis_6_Coverage_of_structure-based_linking.pdf
# output/tables/Analysis_6_coverage_linking.csv
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(scales)

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 6: Structure coverage scatter (records vs unique InChIKeys)
# ==============================================================================
# INTENT
# Sources below the 1:1 diagonal have more records than unique structures:
# the same molecule appears multiple times (e.g., listed under different names
# or in multiple substance group memberships).  The distance from the diagonal
# quantifies the degree of record duplication, which affects workload estimates
# and de-duplication strategies in later analyses.
# ==============================================================================

df6 <- all_substances |>
  group_by(source) |>
  summarise(
    n_records          = n(),
    n_unique_inchikeys = n_distinct(inchikey[!is.na(inchikey)]),
    .groups = "drop"
  ) |>
  mutate(
    # Own addition: linking_ratio = unique structures / total records.
    # A value of 1.0 means every record is a distinct structure (no duplication).
    # A value near 0 means many records share the same structure, or the source
    # has few structure-defined entries.  Avoids division by zero for empty sources.
    linking_ratio = if_else(n_records > 0L,
                            n_unique_inchikeys / n_records,
                            NA_real_)
  )

p6 <- ggplot(df6, aes(x = n_records, y = n_unique_inchikeys, label = source)) +
  # Cleveland (1993): log axes handle right-skewed distributions spanning orders of magnitude
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60") +
  geom_point(colour = "#4a90d9", size = 3, alpha = 0.8) +
  geom_text(vjust = -0.7, hjust = 0.5, size = 3, colour = "grey30") +
  scale_x_continuous(labels = comma, trans = "log10") +
  scale_y_continuous(labels = comma, trans = "log10") +
  labs(
    title    = "Coverage of structure-based linking per source",
    subtitle = "Points below the dashed diagonal: more records than unique structures (duplication / ambiguity)",
    x        = "Number of records (log\u2081\u2080)",
    y        = "Number of unique InChIKeys (log\u2081\u2080)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p6)
ggsave(p6,
       filename = here("output", "figures",
                       "Analysis_6_Coverage_of_structure-based_linking.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns (Analysis_6_coverage_linking.csv):
#   source             — regulatory source identifier
#   n_records          — total number of rows in all_substances for this source,
#                        including entries without an InChIKey
#   n_unique_inchikeys — number of distinct InChIKeys in this source; entries
#                        without an InChIKey are excluded from this count
#   linking_ratio      — n_unique_inchikeys / n_records; 1.0 = no duplication,
#                        < 1.0 = structural redundancy or non-linkable entries;
#                        NA when n_records = 0 (empty source)
write_csv(df6,
          here("output", "tables", "Analysis_6_coverage_linking.csv"))

message("06_coverage_linking.R: analysis 6 completed")
