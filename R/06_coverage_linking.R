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
# OUTPUTS
# -------
# output/figures/Analysis_6_Coverage_of_structure-based_linking.pdf
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)
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
# Own addition: log₁₀ axes to handle the large range of list sizes; the
# dashed diagonal is the theoretical upper bound (records = structures).
# ==============================================================================

df6 <- all_substances |>
  group_by(source) |>
  summarise(
    n_records          = n(),
    n_unique_inchikeys = n_distinct(inchikey[!is.na(inchikey)]),
    .groups = "drop"
  )

p6 <- ggplot(df6, aes(x = n_records, y = n_unique_inchikeys, label = source)) +
  geom_point(colour = "#4a90d9", size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60") +
  geom_text(vjust = -0.7, hjust = 0.5, size = 3, colour = "grey30") +
  scale_x_continuous(labels = comma, trans = "log10") +
  scale_y_continuous(labels = comma, trans = "log10") +
  labs(
    title    = "Coverage of structure-based linking per source",
    subtitle = "Points below the diagonal = more records than unique structures (duplication / ambiguity)",
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

message("06_coverage_linking.R: analysis 6 completed")
