# ==============================================================================
# 10_workload.R
# Analysis 10: Workload estimation — pairwise group-relation matching
#
# PURPOSE
# -------
# The `sommatie_stoffen` source lists substance groups that are defined as
# a sum of members.  Determining the full set of pairwise subset/superset
# relations between these groups (e.g., "are all members of group A also in
# group B?") requires comparing every group against every other group exactly
# once.  This analysis quantifies the resulting workload as a function of the
# number of groups, using a quadratic growth model and an assumed analyst
# throughput rate.
#
# DATA PROVENANCE
# ---------------
# Input: `data/processed/all_substances.rds`
# `n_groups` is the count of distinct substance names in the
# `sommatie_stoffen` source — each unique name represents one substance group
# that must be compared against all others.
#
# METHODOLOGY
# -----------
# Total pairs = n × (n−1) / 2  (combinations without replacement)
# Total days  = ⌈pairs / throughput⌉
# Marginal days per additional group = ⌈(n−1) / throughput⌉
# Own addition: assumed throughput of 100 group-relation assessments per day
# is a professional judgement estimate, not derived from a framework.  It
# should be calibrated against actual analyst performance data if available.
#
# OUTPUTS
# -------
# output/figures/Analysis_10a_workload_curve.pdf
# output/figures/Analysis_10b_workload_marginal.pdf
# output/tables/Analysis_10a_workload_curve.csv
# output/tables/Analysis_10b_workload_marginal.csv
# ==============================================================================

library(dplyr)
library(ggplot2)
library(readr)
library(here)

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 10: Workload model
# ==============================================================================
# INTENT
# Demonstrate that pairwise group-relation assessment scales quadratically
# with the number of groups, making exhaustive manual comparison infeasible
# beyond a relatively small n.  The marginal cost curve (10b) shows how much
# extra work each newly added group imposes, supporting prioritisation
# arguments: it is not cost-neutral to add more groups to scope.
# ==============================================================================

n_groups          <- all_substances |>
  filter(source == "sommatie_stoffen") |>
  distinct(substance_name) |>
  nrow()

# Own addition: 100 relations/day is a professional throughput estimate.
# Calibrate against empirical performance data before citing in reports.
relations_per_day <- 100L

days_required <- function(n, y) ceiling(n * (n - 1L) / 2L / y)
marginal_days <- function(n, y) if (n <= 1L) 0L else ceiling((n - 1L) / y)

n_pairs    <- n_groups * (n_groups - 1L) / 2L
total_days <- days_required(n_groups, relations_per_day)

message(sprintf(
  "Analysis 10: %d groups \u2192 %d unique pairs \u2192 %d person-days (@ %d relations/day)",
  n_groups, n_pairs, total_days, relations_per_day
))

# ------------------------------------------------------------------------------
# 10a: Cumulative days required as a function of n
# ------------------------------------------------------------------------------

p10a_data <- data.frame(
  n    = seq_len(n_groups),
  days = sapply(seq_len(n_groups), days_required, y = relations_per_day)
)

p10a <- ggplot(p10a_data, aes(x = n, y = days)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_vline(xintercept = n_groups,
             linetype = "dashed", colour = "#e05c5c") +
  annotate("label",
           x     = n_groups,
           y     = max(p10a_data$days) * 0.5,
           label = paste0("n = ", n_groups, "\n", total_days, " days"),
           colour = "#e05c5c", size = 3.5, label.size = 0.3) +
  labs(
    title    = paste("Days required at",
                     relations_per_day, "relations per day"),
    subtitle = "Quadratic growth: n\u00d7(n\u22121)/2 pairs",
    x        = "Number of groups (n)",
    y        = "Days required"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p10a)
ggsave(p10a,
       filename = here("output", "figures", "Analysis_10a_workload_curve.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns:
#   n    — number of substance groups (1 … n_groups)
#   days — cumulative person-days required to assess all n×(n−1)/2 pairs at
#           `relations_per_day` throughput; 0 for n = 1 (no pairs possible)
write_csv(p10a_data,
          here("output", "tables", "Analysis_10a_workload_curve.csv"))

# ------------------------------------------------------------------------------
# 10b: Marginal days per additional group
# ------------------------------------------------------------------------------

p10b_data <- data.frame(
  n         = 2L:n_groups,
  marg_days = sapply(2L:n_groups, marginal_days, y = relations_per_day)
)

p10b <- ggplot(p10b_data, aes(x = n, y = marg_days)) +
  geom_line(colour = "firebrick", linewidth = 1) +
  labs(
    title    = paste("Marginal days per additional group at",
                     relations_per_day, "relations/day"),
    subtitle = paste0(
      "Adding 1 group to the existing ", n_groups,
      " groups costs ", marginal_days(n_groups, relations_per_day),
      " person-day(s)"
    ),
    x = "Number of groups (n)",
    y = "Marginal days"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p10b)
ggsave(p10b,
       filename = here("output", "figures",
                       "Analysis_10b_workload_marginal.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns:
#   n         — number of substance groups (2 … n_groups); starts at 2 because
#               a single group has no pairs and therefore zero marginal cost
#   marg_days — additional person-days incurred by adding the n-th group to an
#               existing set of (n−1) groups: ⌈(n−1) / relations_per_day⌉
write_csv(p10b_data,
          here("output", "tables", "Analysis_10b_workload_marginal.csv"))

message("10_workload.R: analysis 10 completed")
