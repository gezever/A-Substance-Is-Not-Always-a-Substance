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
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **Combinatorics: C(n,2) = n(n−1)/2** (standard combinatorial identity) | Exact count of unordered pairs from n elements; determines the total number of assessments required without reference to any assumption. |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | Throughput = 100 relations/day | Professional judgement estimate for the rate at which an analyst can assess one pairwise group-relation; not derived from any published benchmark.  Must be calibrated against empirical performance data before citing in policy reports. |
# | Ceiling function (⌈⌉) for day counts | Partial days of work still require a full working day; ceiling avoids under-estimating the resource requirement. |
# | `substance_name` as the unit of count for non-structured sources (10c) | Each distinct name represents one entity that must be compared against all others; duplicate names within a source would over-count work, so `n_distinct()` is used. |
# | Horizontal bar chart (10c) instead of line chart per source | Multiple overlapping curves are unreadable at this scale; a bar chart sorted by descending workload makes the per-list effort directly comparable without visual clutter. |
#
# INTERPRETATION
# --------------
# **Analyses 10a + 10b — sommatie_stoffen workload**
# The curves grow quadratically: doubling the number of groups quadruples the
# total work.  The red dashed line (10a) and the annotated current-n value show
# where the actual dataset sits on this curve.  If the number of groups grows
# (new substances added to scope), the marginal cost (10b) shows exactly how
# many additional person-days each new group adds.
#
# A steep marginal curve (10b rising toward the current n) means that scope
# expansion is not cost-neutral.  This is the key message for stakeholders
# who propose adding more groups: each addition makes the total problem
# substantially harder.
#
# **Analysis 10c — non-structured substances per ECHA list**
# The bar length is the within-list workload if all non-structured entries were
# assessed pairwise.  The dashed aggregate line is the upper-bound workload if
# all lists were pooled into one assessment — it is almost always far higher
# than any individual list because cross-list unique names accumulate rapidly.
#
# Lists with short bars are feasible for manual harmonisation; lists with bars
# approaching the aggregate dashed line dominate the total work and should be
# prioritised for automated or pre-structured approaches.
#
# OUTPUTS
# -------
# output/figures/Analysis_10a_workload_curve.pdf
# output/figures/Analysis_10b_workload_marginal.pdf
# output/figures/Analysis_10c_workload_non_structured.pdf
# output/tables/Analysis_10a_workload_curve.csv
# output/tables/Analysis_10b_workload_marginal.csv
# output/tables/Analysis_10c_workload_non_structured.csv
# ==============================================================================

library(dplyr)
library(ggplot2)
library(readr)
library(scales)
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

# ==============================================================================
# Analysis 10c: Workload curves — non-structured substances per ECHA list
# ==============================================================================
# INTENT
# The sommatie_stoffen source (10a/10b) is one specific pool of substance
# groups.  The broader question is: how large is the pairwise assessment
# workload if non-structured substances must be assessed within each
# individual regulatory list?  "Non-structured" means any entry without an
# InChIKey — substance groups, UVCBs, mixtures, analytical parameters, and
# unclassified entries.
# sommatie_stoffen is excluded because it is already modelled in 10a/10b.
# Each source gets its own curve (based on unique non-structured names in
# that source); an aggregate curve across all sources combined is added as
# a reference.  Per-source curves illustrate list-level feasibility; the
# aggregate curve shows the upper-bound workload if all lists were treated
# as a single pool.
# ==============================================================================

# Per-source n (unique non-structured names per source, excl. sommatie_stoffen)
ns_per_source <- all_substances |>
  filter(is.na(inchikey), source != "sommatie_stoffen") |>
  group_by(source) |>
  summarise(n = n_distinct(substance_name), .groups = "drop") |>
  arrange(desc(n))

# Aggregate n across all included sources (unique names, not sum of per-source)
n_non_structured <- all_substances |>
  filter(is.na(inchikey), source != "sommatie_stoffen") |>
  distinct(substance_name) |>
  nrow()

message(sprintf(
  "Analysis 10c: aggregate %d unique non-structured names across %d sources",
  n_non_structured, nrow(ns_per_source)
))

# Total days per source at the actual n (endpoint of the workload curve)
p10c_data <- ns_per_source |>
  mutate(total_days = days_required(n, relations_per_day))

# Aggregate total days across all included sources combined
total_days_agg <- days_required(n_non_structured, relations_per_day)

p10c <- ggplot(p10c_data,
               aes(x = reorder(source, total_days), y = total_days)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(aes(label = scales::comma(total_days)),
            hjust = -0.1, size = 3.2) +
  geom_hline(yintercept = total_days_agg,
             linetype = "dashed", colour = "#e05c5c", linewidth = 0.8) +
  annotate("label",
           x     = 1,
           y     = total_days_agg,
           label = paste0("Aggregate\n", scales::comma(total_days_agg), " days"),
           colour = "#e05c5c", size = 3, label.size = 0.3, hjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::comma) +
  coord_flip() +
  labs(
    title    = paste("Total days required at", relations_per_day,
                     "relations per day \u2014 non-structured substances per ECHA list"),
    subtitle = "Bar = days for within-list pairwise assessment; dashed line = aggregate if all lists treated as one pool",
    x        = NULL,
    y        = "Person-days required"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p10c)
ggsave(p10c,
       filename = here("output", "figures",
                       "Analysis_10c_workload_non_structured.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# Columns:
#   source     — regulatory source name (sommatie_stoffen excluded)
#   n          — number of unique non-structured substance names in that source
#   total_days — person-days required for within-list pairwise assessment:
#                ⌈n×(n−1)/2 / relations_per_day⌉
write_csv(p10c_data,
          here("output", "tables", "Analysis_10c_workload_non_structured.csv"))

message("10_workload.R: analysis 10 completed")
