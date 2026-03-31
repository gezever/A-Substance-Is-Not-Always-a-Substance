# ==============================================================================
# 11_ambition_fte.R
# Analysis 11: Ambition vs. capacity — FTE required to meet the 2030 target
#
# PURPOSE
# -------
# Analysis 10 established that exhaustive pairwise group-relation assessment
# of the `sommatie_stoffen` groups requires a large number of person-days.
# This analysis translates that total into a staffing requirement (FTE) by
# dividing the work by the working days remaining until a notional 2030
# policy deadline.  The result provides a concrete capacity figure for
# resource planning discussions.
#
# DATA PROVENANCE
# ---------------
# Input: `data/processed/all_substances.rds`
# `n_groups` and `total_days` are recomputed from scratch using the same
# parameters as `10_workload.R` to avoid cross-script state dependencies.
#
# METHODOLOGY
# -----------
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **FTE definition (Eurostat / EU labour statistics)** | One FTE = one person working full-time (100 %) on the task for the full period; fractional FTE = part-time or shared allocation. Used here to express the staffing requirement in a unit that is directly actionable for resource planning. |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | 200 working days/year | Standard European approximation: 52 weeks × 5 days = 260 days, minus ~60 days for leave and public holidays.  Adjust if organisation-specific norms (e.g., shorter leave entitlements) apply. |
# | Deadline = 2030-01-01 | Notional policy target; the actual deadline should be substituted when known.  The FTE figure is sensitive to this choice: a nearer deadline produces a higher FTE requirement. |
# | Ceiling function for FTE | A fractional FTE below 1.0 still requires at least one dedicated analyst; ceiling prevents under-reporting the staffing need. |
# | Throughput = 100 relations/day | Inherited from Analysis 10; same calibration caveat applies. |
#
# INTERPRETATION
# --------------
# The output is a single number: the FTE required to complete all pairwise
# group-relation assessments by the 2030 deadline, given the current group
# count and throughput assumption.
#
# **Sensitivity:** the FTE figure is proportional to n² (through `total_days`)
# and inversely proportional to remaining time.  Both drivers move in an
# unfavourable direction over time: the group list may grow (more groups →
# more pairs), and every day that passes reduces the remaining capacity.
# Re-running this script closer to 2030 will produce a higher FTE estimate.
#
# **Actionable reading:** an FTE > 1 means one analyst cannot complete the
# work alone by the deadline at the assumed throughput.  An FTE of, say, 3
# means the task requires three full-time analysts, or a combination of more
# analysts at lower throughput.  If the estimated FTE is implausibly high,
# the conclusion is that exhaustive pairwise assessment is infeasible and a
# prioritised or automated approach is needed (see Analysis 14).
#
# OUTPUTS
# -------
# Console message with FTE estimate.
# output/tables/Analysis_11_fte_estimate.csv
# ==============================================================================

library(dplyr)
library(readr)
library(here)

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 11: FTE requirement to 2030
# ==============================================================================
# INTENT
# Convert the quadratic workload estimate (analysis 10) into a staffing
# figure that is directly actionable: "to complete all pairwise assessments
# by 2030 with the current group list size, X FTE must be dedicated to this
# task starting today."  The figure is proportional to n² and inversely
# proportional to remaining time, so it will grow rapidly if either the group
# list expands or the deadline is not met.
# ==============================================================================

# Recompute workload parameters (same as 10_workload.R)
n_groups          <- all_substances |>
  filter(source == "sommatie_stoffen") |>
  distinct(substance_name) |>
  nrow()

# Own addition: 100 relations/day throughput assumption (see 10_workload.R)
relations_per_day <- 100L
total_days        <- ceiling(n_groups * (n_groups - 1L) / 2L / relations_per_day)

# Capacity calculation
end_date             <- as.Date("2030-01-01")
start_date           <- Sys.Date()
remaining_years      <- as.numeric(difftime(end_date, start_date,
                                            units = "days")) / 365.25
# Own addition: 200 working days/year is a standard European estimate.
working_days_to_2030 <- ceiling(remaining_years * 200)

fte_required <- ceiling(total_days / working_days_to_2030)

message(sprintf(
  "Analysis 11: %d working days until 2030 | %d days of work | %d FTE required",
  working_days_to_2030, total_days, fte_required
))

# ==============================================================================
# Serialise results
# ==============================================================================
# INTENT
# Export all model inputs and the derived FTE estimate as a single-row CSV so
# that the result is reproducible and can be consumed by reporting tools or
# policy documents without re-running R.  The run_date column records when the
# estimate was produced; because working_days_to_2030 decreases over time, the
# FTE figure will increase if re-run closer to the deadline.
# Columns:
#   run_date             — date on which this estimate was produced (ISO 8601)
#   deadline             — policy deadline used for the calculation (ISO 8601)
#   n_groups             — number of distinct substance groups in sommatie_stoffen
#   n_pairs              — total pairwise comparisons: n×(n−1)/2
#   relations_per_day    — assumed analyst throughput (group-relation assessments/day)
#   total_days           — total person-days required: ⌈n_pairs / relations_per_day⌉
#   working_days_to_2030 — working days remaining until deadline: ⌈remaining_years×200⌉
#   fte_required         — FTE needed to finish by deadline: ⌈total_days / working_days_to_2030⌉
# ==============================================================================

fte_result <- tibble(
  run_date             = as.character(start_date),
  deadline             = as.character(end_date),
  n_groups             = n_groups,
  n_pairs              = n_groups * (n_groups - 1L) / 2L,
  relations_per_day    = relations_per_day,
  total_days           = total_days,
  working_days_to_2030 = working_days_to_2030,
  fte_required         = fte_required
)

write_csv(fte_result,
          here("output", "tables", "Analysis_11_fte_estimate.csv"))

message("11_ambition_fte.R: analysis 11 completed")
