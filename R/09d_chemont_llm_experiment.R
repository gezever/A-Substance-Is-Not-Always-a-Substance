# ==============================================================================
# 09d_chemont_llm_experiment.R
# ChemOnt classification via LLM — same validation setup as 09c
#
# PURPOSE
# -------
# Tests whether an LLM (Claude via Anthropic API) can classify chemical
# substances into the ChemOnt taxonomy from substance names alone — the same
# task as 09c but using LLM reasoning instead of embedding cosine similarity.
#
# ChemOnt labels (ground-truth levels only, no definitions) are passed as a
# cached system prompt. The LLM classifies each substance at 5 levels in one
# call: kingdom / superclass / class / subclass / direct_parent.
#
# METHODOLOGY
# -----------
# Ground truth and metrics are identical to 09c (tab:chemont_validation):
#   - Ground truth: structure-defined substances with ClassyFire ChemOnt class
#   - Metrics: Hit@k (k=1,3,5,10) and MRR against chemont_label (direct parent)
# The LLM's 5 outputs form an implicit ranked list (direct_parent first,
# kingdom last), making Hit@k directly comparable to embedding models.
# Results include all 4 embedding baselines from 09c.
#
# COST ESTIMATE (claude-sonnet-4-6, PILOT_N = 500)
# -----------------------------------------------
# System prompt: ~8K tokens (labels only, cached after first call)
# Cache write:   8K × $3/Mtok  ≈ $0.02
# Cache reads:  500 × 8K × $0.30/Mtok ≈ $1.20
# Output:       500 × 300 tokens × $15/Mtok ≈ $2.25
# Total:        ≈ $3.50
#
# DATA PROVENANCE
# ---------------
# Input: data/processed/chemont_ground_truth.csv (same as 09c)
# Cache: data/cache/llm_chemont/ (one JSON file per substance, keyed by SHA1)
#
# OUTPUTS
# -------
# output/tables/Analysis_9d_llm_experiment.csv
# output/figures/Analysis_9d_llm_experiment.pdf
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(purrr)
library(readr)
library(tidyr)
library(stringr)
library(patchwork)
library(jsonlite)
library(reticulate)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

PILOT_N       <- 500L    # substances to evaluate; set NULL for full validation set
PILOT_PRINT   <- FALSE   # set TRUE to print prompts/responses (useful for 20-item debug runs)
CACHE_DIR     <- here("data", "cache", "llm_chemont")
HIERARCHY_TXT <- here("data", "processed", "chemont_hierarchy.txt")

dir.create(CACHE_DIR,                          showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "tables"),           showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "figures"),          showWarnings = FALSE, recursive = TRUE)

if (!nzchar(Sys.getenv("ANTHROPIC_API_KEY"))) {
  stop("ANTHROPIC_API_KEY environment variable not set.\n",
       "Set it with: Sys.setenv(ANTHROPIC_API_KEY = 'your-key') or via .Renviron")
}

# ------------------------------------------------------------------------------
# Step 1: Build ChemOnt hierarchy text for the LLM system prompt
# Uses label + definition from ChemOnt_2_1.rds (no SKOS parsing needed).
# Levels (kingdom/superclass) are inferred from the ground truth.
# ------------------------------------------------------------------------------

build_hierarchy_text <- function(ground_truth) {
  # Labels-only, only values that appear in ground truth.
  # Keeps the system prompt under ~8K tokens to stay within the 30K/min rate limit.
  kingdoms       <- sort(unique(na.omit(ground_truth$kingdom)))
  superclass     <- sort(unique(na.omit(ground_truth$superclass)))
  classes        <- sort(unique(na.omit(ground_truth$class)))
  subclasses     <- sort(unique(na.omit(ground_truth$subclass)))
  # chemont_label is the direct parent (skos:broader); often equals subclass but can be more specific
  direct_parents <- sort(unique(na.omit(ground_truth$chemont_label)))
  # Merge subclass + direct parents — deduplicated set covers both evaluation levels
  sub_and_direct <- sort(unique(c(subclasses, direct_parents)))

  message(sprintf(
    "  Hierarchy labels: %d kingdoms, %d superclasses, %d classes, %d subclasses+direct_parents",
    length(kingdoms), length(superclass), length(classes), length(sub_and_direct)
  ))

  c(
    "=== ChemOnt 2.1 Chemical Taxonomy ===",
    "",
    "Classify substances using ONLY the exact labels listed below.",
    "Taxonomy levels (broadest to most specific):",
    "  kingdom → superclass → class → subclass → direct_parent",
    "Use null when uncertain rather than guessing a wrong label.",
    "",
    "--- KINGDOM ---",
    paste0("  ", kingdoms),
    "",
    "--- SUPERCLASS ---",
    paste0("  ", superclass),
    "",
    "--- CLASS ---",
    paste0("  ", classes),
    "",
    "--- SUBCLASS / DIRECT PARENT ---",
    "  (direct_parent is the most specific level; often equals subclass)",
    paste0("  ", sub_and_direct)
  ) |> paste(collapse = "\n")
}

# Ground truth needed to build hierarchy label list
system("bash bash/chemont_broader_substance.sh")
ground_truth_file <- here("data", "processed", "chemont_ground_truth.csv")
if (!file.exists(ground_truth_file)) {
  stop("Ground truth not found. Run bash/chemont_broader_substance.sh first.")
}
ground_truth <- read_csv(ground_truth_file, show_col_types = FALSE) |>
  filter(!is.na(inchikey), !is.na(chemont_label))

if (!file.exists(HIERARCHY_TXT)) {
  message("Building ChemOnt hierarchy text ...")
  hierarchy_text <- build_hierarchy_text(ground_truth)
  writeLines(hierarchy_text, HIERARCHY_TXT)
  n_tokens_est <- nchar(hierarchy_text) %/% 4L
  message(sprintf("  Written: %s chars (~%s tokens)",
                  format(nchar(hierarchy_text), big.mark = ","),
                  format(n_tokens_est, big.mark = ",")))
  if (n_tokens_est > 25000L)
    warning("Hierarchy text may exceed the 30K token/min rate limit. Consider reducing labels.")
} else {
  message("Hierarchy text already exists: ", HIERARCHY_TXT)
}

# ------------------------------------------------------------------------------
# Step 2: Build validation set (identical to 09c)
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances_classified.rds"))

validation_set <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(substance_name, inchikey) |>
  inner_join(
    ground_truth |> select(-any_of("substance_name")),
    by = "inchikey",
    relationship = "many-to-many"
  ) |>
  filter(!is.na(substance_name), nzchar(substance_name))

message(sprintf("Validation set: %d substance-class pairs", nrow(validation_set)))

if (!is.null(PILOT_N) && nrow(validation_set) > PILOT_N) {
  set.seed(42L)
  validation_set <- validation_set |> slice_sample(n = PILOT_N)
  message(sprintf("Pilot mode: sampled %d substances", PILOT_N))
}

# ------------------------------------------------------------------------------
# Step 3: Call LLM classifier (Python script)
# ------------------------------------------------------------------------------

py_script       <- here("python", "llm_chemont_classify.py")
python          <- py_exe()
substances_json <- tempfile(fileext = ".json")
output_json     <- tempfile(fileext = ".json")

write_json(unique(validation_set$substance_name), substances_json, auto_unbox = FALSE)

pilot_flag <- if (isTRUE(PILOT_PRINT)) "--pilot" else character(0)

message("Calling LLM classifier ...")
ret <- system2(
  python,
  args   = c(py_script, substances_json, HIERARCHY_TXT, CACHE_DIR, output_json, pilot_flag),
  stdout = "",
  stderr = ""
)
if (ret != 0L) stop("Python LLM classifier failed (exit code ", ret, ")")

llm_raw <- fromJSON(output_json, simplifyDataFrame = TRUE) |>
  as_tibble() |>
  rename(class = class_)   # match ground truth column name

file.remove(substances_json, output_json)

message(sprintf("LLM results: %d rows", nrow(llm_raw)))

# Token usage summary
if ("usage" %in% names(llm_raw) && is.data.frame(llm_raw$usage)) {
  u <- llm_raw$usage
  message(sprintf(
    "Tokens — input: %s | cache_creation: %s | cache_read: %s | output: %s",
    format(sum(u$input_tokens,                   na.rm = TRUE), big.mark = ","),
    format(sum(u$cache_creation_input_tokens,     na.rm = TRUE), big.mark = ","),
    format(sum(u$cache_read_input_tokens,         na.rm = TRUE), big.mark = ","),
    format(sum(u$output_tokens,                   na.rm = TRUE), big.mark = ",")
  ))
}

# ------------------------------------------------------------------------------
# Step 4: Evaluate — Hit@k and MRR comparable to 09c table
#
# The LLM outputs 5 ordered levels (direct_parent → subclass → class →
# superclass → kingdom), treated as an implicit ranked list of up to 5
# candidates for the correct chemont_label (= direct parent in ground truth).
# This makes Hit@1–5 and MRR directly comparable to the embedding results.
# Hit@10 is set equal to Hit@5 (max rank for LLM is 5).
#
# Per-level Hit@1 (kingdom/superclass/class/subclass) are also computed for
# interpretability — they are not in the 09c table but complement it.
# ------------------------------------------------------------------------------

eval_data <- validation_set |>
  left_join(
    llm_raw |> select(substance_name,
                      llm_kingdom       = kingdom,
                      llm_superclass    = superclass,
                      llm_class         = class,
                      llm_subclass      = subclass,
                      llm_direct_parent = direct_parent,
                      confidence),
    by = "substance_name"
  )

# Rank of chemont_label in LLM's implicit 5-level ranking
llm_ranks <- map_int(seq_len(nrow(eval_data)), function(i) {
  truth      <- eval_data$chemont_label[i]
  candidates <- c(eval_data$llm_direct_parent[i],
                  eval_data$llm_subclass[i],
                  eval_data$llm_class[i],
                  eval_data$llm_superclass[i],
                  eval_data$llm_kingdom[i])
  idx <- which(!is.na(candidates) & candidates == truth)
  if (length(idx) == 0L) NA_integer_ else min(idx)
})

hit1    <- function(pred, truth) mean(!is.na(pred) & pred == truth, na.rm = TRUE)
prop_se <- function(p, n) sqrt(p * (1 - p) / n)   # SE for a proportion
z95     <- 1.96                                     # 95% CI multiplier

n_total          <- nrow(eval_data)
reciprocal_ranks <- ifelse(is.na(llm_ranks), 0, 1 / llm_ranks)

metrics_llm <- tibble(
  model_id     = "Claude Sonnet 4.6",
  # Strict Hit@k against chemont_label — NA counts as miss (rank = Inf)
  hit_at_1     = mean(!is.na(llm_ranks) & llm_ranks == 1L),
  hit_at_3     = mean(!is.na(llm_ranks) & llm_ranks <= 3L),
  hit_at_5     = mean(!is.na(llm_ranks) & llm_ranks <= 5L),
  hit_at_10    = mean(!is.na(llm_ranks) & llm_ranks <= 5L),  # max rank = 5
  mrr          = mean(reciprocal_ranks),
  # Per-level Hit@1 — comparable to 09c per-level columns
  hit_at_1_sub = hit1(eval_data$llm_subclass,   eval_data$subclass),
  hit_at_1_cls = hit1(eval_data$llm_class,      eval_data$class),
  hit_at_1_sup = hit1(eval_data$llm_superclass, eval_data$superclass),
  hit_at_1_kin = hit1(eval_data$llm_kingdom,    eval_data$kingdom),
  n_eval       = n_total
) |>
  mutate(
    # 95% CI (normal approximation) for all proportion metrics
    across(
      c(hit_at_1, hit_at_3, hit_at_5, hit_at_10,
        hit_at_1_sub, hit_at_1_cls, hit_at_1_sup, hit_at_1_kin),
      list(
        lo = \(p) pmax(0, p - z95 * prop_se(p, n_eval)),
        hi = \(p) pmin(1, p + z95 * prop_se(p, n_eval))
      ),
      .names = "{.col}_{.fn}"
    ),
    # 95% CI for MRR (SE of a mean, not a proportion)
    mrr_lo = pmax(0, mrr - z95 * sd(reciprocal_ranks) / sqrt(n_eval)),
    mrr_hi = pmin(1, mrr + z95 * sd(reciprocal_ranks) / sqrt(n_eval))
  )

print(metrics_llm)

# Append all embedding baselines from 09c — produces a complete supplementary table
# Column order matches tab:chemont_validation in 30_results.tex
baseline_file <- here("output", "tables", "Analysis_9c_model_validation.csv")
metrics_combined <- metrics_llm
if (file.exists(baseline_file)) {
  baseline <- read_csv(baseline_file, show_col_types = FALSE) |>
    select(model_id, hit_at_1, hit_at_3, hit_at_5, hit_at_10, mrr,
           hit_at_1_sub, hit_at_1_cls, hit_at_1_sup, hit_at_1_kin, n_eval)
  metrics_combined <- bind_rows(baseline, metrics_llm)
  message(sprintf("Appended %d embedding baselines from 09c", nrow(baseline)))
}

write_csv(metrics_combined, here("output", "tables", "Analysis_9d_llm_experiment.csv"))

# ------------------------------------------------------------------------------
# Plot: two panels
# Panel 1 — strict Hit@k: LLM vs all embedding models (comparable to 09c table)
# Panel 2 — per-level Hit@1: LLM only (hierarchy breakdown, unique to 09d)
# ------------------------------------------------------------------------------

# Panel 1: strict Hit@k
hitatk_long <- metrics_combined |>
  select(model_id, `Hit@1` = hit_at_1, `Hit@3` = hit_at_3,
         `Hit@5` = hit_at_5, `Hit@10` = hit_at_10) |>
  pivot_longer(-model_id, names_to = "k", values_to = "value") |>
  mutate(k = factor(k, levels = c("Hit@1", "Hit@3", "Hit@5", "Hit@10")))

p1 <- ggplot(hitatk_long, aes(x = k, y = value, fill = model_id)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", value)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3) +
  scale_y_continuous(limits = c(0, 1.1), expand = expansion(mult = c(0, 0))) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Strict Hit@k (direct ChemOnt parent, chemont_label)",
    subtitle = "LLM: 5-level implicit ranking; embedding: ranked over 4,825 classes",
    x = NULL, y = "Fraction correct"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

# Panel 2: per-level Hit@1 for LLM only
levels_long <- tibble(
  level    = factor(c("subclass", "class", "superclass", "kingdom"),
                    levels = c("subclass", "class", "superclass", "kingdom")),
  hit_at_1 = c(metrics_llm$hit_at_1_sub, metrics_llm$hit_at_1_cls,
               metrics_llm$hit_at_1_sup, metrics_llm$hit_at_1_kin),
  lo       = c(metrics_llm$hit_at_1_sub_lo, metrics_llm$hit_at_1_cls_lo,
               metrics_llm$hit_at_1_sup_lo, metrics_llm$hit_at_1_kin_lo),
  hi       = c(metrics_llm$hit_at_1_sub_hi, metrics_llm$hit_at_1_cls_hi,
               metrics_llm$hit_at_1_sup_hi, metrics_llm$hit_at_1_kin_hi)
)

p2 <- ggplot(levels_long, aes(x = level, y = hit_at_1)) +
  geom_col(fill = "#66C2A5", width = 0.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.3f", hit_at_1)), vjust = -0.4, size = 3.5) +
  scale_y_continuous(limits = c(0, 1.1), expand = expansion(mult = c(0, 0))) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "LLM Hit@1 per hierarchy level (most → least specific)",
    subtitle = sprintf("n = %d; higher levels (kingdom) should be easier than subclass",
                       nrow(eval_data)),
    x = NULL, y = "Hit@1"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

p_final <- p1 / p2 + plot_layout(heights = c(2, 1))

ggsave(p_final,
       filename = here("output", "figures", "Analysis_9d_llm_experiment.pdf"),
       device = "pdf", height = 25, width = 30, units = "cm")

message("09d_chemont_llm_experiment.R: completed")
