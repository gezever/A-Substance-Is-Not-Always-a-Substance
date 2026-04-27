# ==============================================================================
# 09e_llm_nonstructure.R
# LLM-based ChemOnt classification of non-structure substances
#
# PURPOSE
# -------
# Applies the LLM classifier from 09d to non-structure substances (no InChIKey)
# that cannot be classified by ClassyFire.  The embedding-based approach in
# 09_embedding_chemont.R assigned ChemOnt classes to these substances via cosine
# similarity; this script asks whether an LLM can do the same — and whether
# the two methods agree.
#
# Because there is no structural ground truth for non-structure substances,
# quality is assessed via:
#   (a) the LLM's own confidence signal (validated as reliable in 09d)
#   (b) agreement with the embedding match from 09_embedding_chemont.R
#       where both methods produce a result
#   (c) null rates per hierarchy level (appropriate uncertainty)
#
# MATCHABLE FILTER
# ----------------
# Reaction masses, petroleum fractions, and UVCB substances are excluded before
# passing names to the classifier.  This filter is method-independent: a name
# that does not describe a single chemical entity cannot yield a meaningful
# class assignment regardless of the classifier used.  See Analysis 9d in
# 09_embedding_chemont.R for the canonical definition of this filter.
#
# METHODOLOGY
# -----------
# Identical to 09d: the ChemOnt label list is passed as a cached system prompt;
# the LLM returns kingdom / superclass / class / subclass / direct_parent in a
# single call per substance.  Results are cached per substance (SHA-1 key) to
# avoid repeated API calls.  The same Python classifier and cache format are
# used; only the cache directory differs.
#
# COST ESTIMATE (claude-sonnet-4-6)
# -----------------------------------
# System prompt: ~8K tokens (cached after first call)
# Per PILOT_N = 500 substances:
#   Cache write:  8K × $3/Mtok              ≈ $0.02
#   Cache reads:  500 × 8K × $0.30/Mtok    ≈ $1.20
#   Output:       500 × ~156tok × $15/Mtok  ≈ $1.17
#   Total pilot:                             ≈ $2.40
# Full run (9,547 matchable substances):     ≈ $45
# Set PILOT_N <- NULL to run the full set.
#
# DATA PROVENANCE
# ---------------
# Input 1: data/processed/all_substances_classified.rds
# Input 2: data/processed/chemont_hierarchy.txt  (built by 09d)
# Input 3: output/tables/final_matches.csv       (embedding matches from 09)
# Cache:   data/cache/llm_chemont_nonstructure/  (separate from 09d cache)
#
# OUTPUTS
# -------
# output/tables/Analysis_9e_llm_nonstructure.csv
# output/tables/llm_chemont_nonstructure_cache.csv
# output/figures/Analysis_9e_llm_nonstructure.pdf
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(jsonlite)
library(patchwork)
library(purrr)
library(readr)
library(reticulate)
library(stringr)
library(tidyr)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

PILOT_N      <- 500L    # set NULL for full run (~9,547 substances, ≈ $45)
PILOT_PRINT  <- FALSE   # set TRUE to print prompts/responses for debugging
CACHE_DIR    <- here("data", "cache", "llm_chemont_nonstructure")
HIERARCHY_TXT <- here("data", "processed", "chemont_hierarchy.txt")

dir.create(CACHE_DIR,                showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "figures"),showWarnings = FALSE, recursive = TRUE)

if (!nzchar(Sys.getenv("ANTHROPIC_API_KEY"))) {
  stop("ANTHROPIC_API_KEY environment variable not set.\n",
       "Set it with: Sys.setenv(ANTHROPIC_API_KEY = 'your-key') or via .Renviron")
}

if (!file.exists(HIERARCHY_TXT)) {
  stop("ChemOnt hierarchy text not found: ", HIERARCHY_TXT,
       "\nRun 09d_chemont_llm_experiment.R first to build it.")
}

# ------------------------------------------------------------------------------
# Step 1: Build non-structure substance list with matchable filter
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances_classified.rds"))

non_structure <- all_substances |>
  filter(is.na(inchikey)) |>
  filter(!entity_type %in% c("Regulatory entry")) |>
  distinct(substance_name) |>
  filter(!is.na(substance_name), nzchar(substance_name)) |>
  mutate(
    matchable = !str_detect(
      substance_name,
      regex("reaction mass|petroleum|distillate|UVCB", ignore_case = TRUE)
    )
  )

matchable_substances <- non_structure |> filter(matchable)

message(sprintf(
  "Non-structure substances: %d total, %d matchable, %d excluded by filter",
  nrow(non_structure), nrow(matchable_substances),
  nrow(non_structure) - nrow(matchable_substances)
))

if (!is.null(PILOT_N) && nrow(matchable_substances) > PILOT_N) {
  set.seed(42L)
  matchable_substances <- matchable_substances |> slice_sample(n = PILOT_N)
  message(sprintf("Pilot mode: sampled %d substances", PILOT_N))
}

# ------------------------------------------------------------------------------
# Step 2: Call LLM classifier
# ------------------------------------------------------------------------------

py_script       <- here("python", "llm_chemont_classify.py")
python          <- py_exe()
substances_json <- tempfile(fileext = ".json")
output_json     <- tempfile(fileext = ".json")

write_json(matchable_substances$substance_name, substances_json, auto_unbox = FALSE)

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
  rename(class = class_)

file.remove(substances_json, output_json)

message(sprintf("LLM results: %d rows", nrow(llm_raw)))

# Token usage summary
if ("usage" %in% names(llm_raw) && is.data.frame(llm_raw$usage)) {
  u <- llm_raw$usage
  message(sprintf(
    "Tokens — input: %s | cache_creation: %s | cache_read: %s | output: %s",
    format(sum(u$input_tokens,               na.rm = TRUE), big.mark = ","),
    format(sum(u$cache_creation_input_tokens, na.rm = TRUE), big.mark = ","),
    format(sum(u$cache_read_input_tokens,     na.rm = TRUE), big.mark = ","),
    format(sum(u$output_tokens,               na.rm = TRUE), big.mark = ",")
  ))
}

# ------------------------------------------------------------------------------
# Step 3: Export cache as flat CSV
# ------------------------------------------------------------------------------

cache_files <- list.files(CACHE_DIR, pattern = "\\.json$", full.names = TRUE)

cache_df <- map(cache_files, \(f) {
  parsed        <- fromJSON(f)
  usage         <- as_tibble(parsed$usage)
  parsed$usage  <- NULL
  parsed        <- map(parsed, \(x) if (is.null(x)) NA_character_ else x)
  bind_cols(as_tibble(parsed), usage)
}) |>
  list_rbind() |>
  rename(class = class_)

write_csv(cache_df, here("output", "tables", "llm_chemont_nonstructure_cache.csv"))
message(sprintf(
  "Cache exported: %d rows → output/tables/llm_chemont_nonstructure_cache.csv",
  nrow(cache_df)
))

# ------------------------------------------------------------------------------
# Step 4: Compare with embedding matches from 09_embedding_chemont.R
# ------------------------------------------------------------------------------

embedding_file <- here("output", "tables", "final_matches.csv")

results <- matchable_substances |>
  left_join(
    llm_raw |> select(
      substance_name,
      llm_kingdom       = kingdom,
      llm_superclass    = superclass,
      llm_class         = class,
      llm_subclass      = subclass,
      llm_direct_parent = direct_parent,
      confidence
    ),
    by = "substance_name"
  )

if (file.exists(embedding_file)) {
  embedding_matches <- read_csv(embedding_file, show_col_types = FALSE) |>
    select(substance_name, emb_match = best_match, emb_score = score)

  results <- results |>
    left_join(embedding_matches, by = "substance_name") |>
    mutate(
      # Agreement: LLM direct_parent matches embedding best_match
      agree      = !is.na(llm_direct_parent) & !is.na(emb_match) &
                   llm_direct_parent == emb_match,
      both_match = !is.na(llm_direct_parent) & !is.na(emb_match)
    )

  n_both  <- sum(results$both_match, na.rm = TRUE)
  n_agree <- sum(results$agree,      na.rm = TRUE)
  message(sprintf(
    "Both methods matched: %d substances; agreement: %d (%.1f%%)",
    n_both, n_agree, 100 * n_agree / max(n_both, 1)
  ))
} else {
  warning("Embedding matches file not found: ", embedding_file,
          "\nSkipping comparison with embedding results.")
  results <- results |>
    mutate(emb_match = NA_character_, emb_score = NA_real_,
           agree = NA, both_match = NA)
}

# ------------------------------------------------------------------------------
# Step 5: Quality summary
# ------------------------------------------------------------------------------

confidence_summary <- results |>
  count(confidence) |>
  mutate(pct = n / sum(n))

null_rates <- results |>
  summarise(
    kingdom    = mean(is.na(llm_kingdom)),
    superclass = mean(is.na(llm_superclass)),
    class      = mean(is.na(llm_class)),
    subclass   = mean(is.na(llm_subclass)),
    direct_parent = mean(is.na(llm_direct_parent))
  ) |>
  pivot_longer(everything(), names_to = "level", values_to = "null_rate") |>
  mutate(level = factor(level,
    levels = c("kingdom", "superclass", "class", "subclass", "direct_parent")))

message("\nConfidence distribution:")
print(confidence_summary)
message("\nNull rates per level:")
print(null_rates)

# ------------------------------------------------------------------------------
# Step 6: Write results CSV
# ------------------------------------------------------------------------------

write_csv(results, here("output", "tables", "Analysis_9e_llm_nonstructure.csv"))
message("Results written → output/tables/Analysis_9e_llm_nonstructure.csv")

# ------------------------------------------------------------------------------
# Step 7: Figures
#
# Panel 1 — Confidence distribution
# Panel 2 — Null rate per hierarchy level
# Panel 3 — LLM vs. embedding agreement (only for substances with both results)
# ------------------------------------------------------------------------------

conf_levels <- c("high", "medium", "low")

p1 <- confidence_summary |>
  filter(!is.na(confidence)) |>
  mutate(confidence = factor(confidence, levels = conf_levels)) |>
  ggplot(aes(x = confidence, y = pct, fill = confidence)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(n = %d)", 100 * pct, n)),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c(high = "#2166ac", medium = "#f4a582", low = "#d6604d"),
                    guide = "none") +
  theme_minimal(base_size = 12) +
  labs(
    title    = "LLM confidence distribution",
    subtitle = sprintf("n = %d matchable non-structure substances", nrow(results)),
    x = NULL, y = "Fraction"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

p2 <- null_rates |>
  ggplot(aes(x = level, y = null_rate)) +
  geom_col(fill = "#878787", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * null_rate)),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Null rate per hierarchy level",
    subtitle = "Higher null rate at specific levels reflects appropriate uncertainty",
    x = NULL, y = "Fraction null"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

if (file.exists(embedding_file) && any(results$both_match, na.rm = TRUE)) {
  agree_summary <- results |>
    filter(both_match) |>
    count(confidence, agree) |>
    group_by(confidence) |>
    mutate(pct = n / sum(n),
           confidence = factor(confidence, levels = conf_levels))

  p3 <- agree_summary |>
    filter(agree) |>
    ggplot(aes(x = confidence, y = pct, fill = confidence)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * pct)), vjust = -0.3, size = 3.5) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c(high = "#2166ac", medium = "#f4a582", low = "#d6604d"),
                      guide = "none") +
    theme_minimal(base_size = 12) +
    labs(
      title    = "Agreement with embedding match, per confidence level",
      subtitle = sprintf(
        "Substances with both LLM and embedding result: n = %d",
        sum(results$both_match, na.rm = TRUE)
      ),
      x = NULL, y = "Fraction agreeing"
    ) +
    theme(plot.subtitle = element_text(colour = "grey40"))
} else {
  p3 <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Embedding matches not available\n(final_matches.csv not found)",
             colour = "grey50", size = 5) +
    theme_void()
}

p_final <- (p1 | p2 | p3) +
  plot_layout(widths = c(1, 1, 1))

ggsave(p_final,
       filename = here("output", "figures", "Analysis_9e_llm_nonstructure.pdf"),
       device = "pdf", height = 20, width = 40, units = "cm")

message("09e_llm_nonstructure.R: completed")
