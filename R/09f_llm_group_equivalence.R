# ==============================================================================
# 09f_llm_group_equivalence.R
# LLM-based ChemOnt equivalence assessment for non-structure substances
#
# PURPOSE
# -------
# Experiments 09d and 09e asked "what ChemOnt class does this substance belong
# to?" (membership).  This script asks a different question: "is this substance
# name itself equivalent to a ChemOnt node?" (equivalence).
#
# When a legislative entry such as "cresols", "fatty alcohols", or "nickel and
# its compounds" maps to a ChemOnt node, the taxonomy can be used to enumerate
# all child substances automatically covered by that entry.
#
# Three relationship types are distinguished:
#   exact   — substance name IS the ChemOnt class  ("cresols" ≡ "Cresols")
#   subset  — name covers a specific subset         ("fatty acids C8-10" ⊂ ...)
#   broader — name spans multiple ChemOnt nodes     ("heavy metals" > any node)
#
# METHODOLOGY COMPARISON
# ----------------------
# entity_type = "Substance group" in all_substances_classified.rds is a
# rule-based (regex) classification of group entries.  The LLM's is_group
# field provides an independent assessment.  The cross-tabulation of these two
# reveals regex over- and under-coverage.
#
# COST ESTIMATE (claude-sonnet-4-6)
# ----------------------------------
# System prompt: ~31K tokens from full SKOS (cached after first call)
# Per PILOT_N = 500 substances:
#   Cache write:  31K × $3/Mtok              ≈ $0.09
#   Cache reads:  500 × 31K × $0.30/Mtok    ≈ $4.65
#   Output:       500 × ~150tok × $15/Mtok   ≈ $1.13
#   Total pilot:                              ≈ $6
# Full run (~9,547 matchable substances):     ≈ $130
# Set PILOT_N <- NULL to run the full set.
#
# DATA PROVENANCE
# ---------------
# Input 1: data/processed/all_substances_classified.rds
# Input 2: data/processed/chemont_hierarchy_full.txt  (full SKOS, 4824 labels)
# Input 3: output/tables/final_matches.csv  (embedding matches from 09)
# Cache:   data/cache/llm_chemont_equivalence/
#
# OUTPUTS
# -------
# output/tables/Analysis_9f_group_equivalence.csv
# output/tables/llm_chemont_equivalence_cache.csv
# output/figures/Analysis_9f_group_equivalence.pdf
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

PILOT_N       <- 5000L   # set NULL for full run (~9,547 substances, ≈ $130)
PILOT_PRINT   <- FALSE   # set TRUE to print prompts/responses for debugging
CACHE_DIR     <- here("data", "cache", "llm_chemont_equivalence")
HIERARCHY_TXT <- here("data", "processed", "chemont_hierarchy_full.txt")

dir.create(CACHE_DIR,                showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "figures"),showWarnings = FALSE, recursive = TRUE)

if (!nzchar(Sys.getenv("ANTHROPIC_API_KEY"))) {
  stop("ANTHROPIC_API_KEY environment variable not set.\n",
       "Set it with: Sys.setenv(ANTHROPIC_API_KEY = 'your-key') or via .Renviron")
}

if (!file.exists(HIERARCHY_TXT)) {
  stop("ChemOnt hierarchy text not found: ", HIERARCHY_TXT,
       "\nRun: python3 python/build_chemont_hierarchy.py ",
       "data/source/chemont/ChemOnt_2_1_skos.ttl ",
       "data/processed/chemont_hierarchy_full.txt")
}

# ------------------------------------------------------------------------------
# Step 1: Build substance list
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances_classified.rds"))

# entity_type lookup: one row per substance_name, preferring "Substance group"
# when a name appears under multiple entity_types
entity_type_lookup <- all_substances |>
  filter(is.na(inchikey)) |>
  filter(!entity_type %in% c("Regulatory entry")) |>
  distinct(substance_name, entity_type) |>
  group_by(substance_name) |>
  arrange(desc(entity_type == "Substance group")) |>
  slice(1) |>
  ungroup()

non_structure <- entity_type_lookup |>
  filter(!is.na(substance_name), nzchar(substance_name)) |>
  mutate(
    matchable = !str_detect(
      substance_name,
      regex("reaction mass|petroleum|distillate|UVCB", ignore_case = TRUE)
    )
  )

matchable_substances <- non_structure |>
  filter(matchable) |>
  filter(entity_type != "Analytical parameter")

message(sprintf(
  "Non-structure substances: %d total, %d matchable, %d excluded by matchable filter",
  nrow(non_structure),
  nrow(non_structure |> filter(matchable)),
  nrow(non_structure) - nrow(non_structure |> filter(matchable))
))

message(sprintf(
  "entity_type breakdown (matchable): %s",
  paste(
    matchable_substances |>
      count(entity_type) |>
      mutate(s = sprintf("%s: %d", entity_type, n)) |>
      pull(s),
    collapse = " | "
  )
))

if (!is.null(PILOT_N) && nrow(matchable_substances) > PILOT_N) {
  set.seed(42L)
  matchable_substances <- matchable_substances |> slice_sample(n = PILOT_N)
  message(sprintf("Pilot mode: sampled %d substances", PILOT_N))
}

# ------------------------------------------------------------------------------
# Step 2: Call LLM equivalence classifier
# ------------------------------------------------------------------------------

py_script       <- here("python", "llm_chemont_equivalence.py")
python          <- py_exe()
substances_json <- tempfile(fileext = ".json")
output_json     <- tempfile(fileext = ".json")

write_json(matchable_substances$substance_name, substances_json, auto_unbox = FALSE)

pilot_flag <- if (isTRUE(PILOT_PRINT)) "--pilot" else character(0)

message("Calling LLM equivalence classifier ...")
ret <- system2(
  python,
  args   = c(py_script, substances_json, HIERARCHY_TXT, CACHE_DIR, output_json, pilot_flag),
  stdout = "",
  stderr = ""
)
if (ret != 0L) stop("Python LLM classifier failed (exit code ", ret, ")")

llm_raw <- fromJSON(output_json, simplifyDataFrame = TRUE) |>
  as_tibble()

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
  parsed       <- fromJSON(f)
  usage        <- as_tibble(parsed$usage)
  parsed$usage <- NULL
  parsed       <- map(parsed, \(x) if (is.null(x)) NA_character_ else as.character(x))
  bind_cols(as_tibble(parsed), usage)
}) |>
  list_rbind() |>
  mutate(is_group = as.logical(is_group))

write_csv(cache_df, here("output", "tables", "llm_chemont_equivalence_cache.csv"))
message(sprintf(
  "Cache exported: %d rows → output/tables/llm_chemont_equivalence_cache.csv",
  nrow(cache_df)
))

# ------------------------------------------------------------------------------
# Step 4: Build results table
# ------------------------------------------------------------------------------

results <- matchable_substances |>
  left_join(
    llm_raw |> select(
      substance_name,
      is_group,
      equivalent_node,
      equivalent_level,
      relationship_type,
      confidence
    ),
    by = "substance_name"
  )

# Validate against embedding matches (final_matches.csv)
embedding_file <- here("output", "tables", "final_matches.csv")

if (file.exists(embedding_file)) {
  embedding_matches <- read_csv(embedding_file, show_col_types = FALSE) |>
    select(substance_name, emb_match = best_match, emb_score = score)

  results <- results |>
    left_join(embedding_matches, by = "substance_name") |>
    mutate(
      emb_agrees = !is.na(equivalent_node) & !is.na(emb_match) &
                   str_to_lower(equivalent_node) == str_to_lower(emb_match)
    )
} else {
  warning("Embedding matches file not found: ", embedding_file,
          "\nSkipping validation against embedding results.")
  results <- results |>
    mutate(emb_match = NA_character_, emb_score = NA_real_, emb_agrees = NA)
}

# ------------------------------------------------------------------------------
# Step 5: Quality summary
# ------------------------------------------------------------------------------

message("\n=== is_group distribution ===")
results |>
  count(is_group) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  print()

message("\n=== relationship_type for is_group == TRUE ===")
results |>
  filter(is_group == TRUE) |>
  count(relationship_type, confidence) |>
  arrange(relationship_type, confidence) |>
  print()

message("\n=== equivalent_level for is_group == TRUE ===")
results |>
  filter(is_group == TRUE) |>
  count(equivalent_level) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  print()

message("\n=== LLM vs regex: cross-tabulation ===")
results |>
  mutate(
    regex_group = entity_type == "Substance group",
    llm_group   = is_group == TRUE
  ) |>
  count(regex_group, llm_group) |>
  print()

if (file.exists(embedding_file)) {
  n_both  <- sum(!is.na(results$emb_match), na.rm = TRUE)
  n_agree <- sum(results$emb_agrees, na.rm = TRUE)
  message(sprintf(
    "\nEmbedding validation: %d substances with embedding match; LLM agrees on node: %d (%.1f%%)",
    n_both, n_agree, 100 * n_agree / max(n_both, 1)
  ))
}

# ------------------------------------------------------------------------------
# Step 6: Write results CSV
# ------------------------------------------------------------------------------

write_csv(results, here("output", "tables", "Analysis_9f_group_equivalence.csv"))
message("Results written → output/tables/Analysis_9f_group_equivalence.csv")

# ------------------------------------------------------------------------------
# Step 7: Figures
#
# Panel 1 — is_group distribution with confidence overlay
# Panel 2 — relationship_type for is_group == TRUE
# Panel 3 — LLM vs regex cross-tabulation (methodology comparison)
# ------------------------------------------------------------------------------

conf_levels <- c("high", "medium", "low")
rel_levels  <- c("exact", "subset", "broader")

# Panel 1: overall is_group × confidence
p1_data <- results |>
  mutate(
    is_group_label = case_when(
      isTRUE(is_group)  ~ "Group",
      isFALSE(is_group) ~ "Not a group",
      TRUE              ~ "Undetermined"
    ),
    confidence = factor(confidence, levels = conf_levels)
  ) |>
  count(is_group_label, confidence)

p1 <- p1_data |>
  ggplot(aes(x = is_group_label, y = n, fill = confidence)) +
  geom_col(width = 0.6) +
  scale_fill_manual(
    values = c(high = "#2166ac", medium = "#f4a582", low = "#d6604d"),
    na.value = "grey70", name = "Confidence"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "LLM group assessment",
    subtitle = sprintf("n = %d matchable non-structure substances", nrow(results)),
    x = NULL, y = "Count"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

# Panel 2: relationship_type for is_group == TRUE
p2_data <- results |>
  filter(is_group == TRUE) |>
  count(relationship_type, confidence) |>
  mutate(
    relationship_type = factor(relationship_type, levels = rel_levels),
    confidence        = factor(confidence, levels = conf_levels)
  )

p2 <- p2_data |>
  ggplot(aes(x = relationship_type, y = n, fill = confidence)) +
  geom_col(width = 0.6) +
  scale_fill_manual(
    values = c(high = "#2166ac", medium = "#f4a582", low = "#d6604d"),
    na.value = "grey70", name = "Confidence"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Equivalence relationship type",
    subtitle = sprintf(
      "Substances identified as groups: n = %d",
      sum(isTRUE(results$is_group), na.rm = TRUE)
    ),
    x = NULL, y = "Count"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

# Panel 3: LLM vs regex cross-tabulation
p3_data <- results |>
  mutate(
    regex_label = if_else(entity_type == "Substance group", "Regex: group", "Regex: other"),
    llm_label   = case_when(
      isTRUE(is_group)  ~ "LLM: group",
      isFALSE(is_group) ~ "LLM: not group",
      TRUE              ~ "LLM: undetermined"
    )
  ) |>
  count(regex_label, llm_label)

p3 <- p3_data |>
  ggplot(aes(x = regex_label, y = n, fill = llm_label)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(
    values = c(
      "LLM: group"        = "#2166ac",
      "LLM: not group"    = "#d6604d",
      "LLM: undetermined" = "grey70"
    ),
    name = NULL
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "LLM vs. regex methodology comparison",
    subtitle = "entity_type 'Substance group' (regex) vs. LLM is_group",
    x = NULL, y = "Count"
  ) +
  theme(
    plot.subtitle = element_text(colour = "grey40"),
    legend.position = "bottom"
  )

p_final <- (p1 | p2 | p3) + plot_layout(widths = c(1, 1, 1))

ggsave(p_final,
       filename = here("output", "figures", "Analysis_9f_group_equivalence.pdf"),
       device = "pdf", height = 20, width = 40, units = "cm")

message("09f_llm_group_equivalence.R: completed")
