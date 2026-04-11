# ==============================================================================
# 09b_chemont_model_comparison.R
# Compares ChemOnt matching quality across multiple sentence transformer models
#
# PURPOSE
# -------
# 09_embedding_chemont.R uses all-MiniLM-L6-v2, a compact general-purpose
# model.  Manual inspection of the top-133 matches revealed that a majority
# are chemically incorrect, with high scores driven by surface-form overlap
# rather than semantic understanding.  This script re-runs the same matching
# pipeline for three alternative models and compares:
#   (a) Score distributions
#   (b) Number of entries retained above threshold
#   (c) Overlap in retained matches between models
#   (d) Top matches per model (for manual inspection)
#
# MODELS COMPARED
# ---------------
# | ID       | Model name                          | Rationale |
# |----------|-------------------------------------|-----------|
# | miniLM   | all-MiniLM-L6-v2                    | Baseline (current) |
# | mpnet    | multi-qa-mpnet-base-v2              | Higher quality, query-doc matching |
# | scibert  | allenai/scibert_scivocab_uncased    | Scientific text domain |
# | biobert  | dmis-lab/biobert-v1.1               | Biomedical/chemical names |
#
# DATA PROVENANCE
# ---------------
# Input:  data/processed/all_substances_classified.rds  (from 04)
#         data/source/ChemOnt_2_1.rds
# Cache:  data/cache/embeddings/<model_id>_echa.rds
#         data/cache/embeddings/<model_id>_chemont.rds
#
# OUTPUTS
# -------
# output/figures/Analysis_9b_score_distributions.pdf
# output/figures/Analysis_9b_retention_comparison.pdf
# output/tables/Analysis_9b_model_comparison.csv
# output/tables/Analysis_9b_top20_<model_id>.csv   (one per model)
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(reticulate)
library(purrr)
library(readr)
library(patchwork)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

models <- list(
  list(id = "miniLM",  name = "all-MiniLM-L6-v2"),
  list(id = "mpnet",   name = "all-mpnet-base-v2"),
  list(id = "scibert", name = "allenai/scibert_scivocab_uncased"),
  list(id = "biobert", name = "dmis-lab/biobert-v1.1")
)

score_threshold <- 0.80

dir.create(here("data", "cache", "embeddings"), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances_classified.rds"))

non_structure <- all_substances |>
  filter(is.na(inchikey)) |>
  filter(!entity_type %in% c("Regulatory entry")) |>
  distinct(substance_name) |>
  filter(!is.na(substance_name), nzchar(substance_name))

chemont     <- readRDS(here("data", "source", "ChemOnt_2_1.rds"))
chemont$substance_name <- chemont$label
chemont <- chemont |> select("substance", "substance_name")

normalize <- function(x) x / sqrt(rowSums(x^2))

# ------------------------------------------------------------------------------
# Helper: embed a character vector with a given model, with file-based cache
# ------------------------------------------------------------------------------

embed_with_model <- function(texts, model_name, cache_file) {
  if (file.exists(cache_file)) {
    message(sprintf("  Loading cached embeddings from %s", basename(cache_file)))
    return(readRDS(cache_file))
  }

  message(sprintf("  Embedding %d texts with %s", length(texts), model_name))
  st    <- reticulate::import("sentence_transformers")
  model <- st$SentenceTransformer(model_name)
  emb   <- model$encode(texts, show_progress_bar = TRUE, convert_to_numpy = TRUE)
  saveRDS(emb, cache_file)
  emb
}

# ------------------------------------------------------------------------------
# Helper: run matching pipeline for one model, return tidy results
# ------------------------------------------------------------------------------

run_matching <- function(model_id, model_name) {
  message(sprintf("\n=== Model: %s (%s) ===", model_id, model_name))

  emb_echa <- embed_with_model(
    texts      = non_structure$substance_name,
    model_name = model_name,
    cache_file = here("data", "cache", "embeddings",
                      sprintf("%s_echa.rds", model_id))
  )

  emb_chemont <- embed_with_model(
    texts      = chemont$substance_name,
    model_name = model_name,
    cache_file = here("data", "cache", "embeddings",
                      sprintf("%s_chemont.rds", model_id))
  )

  echa_norm    <- normalize(as.matrix(emb_echa))
  chemont_norm <- normalize(as.matrix(emb_chemont))
  sim_matrix   <- echa_norm %*% t(chemont_norm)

  best_idx   <- max.col(sim_matrix, ties.method = "first")
  best_score <- sim_matrix[cbind(seq_len(nrow(sim_matrix)), best_idx)]

  tibble(
    model_id       = model_id,
    model_name     = model_name,
    substance_name = non_structure$substance_name,
    best_match     = chemont$substance_name[best_idx],
    score_1        = best_score
  )
}

# ------------------------------------------------------------------------------
# Run all models
# ------------------------------------------------------------------------------

results <- map_dfr(models, ~ run_matching(.x$id, .x$name))

# ------------------------------------------------------------------------------
# Analysis 9b-i: Score distributions per model
# ------------------------------------------------------------------------------

p_dist <- ggplot(results, aes(x = score_1, fill = model_id)) +
  geom_histogram(bins = 60, colour = "white", alpha = 0.85) +
  geom_vline(xintercept = score_threshold, linetype = "dashed",
             colour = "#e05c5c", linewidth = 0.8) +
  facet_wrap(~ model_id, ncol = 2, scales = "free_y") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top-1 cosine similarity distributions by model",
    subtitle = sprintf("Dashed line: threshold s\u2081 > %.2f", score_threshold),
    x = "Top-1 cosine similarity (score\u2081)",
    y = "Number of substance names"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

ggsave(p_dist,
       filename = here("output", "figures", "Analysis_9b_score_distributions.pdf"),
       device = "pdf", height = 20, width = 30, units = "cm")

# ------------------------------------------------------------------------------
# Analysis 9b-ii: Retention counts above threshold
# ------------------------------------------------------------------------------

retention <- results |>
  group_by(model_id, model_name) |>
  summarise(
    n_total    = n(),
    n_retained = sum(score_1 > score_threshold),
    pct        = round(n_retained / n_total * 100, 2),
    median_score = round(median(score_1), 3),
    .groups = "drop"
  )

p_retention <- ggplot(retention,
                      aes(x = reorder(model_id, n_retained), y = n_retained,
                          fill = model_id)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n_retained, pct)),
            hjust = -0.1, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = sprintf("Entries retained above threshold (s\u2081 > %.2f) per model",
                    score_threshold),
    x = NULL, y = "Number of retained matches"
  )

ggsave(p_retention,
       filename = here("output", "figures", "Analysis_9b_retention_comparison.pdf"),
       device = "pdf", height = 12, width = 20, units = "cm")

# ------------------------------------------------------------------------------
# Analysis 9b-iii: Top-20 matches per model (for manual inspection)
# ------------------------------------------------------------------------------

walk(models, function(m) {
  results |>
    filter(model_id == m$id, score_1 > score_threshold) |>
    arrange(desc(score_1)) |>
    head(20) |>
    select(substance_name, best_match, score_1) |>
    write_csv(here("output", "tables",
                   sprintf("Analysis_9b_top20_%s.csv", m$id)))
})

# ------------------------------------------------------------------------------
# Summary table
# Columns:
#   model_id     — short identifier
#   model_name   — full HuggingFace model name
#   n_total      — total number of substance-class pairs evaluated
#   n_retained   — entries with score_1 > threshold
#   pct          — percentage retained
#   median_score — median top-1 score across all entries
# ------------------------------------------------------------------------------

write_csv(retention,
          here("output", "tables", "Analysis_9b_model_comparison.csv"))

print(retention)

message("09b_chemont_model_comparison.R: comparison completed")
