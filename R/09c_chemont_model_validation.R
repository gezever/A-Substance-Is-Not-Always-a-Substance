# ==============================================================================
# 09c_chemont_model_validation.R
# Model validation via structure-defined substances as ground truth
#
# PURPOSE
# -------
# Structure-defined substances have a known ChemOnt class via ClassyFire
# (encoded in data/processed/rdf/substances_taxonomy.ttl as skos:broader).
# This provides an internal ground truth: for each such substance, we know
# which ChemOnt class is correct.
#
# By embedding the *names* of structure-defined substances and ranking ChemOnt
# classes by cosine similarity, we can measure how often each model places the
# correct class in the top-1, top-3, and top-5 positions — a standard
# information retrieval metric (Hit@k).
#
# This is a principled comparison that does not require external labelled data:
# the ground truth is derived from structure-based classification (ClassyFire),
# which is independent of the embedding models being evaluated.
#
# METHODOLOGY
# -----------
# 1. Parse substances_taxonomy.ttl to extract (inchikey → ChemOnt class label)
# 2. Join with all_substances to get (substance_name → correct_chemont_label)
# 3. For each model: embed substance names → compute cosine similarity with
#    ChemOnt class label embeddings → rank classes per substance
# 4. Compute Hit@1, Hit@3, Hit@5 and mean reciprocal rank (MRR)
#
# OUTPUTS
# -------
# output/figures/Analysis_9c_model_validation.pdf
# output/tables/Analysis_9c_model_validation.csv
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(reticulate)
library(purrr)
library(readr)
library(tidyr)
library(stringr)
library(patchwork)


# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

models <- list(
  list(id = "miniLM",  name = "all-MiniLM-L6-v2"),
  list(id = "mpnet",   name = "multi-qa-mpnet-base-v2"),
  list(id = "scibert", name = "allenai/scibert_scivocab_uncased"),
  list(id = "biobert", name = "dmis-lab/biobert-v1.1")
)

dir.create(here("data", "cache", "embeddings"), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Step 1: Load ground truth from SPARQL output
# Run bash/chemont_broader_substance.sh first to generate this file.
# Columns: inchikey, substance_name, chemont_label (direct parent),
#          chemont_labels (all ancestors, pipe-separated)
# ------------------------------------------------------------------------------
system("bash bash/chemont_broader_substance.sh")


ground_truth_file <- here("data", "processed", "chemont_ground_truth.csv")

if (!file.exists(ground_truth_file)) {
  stop(paste(
    "Ground truth file not found:", ground_truth_file,
    "\nRun bash/chemont_broader_substance.sh first."
  ))
}

ground_truth <- read_csv(ground_truth_file, show_col_types = FALSE) |>
  rename(
    inchikey        = inchikey,
    chemont_label   = chemont_label,
    chemont_labels  = chemont_labels
  ) |>
  mutate(
    # All ancestor labels as a character vector per row
    chemont_ancestors = str_split(chemont_labels, "\\|")
  ) |>
  filter(!is.na(inchikey), !is.na(chemont_label))

message(sprintf("Ground truth: %d substance-class pairs loaded",
                nrow(ground_truth)))

# ------------------------------------------------------------------------------
# Step 2: Join with substance names
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances_classified.rds"))

validation_set <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(substance_name, inchikey) |>
  inner_join(ground_truth, by = "inchikey") |>
  filter(!is.na(substance_name), nzchar(substance_name))

message(sprintf("Validation set: %d (substance_name, correct_chemont_label) pairs",
                nrow(validation_set)))

# ------------------------------------------------------------------------------
# Step 3: Load ChemOnt labels and helper functions
# ------------------------------------------------------------------------------

chemont     <- readRDS(here("data", "source", "ChemOnt_2_1.rds"))
chemont$substance_name <- chemont$label
chemont <- chemont |> select("substance", "substance_name")

normalize <- function(x) x / sqrt(rowSums(x^2))

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
# Step 4: Evaluate each model via Hit@k and MRR
# ------------------------------------------------------------------------------

evaluate_model <- function(model_id, model_name) {
  message(sprintf("\n=== Evaluating: %s (%s) ===", model_id, model_name))

  emb_substances <- embed_with_model(
    texts      = validation_set$substance_name,
    model_name = model_name,
    cache_file = here("data", "cache", "embeddings",
                      sprintf("%s_validation.rds", model_id))
  )

  emb_chemont <- embed_with_model(
    texts      = chemont$substance_name,
    model_name = model_name,
    cache_file = here("data", "cache", "embeddings",
                      sprintf("%s_chemont.rds", model_id))
  )

  sub_norm <- normalize(as.matrix(emb_substances))
  che_norm <- normalize(as.matrix(emb_chemont))
  sim      <- sub_norm %*% t(che_norm)

  # Strict: rank of the direct parent class
  correct_idx_strict <- match(validation_set$chemont_label, chemont$substance_name)

  ranks_strict <- map_int(seq_len(nrow(sim)), function(i) {
    if (is.na(correct_idx_strict[i])) return(NA_integer_)
    rank(-sim[i, ], ties.method = "min")[correct_idx_strict[i]]
  })

  # Lenient: best rank among all ancestor classes
  ranks_lenient <- map_int(seq_len(nrow(sim)), function(i) {
    ancestors <- validation_set$chemont_ancestors[[i]]
    idxs      <- match(ancestors, chemont$substance_name)
    idxs      <- idxs[!is.na(idxs)]
    if (length(idxs) == 0L) return(NA_integer_)
    min(rank(-sim[i, ], ties.method = "min")[idxs])
  })

  tibble(
    model_id       = model_id,
    model_name     = model_name,
    # Strict metrics (direct parent only)
    hit_at_1       = mean(ranks_strict == 1L,  na.rm = TRUE),
    hit_at_3       = mean(ranks_strict <= 3L,  na.rm = TRUE),
    hit_at_5       = mean(ranks_strict <= 5L,  na.rm = TRUE),
    hit_at_10      = mean(ranks_strict <= 10L, na.rm = TRUE),
    mrr            = mean(1 / ranks_strict,    na.rm = TRUE),
    # Lenient metrics (any ancestor class)
    hit_at_1_len   = mean(ranks_lenient == 1L,  na.rm = TRUE),
    hit_at_3_len   = mean(ranks_lenient <= 3L,  na.rm = TRUE),
    hit_at_5_len   = mean(ranks_lenient <= 5L,  na.rm = TRUE),
    hit_at_10_len  = mean(ranks_lenient <= 10L, na.rm = TRUE),
    mrr_len        = mean(1 / ranks_lenient,    na.rm = TRUE),
    n_eval         = sum(!is.na(ranks_strict))
  )
}

metrics <- map_dfr(models, ~ evaluate_model(.x$id, .x$name))

print(metrics)

# Columns (Analysis_9c_model_validation.csv):
#   model_id  — short model identifier
#   hit_at_1  — fraction where correct class is ranked #1
#   hit_at_3  — fraction where correct class is in top 3
#   hit_at_5  — fraction where correct class is in top 5
#   hit_at_10 — fraction where correct class is in top 10
#   mrr       — mean reciprocal rank (higher = better)
#   n_eval    — number of substances evaluated (correct class found in ChemOnt)
write_csv(metrics,
          here("output", "tables", "Analysis_9c_model_validation.csv"))

# ------------------------------------------------------------------------------
# Plot: Hit@k per model
# ------------------------------------------------------------------------------

metrics_long <- metrics |>
  select(model_id, hit_at_1, hit_at_3, hit_at_5, hit_at_10) |>
  pivot_longer(-model_id, names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric,
                         levels = c("hit_at_1", "hit_at_3", "hit_at_5", "hit_at_10"),
                         labels = c("Hit@1", "Hit@3", "Hit@5", "Hit@10")))

p_hits <- ggplot(metrics_long, aes(x = metric, y = value, fill = model_id)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set2", name = "Model") +
  theme_minimal(base_size = 12) +
  labs(
    title    = "ChemOnt class retrieval performance by model",
    subtitle = "Evaluated on structure-defined substances with known ClassyFire class (ground truth)",
    x        = NULL,
    y        = "Fraction correct"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

p_mrr <- ggplot(metrics, aes(x = reorder(model_id, mrr), y = mrr, fill = model_id)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = round(mrr, 3)), hjust = -0.1, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Mean Reciprocal Rank (MRR) by model",
    x     = NULL,
    y     = "MRR"
  )

p_val <- (p_hits / p_mrr) +
  plot_layout(heights = c(2, 1))

ggsave(p_val,
       filename = here("output", "figures", "Analysis_9c_model_validation.pdf"),
       device = "pdf", height = 25, width = 30, units = "cm")

message("09c_chemont_model_validation.R: validation completed")
