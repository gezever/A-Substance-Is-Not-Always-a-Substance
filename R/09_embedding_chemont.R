# ==============================================================================
# 09_embedding_chemont.R
# Analysis 9: Cosine similarity matching of substance names to ChemOnt classes
#
# PURPOSE
# -------
# Non-structure substances (no InChIKey) cannot be classified via
# ClassyFire's structure-based API.  This analysis asks whether the *names*
# of these substances can be matched to ChemOnt ontology classes via semantic
# similarity of sentence embeddings.  A high-confidence match allows a
# substance group or UVCB entry to be tentatively assigned to a chemical class
# without a known structure.
#
# DATA PROVENANCE
# ---------------
# Input 1: `data/processed/embeddings_echa.rds`
#          Embeddings of non-structure substance names; produced and cached
#          by `08_embedding_clustering.R`.
# Input 2: `data/source/ChemOnt_2_1.rds`
#          ChemOnt class labels (label column used as text to embed).
# Input 3: `data/processed/all_substances_classified.rds`
#          Used to reconstruct the `non_structure` name list with the same
#          filter applied in analysis 8.
#
# METHODOLOGY
# -----------
# 1.  ChemOnt class labels are embedded with the same all-MiniLM-L6-v2 model
#     used for substance names (ensuring compatible vector spaces).
# 2.  Cosine similarity matrix is computed as the matrix product of L2-
#     normalised embedding matrices: sim(i,j) = ê_i · ĉ_j.
# 3.  For each substance, the top-3 ChemOnt matches by cosine score are
#     retained.
# 4.  A composite confidence score is derived from the absolute top-1 score,
#     the gap between ranks 1 and 2, and the gap between ranks 2 and 3.
#     Own addition: the 0.6/0.3/0.1 weight allocation is not mandated by any
#     framework; it reflects the empirical observation that score_1 is the
#     dominant signal and the gap_12 a secondary discriminator.
# 5.  Quality filters: score_1 > 0.8 and gap_12 > 0.05 retain only
#     high-confidence, unambiguous matches.  Own addition: thresholds chosen
#     after visual inspection of the score and gap distributions (9e, 9n).
#
# OUTPUTS
# -------
# output/figures/Analysis_9e_distribution_only_label.pdf
# output/figures/Analysis_9f_coverage_threshold.pdf
# output/figures/Analysis_9g_precision_estimation.pdf
# output/figures/Analysis_9n_quality_landscape.pdf
# output/tables/final_matches.csv
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(reticulate)
library(purrr)
library(stringr)

# ------------------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed",
                               "all_substances_classified.rds"))

# Reconstruct non_structure with the same filter as in 08_embedding_clustering.R
non_structure <- all_substances |>
  filter(is.na(inchikey)) |>
  filter(!entity_type %in% c("Regulatory entry")) |>
  distinct(substance_name) |>
  filter(!is.na(substance_name), nzchar(substance_name))

embeddings_echa <- readRDS(here("data", "processed", "embeddings_echa.rds"))

chemont_rds <- here("data", "source", "ChemOnt_2_1.rds")
chemont     <- readRDS(chemont_rds)
chemont$substance_name <- chemont$label
chemont <- chemont |> select("substance", "substance_name")

# ==============================================================================
# Analysis 9a: Embed ChemOnt class labels
# ==============================================================================
# INTENT
# Embed ChemOnt class labels in the same 384-dimensional space as the
# substance names, enabling cosine similarity to measure semantic proximity
# between a substance name and a class label.  Using labels only (not
# definitions or synonyms) keeps the vector space consistent with the
# substance-name embeddings.
# ==============================================================================

chemont_embeddings_file <- here("data", "processed",
                                "embeddings_chemont_only_label.rds")

if (!file.exists(chemont_embeddings_file)) {
  message(sprintf("Analysis 9: embedding %d ChemOnt class labels",
                  nrow(chemont)))

  st           <- reticulate::import("sentence_transformers")
  model_chemont <- st$SentenceTransformer("all-MiniLM-L6-v2")

  embeddings_chemont <- model_chemont$encode(
    chemont$substance_name,
    show_progress_bar = TRUE,
    convert_to_numpy  = TRUE
  )

  saveRDS(embeddings_chemont, chemont_embeddings_file)
}

embeddings_chemont  <- readRDS(chemont_embeddings_file)
emb_chemont_matrix  <- as.matrix(embeddings_chemont)

# ==============================================================================
# Analysis 9b/9c: Cosine similarity matrix and best match per substance
# ==============================================================================
# INTENT
# Compute all pairwise cosine similarities between substance names and
# ChemOnt class labels.  For each substance, identify the highest-scoring
# class (best_idx) and record the score.  L2 normalisation before matrix
# multiplication produces cosine similarity directly.
# ==============================================================================

normalize <- function(x) x / sqrt(rowSums(x^2))

chemont_norm <- normalize(as.matrix(embeddings_chemont))
echa_norm    <- normalize(as.matrix(embeddings_echa))

# Cosine similarity matrix: (n_substances × n_chemont_classes)
sim_matrix <- echa_norm %*% t(chemont_norm)

best_idx   <- max.col(sim_matrix)
best_score <- sim_matrix[cbind(seq_len(nrow(sim_matrix)), best_idx)]

matches <- tibble(
  substance_name = non_structure$substance_name,
  chemont_label  = chemont$substance_name[best_idx],
  score          = best_score
)

# ==============================================================================
# Analysis 9d: Exclude inherently unmatchable substance types
# ==============================================================================
# INTENT
# Reaction masses, petroleum fractions, and UVCB substances have names that do
# not correspond to single chemical classes; matching them to ChemOnt would
# produce spurious high-similarity results based on incidental word overlap.
# Flagging these allows downstream filtering without dropping them from the
# score distribution plots.
# ==============================================================================

matches <- matches |>
  mutate(
    matchable = !str_detect(
      substance_name,
      regex("reaction mass|petroleum|distillate|UVCB", ignore_case = TRUE)
    )
  )

# ==============================================================================
# Analysis 9e: Score distribution
# ==============================================================================
# INTENT
# Visualise the distribution of best-match cosine scores across all substance
# names.  A bimodal distribution (one mode near 0.5–0.6, one near 0.8–0.9)
# would suggest a natural threshold separating confident from uncertain matches.
# This plot informs the threshold choice in 9h.
# ==============================================================================

p9e <- ggplot(matches, aes(score)) +
  geom_histogram(bins = 60, fill = "#4a90d9", colour = "white") +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Analysis 9e: Cosine similarity score distribution",
    subtitle = "Distribution of best-match scores between substance names and ChemOnt labels",
    x        = "Cosine similarity (best match)",
    y        = "Number of substance names"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p9e)
ggsave(p9e,
       filename = here("output", "figures",
                       "Analysis_9e_distribution_only_label.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")

# ==============================================================================
# Analysis 9f: Coverage vs threshold
# ==============================================================================
# INTENT
# The threshold that defines a "match" trades off coverage (fraction matched)
# against precision (fraction of matches that are correct).  This curve shows
# how rapidly coverage drops as the threshold rises, supporting the threshold
# choice at 0.7 (retaining a practical fraction of matches).
# Own addition: threshold range 0.3–0.9 and step size 0.02 are exploratory
# choices; the curve justifies the final threshold selected in 9h.
# ==============================================================================

thresholds <- seq(0.3, 0.9, by = 0.02)

coverage <- map_df(thresholds, function(t) {
  tibble(
    threshold = t,
    n_matches = sum(matches$score >= t),
    coverage  = mean(matches$score >= t)
  )
})

p9f <- ggplot(coverage, aes(threshold, coverage)) +
  geom_line(colour = "#4a90d9") +
  geom_point(size = 2, colour = "#4a90d9") +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Analysis 9f: Coverage vs threshold",
    subtitle = "Fraction of substance names matched at each cosine similarity threshold",
    x        = "Similarity threshold",
    y        = "Fraction matched"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p9f)
ggsave(p9f,
       filename = here("output", "figures",
                       "Analysis_9f_coverage_threshold.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")

# ==============================================================================
# Analysis 9g: ECDF of similarity scores (precision estimation)
# ==============================================================================
# INTENT
# The empirical cumulative distribution function of similarity scores gives
# the fraction of substance names with a score at or below any given value.
# Reading the ECDF at a candidate threshold shows what fraction would be
# excluded — complementary to the coverage curve (9f).
# ==============================================================================

p9g <- ggplot(matches, aes(x = score)) +
  stat_ecdf(colour = "#4a90d9") +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Analysis 9g: ECDF of similarity scores",
    subtitle = "Cumulative fraction of substance names at or below each score",
    x        = "Cosine similarity",
    y        = "Cumulative fraction"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p9g)
ggsave(p9g,
       filename = here("output", "figures",
                       "Analysis_9g_precision_estimation.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")

# ==============================================================================
# Analysis 9h–9l: Top-3 matches, confidence scoring, and filtering
# ==============================================================================
# INTENT
# Retaining the top-3 ChemOnt candidates per substance (rather than just the
# best) allows a confidence score that penalises ambiguous assignments where
# two or more classes score almost equally.  The composite confidence score
# rewards high absolute similarity (score_1), clear separation from the second
# candidate (gap_12), and stability further down the ranking (gap_23).
# ==============================================================================

# 9i: Retrieve top-3 matches per substance
top_k <- 3L

top_idx <- t(apply(sim_matrix, 1, function(x)
  order(x, decreasing = TRUE)[seq_len(top_k)]))
top_scores <- t(apply(sim_matrix, 1, function(x)
  sort(x, decreasing = TRUE)[seq_len(top_k)]))

top_matches <- map_dfr(seq_len(nrow(top_idx)), function(i) {
  tibble(
    substance_name = non_structure$substance_name[i],
    chemont_label  = chemont$substance_name[top_idx[i, ]],
    score          = top_scores[i, ]
  )
})

# 9j: Wide format with rank columns
top_ranked <- top_matches |>
  group_by(substance_name) |>
  arrange(desc(score)) |>
  mutate(rank = row_number()) |>
  ungroup()

top_wide <- top_ranked |>
  select(substance_name, rank, score, chemont_label) |>
  pivot_wider(
    names_from  = rank,
    values_from = c(score, chemont_label),
    names_sep   = "_"
  )

# 9k: Consistency metrics
top_wide <- top_wide |>
  mutate(
    gap_12    = score_1 - score_2,
    gap_23    = score_2 - score_3,
    mean_top3 = (score_1 + score_2 + score_3) / 3
  )

# 9l: Composite confidence score
# Own addition: weights 0.6 / 0.3 / 0.1 reflect empirical importance of each
# component; validated by inspecting the ranked match quality landscape (9n).
top_wide <- top_wide |>
  mutate(
    confidence =
      score_1 * 0.6 +   # absolute similarity strength
      gap_12  * 0.3 +   # separation between best and second-best match
      gap_23  * 0.1     # stability further down the ranking
  )

# 9m: Quality filter
# Own addition: score_1 > 0.8 (high absolute quality) and gap_12 > 0.05
# (unambiguous best match); thresholds selected after visual inspection of 9n.
filtered <- top_wide |>
  filter(
    score_1 > 0.8,
    gap_12  > 0.05
  )

# ==============================================================================
# Analysis 9n: Match quality landscape
# ==============================================================================
# INTENT
# A scatter of gap_12 (x-axis) versus score_1 (y-axis) reveals the joint
# distribution of match strength and match clarity.  The applied filter
# (score_1 > 0.8, gap_12 > 0.05) corresponds to the upper-right quadrant of
# this plot; visualising it confirms the filter boundary is meaningful and
# not arbitrary relative to the data distribution.
# ==============================================================================

p9n <- ggplot(top_wide, aes(x = gap_12, y = score_1)) +
  geom_point(alpha = 0.4, colour = "#4a90d9", size = 1) +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "#e05c5c") +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "#e05c5c") +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Analysis 9n: Match quality landscape",
    subtitle = "Dashed lines show applied filters: score_1 > 0.8 and gap_12 > 0.05",
    x        = "Gap between top-1 and top-2 score",
    y        = "Top-1 cosine similarity score"
  ) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p9n)
ggsave(p9n,
       filename = here("output", "figures",
                       "Analysis_9n_quality_landscape.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")

# ==============================================================================
# 9o: Final matches export
# ==============================================================================
# INTENT
# Export the filtered, high-confidence substance → ChemOnt assignments for
# use in downstream classification and prioritisation steps.
# Columns:
#   substance_name — original substance name from ECHA dataset
#   best_match     — top-1 ChemOnt class label (label text, not CHEMONT ID)
#   score          — cosine similarity of the top-1 match (0–1)
#   confidence     — composite confidence score (0.6·score + 0.3·gap12 + 0.1·gap23)
# ==============================================================================

final_matches <- filtered |>
  select(
    substance_name,
    best_match = chemont_label_1,
    score      = score_1,
    confidence
  )

write.csv(
  final_matches,
  file      = here("output", "tables", "final_matches.csv"),
  row.names = FALSE
)

message(sprintf(
  "Analysis 9: %d substance names matched to ChemOnt with score > 0.8 and gap > 0.05",
  nrow(final_matches)
))

message("09_embedding_chemont.R: analysis 9 completed")
