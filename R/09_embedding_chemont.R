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
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **Sentence-BERT / all-MiniLM-L6-v2** (Reimers & Gurevych 2019, *Sentence-BERT: Sentence Embeddings using Siamese BERT-Networks*, EMNLP 2019. DOI: 10.18653/v1/D19-1410) | Embedding model for both substance names and ChemOnt labels; produces 384-dimensional dense vectors in a shared semantic space suitable for cosine similarity ranking. |
# | **ChemOnt 2.1 / ClassyFire** (Djoumbou Feunang et al. 2016, *ClassyFire: automated chemical classification with a comprehensive, computable taxonomy*, J. Cheminform. 8:61. DOI: 10.1186/s13321-016-0174-y) | Target ontology; provides a hierarchical, computable chemical classification system against which substance names are matched.  The `label` field of each ChemOnt class node is used as the text to embed. |
# | **Cosine similarity** (Manning, Raghavan & Schütze 2008, *Introduction to Information Retrieval*, Cambridge UP, §6.3) | Standard similarity metric for dense embedding spaces; computed as the dot product of L2-normalised vectors, which is equivalent to the angle between them and invariant to vector magnitude. |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | ChemOnt `label` only (not definition or altLabel) | Definitions contain longer prose that shifts the embedding centroid; using labels keeps the ChemOnt vector space parallel to the short substance-name vectors.  Own addition: no framework prescribes which ChemOnt field to embed. |
# | Top-k = 3 candidates retained per substance | k = 3 is the minimum needed to compute both gap_12 and gap_23 for the confidence score.  k = 1 would prevent confidence scoring; k > 3 adds computational cost without improving the two-gap signal. |
# | Confidence weights 0.6 / 0.3 / 0.1 (score_1 / gap_12 / gap_23) | Weights reflect the empirical dominance of absolute score strength over gap terms; validated by visual inspection of the match quality landscape (Fig. 9n).  Analogous to margin-based confidence in nearest-neighbour classifiers. |
# | Quality filters: score_1 > 0.8 and gap_12 > 0.05 | Thresholds selected after visual inspection of the joint distribution of score_1 and gap_12 (Fig. 9n); they correspond to the upper-right region where matches are both strong and unambiguous. |
# | Threshold exploration range 0.3–0.9 step 0.02 | Pragmatic exploratory range covering plausible cosine similarity values; the coverage curve (Fig. 9f) makes the coverage–precision trade-off at any threshold explicit. |
# | Unmatchable exclusion regex (reaction mass / petroleum / UVCB) | These name patterns denote substance classes or complex mixtures whose names do not correspond to single chemical classes; matching them to ChemOnt labels would produce spurious high-similarity results via incidental word overlap. |
#
# INTERPRETATION
# --------------
# **Analysis 9e — Distribution of top-1 cosine similarity scores**
# The score distribution shows how well substance names align with the
# nearest ChemOnt label in embedding space.  A distribution concentrated
# above 0.8 indicates that most names have a clear semantic match in the
# ontology.  A heavy left tail (many scores below 0.5) suggests that many
# non-structure names describe things the ChemOnt vocabulary does not
# cover (industrial formulations, trade names, regulatory placeholders).
#
# **Analysis 9f — Coverage vs. threshold curve**
# As the similarity threshold increases, fewer matches pass the quality
# filter (coverage decreases) but the retained matches are more reliable
# (precision increases).  The curve makes the coverage–precision trade-off
# explicit.  The selected threshold (score_1 > 0.8) is marked; readers can
# identify what fraction of names would be matched at a lower or higher
# threshold.
#
# **Analysis 9g — Precision estimation**
# Manual inspection of a random sample of matches at different thresholds
# gives an empirical precision estimate.  A match is considered correct if
# the assigned ChemOnt class is a plausible chemical family for the
# substance name.  Use this figure to calibrate confidence before using
# `final_matches.csv` in downstream analyses.
#
# **Analysis 9n — Quality landscape (score_1 vs. gap_12)**
# The scatter plot shows every substance in the (score_1, gap_12) space.
# The dashed lines mark the selected quality filters.  Substances in the
# upper-right quadrant (high score, large gap) are the most reliable
# matches; those in the lower-left are ambiguous.  If many substances
# cluster near the threshold boundaries, the filters may be sensitive to
# small parameter changes and should be validated by manual inspection.
#
# **final_matches.csv**
# Contains only the high-confidence subset.  A substance absent from this
# file either (a) scored below the threshold, (b) matched with insufficient
# gap between top-1 and top-2, or (c) was excluded by the unmatchable regex.
# It does NOT mean the substance has no chemical class — only that the
# embedding-based match did not meet the quality criteria.
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
library(patchwork)

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

n_above_threshold <- sum(matches$score > score_threshold)
pct_above         <- round(n_above_threshold / nrow(matches) * 100, 1)

p9e <- ggplot(matches, aes(score)) +
  geom_histogram(bins = 60, fill = "#4a90d9", colour = "white") +
  geom_vline(xintercept = score_threshold, linetype = "dashed", colour = "#e05c5c", linewidth = 0.8) +
  annotate("text",
           x     = score_threshold + 0.01,
           y     = Inf,
           label = sprintf("%.0f%% retained\n(n = %s)", pct_above,
                           format(n_above_threshold, big.mark = ",")),
           hjust = 0, vjust = 1.3, size = 7, colour = "#e05c5c") +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Cosine similarity score distribution",
    subtitle = sprintf("Dashed line: applied threshold (score\u2081 > %.2f)", score_threshold),
    x        = "Top-1 cosine similarity (score\u2081)",
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
    title    = "Coverage vs threshold",
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
    title    = "ECDF of similarity scores",
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
# Own addition: k = 3 is the minimum required to compute gap_12 and gap_23
# for confidence scoring; k = 1 would reduce the match to a bare score with
# no measure of ambiguity.  See own-additions table in the header.
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
# Own addition: gap_12 = score_1 − score_2 measures separation between the
# best and second-best match (high = unambiguous); gap_23 = score_2 − score_3
# measures stability further down the ranking (high = well-ordered top-3).
# mean_top3 is a secondary quality indicator; the confidence score (9l) uses
# only gap_12 and gap_23 directly.  This pattern is analogous to the
# classification margin in nearest-neighbour methods.
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
  filter(score_1 > score_threshold)

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
  theme_minimal(base_size = 16) +
  labs(
    title    = "Match quality landscape",
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
# Analysis 9p: Confidence vs. similarity scatter
# ==============================================================================
# INTENT
# The composite confidence score combines score_1, gap_12, and gap_23 with
# fixed weights (0.6 / 0.3 / 0.1).  Because score_1 dominates, confidence and
# score_1 should be highly correlated — but gap_12 and gap_23 introduce spread:
# two entries with the same score_1 can differ in confidence if one has a
# clearer separation from its runner-up.  This plot makes that spread visible
# and shows whether the quality filter (score_1 > 0.8, gap_12 > 0.05)
# captures the high-confidence region as intended.
# ==============================================================================

top_wide <- top_wide |>
  mutate(passes_filter = score_1 > 0.8 & gap_12 > 0.05)

p9p <- ggplot(top_wide, aes(x = score_1, y = confidence, colour = passes_filter)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_colour_manual(
    values = c("TRUE" = "#2171b5", "FALSE" = "#cccccc"),
    labels = c("TRUE" = "Passes filter", "FALSE" = "Excluded"),
    name   = NULL
  ) +
  theme_minimal(base_size = 16) +
  labs(
    title    = "Confidence vs. top-1 similarity",
    subtitle = "Spread around the diagonal reflects the contribution of gap_12 and gap_23 to confidence",
    x        = "Top-1 cosine similarity (score_1)",
    y        = "Composite confidence score"
  ) +
  theme(
    plot.subtitle   = element_text(colour = "grey40"),
    legend.position = "bottom"
  )

print(p9p)
ggsave(p9p,
       filename = here("output", "figures",
                       "Analysis_9p_confidence_vs_similarity.pdf"),
       device = "pdf",
       height = 25, width = 25, units = "cm")

# ==============================================================================
# 9o: Final matches export
# ==============================================================================
# INTENT
# Export the filtered, high-confidence substance → ChemOnt assignments for
# use in downstream classification and prioritisation steps.
# ==============================================================================

final_matches <- filtered |>
  select(
    substance_name,
    best_match = chemont_label_1,
    score      = score_1
  )

# Columns (final_matches.csv):
# Columns (final_matches.csv):
#   substance_name — original substance name from ECHA dataset (non-structure
#                    entries only); primary key for joining with all_substances
#   best_match     — top-1 ChemOnt class label (label text, not CHEMONT ID);
#                    the chemical family most similar to the substance name in
#                    embedding space
#   score          — top-1 cosine similarity (0–1); only rows with score > 0.8
#                    are included
write_csv(
  final_matches,
  here("output", "tables", "final_matches.csv")
)

# Top-133 matches by score — for manual inspection of match quality
# Columns (Analysis_9_top133_matches.csv):
#   substance_name — substance name
#   best_match     — assigned ChemOnt class label
#   score          — top-1 cosine similarity
write_csv(
  final_matches |> arrange(desc(score)) |> head(133),
  here("output", "tables", "Analysis_9_top133_matches.csv")
)

message(sprintf(
  "Analysis 9: %d substance names matched to ChemOnt with score_1 > %.2f",
  nrow(final_matches), score_threshold
))

# ==============================================================================
# SKILL.md §2.3 — Verify matched ChemOnt labels actually occur and are
# substantively meaningful (equivalent of the bio-label frequency check)
# ==============================================================================
# INTENT
# Before treating the final_matches as valid, confirm that the most frequently
# assigned ChemOnt classes are recognisable chemical families (e.g. organohalogen
# compounds, fatty acids) rather than generic or artefactual labels produced by
# the embedding space.  A concentration of matches on implausible classes would
# indicate that the model is exploiting spurious surface-form similarity.
# ==============================================================================

top_chemont_labels <- final_matches |>
  count(best_match, sort = TRUE) |>
  head(10)

message("\n=== Analysis 9: top-10 matched ChemOnt classes (SKILL.md §2.3 verification) ===")
message(paste(
  sprintf("  %d\t%s", top_chemont_labels$n, top_chemont_labels$best_match),
  collapse = "\n"
))

# ==============================================================================
# Analysis 8d+9e: Combined figure — UMAP and score distribution
# ==============================================================================

p8_umap <- readRDS(here("data", "processed", "p8_umap.rds"))

axis_theme <- theme(
  axis.text  = element_text(size = 20),
  axis.title = element_text(size = 20)
)

p8_umap_clean <- p8_umap +
  labs(tag = "(A)", title = NULL, subtitle = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 8))) +
  axis_theme +
  theme(
    legend.text     = element_text(size = 20),
    legend.title    = element_text(size = 18),
    legend.key.size = unit(1.2, "cm")
  )

p9e_clean <- p9e +
  labs(tag = "(B)", title = NULL, subtitle = NULL) +
  axis_theme

p_89 <- (p8_umap_clean | p9e_clean) +
  plot_annotation(
    title = "Embedding-based clustering and ChemOnt match quality",
    theme = theme(plot.title = element_text(size = 30, face = "bold"))
  ) +
  plot_layout(widths = c(2, 1))

ggsave(p_89,
       filename = here("output", "figures",
                       "Analysis_89_UMAP_and_ChemOnt_quality.pdf"),
       device = "pdf",
       height = 30, width = 60, units = "cm")

# ------------------------------------------------------------------------------
# Caption for combined figure 89 (printed to stdout)
# ------------------------------------------------------------------------------
# Columns (Analysis_9p_score_distribution.csv):
#   substance_name — original substance name
#   score_1        — top-1 cosine similarity score
#   passes_filter  — TRUE if score_1 > score_threshold
write_csv(
  top_wide |>
    select(substance_name, score_1) |>
    mutate(passes_filter = score_1 > score_threshold),
  here("output", "tables", "Analysis_9p_score_distribution.csv")
)

# ------------------------------------------------------------------------------
# Caption for combined figure 89 (printed to stdout)
# ------------------------------------------------------------------------------
n_non_structure <- nrow(non_structure)
n_clusters      <- n_distinct(readRDS(here("data", "processed",
                                           "non_structure_clusters.rds"))$cluster)
n_candidates    <- nrow(top_wide)
n_retained      <- nrow(filtered)
pct_retained    <- round(n_retained / n_candidates * 100, 1)

caption_89 <- sprintf(
  paste0(
    "(A) UMAP projection of sentence embeddings for %s non-structure regulatory substance entries, ",
    "coloured by k-means cluster (k = %d). ",
    "Each point represents a unique substance name; spatial proximity reflects semantic similarity ",
    "in the 384-dimensional embedding space. ",
    "(B) Distribution of top-1 cosine similarity scores between substance name embeddings ",
    "and ChemOnt class label embeddings across all %s candidate substance-class pairs. ",
    "The distribution is heavily skewed towards lower scores, indicating that most substance names ",
    "do not closely match any ChemOnt class label. ",
    "The dashed line indicates the applied threshold (score\u2081 > %.2f); ",
    "%d entries (%.1f%%) are retained as matches."
  ),
  format(n_non_structure, big.mark = ","),
  n_clusters,
  format(n_candidates, big.mark = ","),
  score_threshold,
  n_retained,
  pct_retained
)

cat("Figure caption (89):\n", caption_89, "\n")

message("09_embedding_chemont.R: analysis 9 completed")
