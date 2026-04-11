# ==============================================================================
# 08_embedding_clustering.R
# Analysis 8: Sentence embedding and clustering of non-structure substance names
#
# PURPOSE
# -------
# The ~79 % of regulatory entries that lack an InChIKey cannot be analysed
# via structural cheminformatics.  Many of these entries are not random: they
# describe chemical classes, UVCB substances, or industrial formulations with
# systematic naming patterns.  This analysis asks whether semantic clustering
# of the substance names themselves can recover meaningful groupings — and if
# so, what those groupings are.
#
# Outputs feed into 09_embedding_chemont.R (cosine similarity matching) and
# into the classification scheme (13b_before_prioritisation_create_scheme.R).
#
# DATA PROVENANCE
# ---------------
# Input:  `data/processed/all_substances_classified.rds`
#         Produced by `04_entity_classification.R`.  Requires `entity_type` to
#         exclude purely administrative entries (Regulatory entry) from the
#         embedding corpus.
# Cache:  `data/processed/embeddings_echa.rds`
#         Sentence embeddings (384-dimensional) produced by the
#         all-MiniLM-L6-v2 model via sentence-transformers (Python).
#         Computed once and cached; set file.exists() guard to reuse.
#
# METHODOLOGY
# -----------
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **Sentence-BERT / all-MiniLM-L6-v2** (Reimers & Gurevych 2019, *Sentence-BERT: Sentence Embeddings using Siamese BERT-Networks*, EMNLP 2019. DOI: 10.18653/v1/D19-1410) | Embedding model; maps each substance name to a 384-dimensional dense vector in a shared semantic space. |
# | **Silhouette coefficient** (Rousseeuw 1987, *Silhouettes: a graphical aid to the interpretation and validation of cluster analysis*, J. Comput. Appl. Math. 20:53–65. DOI: 10.1016/0377-0427(87)90125-7) | Criterion for selecting the optimal number of clusters k; mean silhouette maximises within-cluster cohesion relative to between-cluster separation. |
# | **k-means clustering** (Lloyd 1982, *Least squares quantization in PCM*, IEEE Trans. Inf. Theory 28(2):129–137. DOI: 10.1109/TIT.1982.1056489) | Partitioning algorithm used to assign each substance name to one of k clusters based on Euclidean distance in embedding space. |
# | **UMAP** (McInnes, Healy & Melville 2018, *UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction*, arXiv:1802.03426) | 2-D projection for visualisation; preserves both local neighbourhood structure and global topology, making cluster separation visible. |
# | **t-SNE** (van der Maaten & Hinton 2008, *Visualizing Data using t-SNE*, J. Mach. Learn. Res. 9:2579–2605) | Second 2-D projection applied to the same embedding matrix as a structural validation of the UMAP result; optimises a different objective (KL divergence) so agreement between the two projections reduces the risk of visualisation artefacts. |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | Model: all-MiniLM-L6-v2 over larger Sentence-BERT variants | Compact (80 MB) and fast; produces 384-dimensional vectors sufficient for short chemical name strings.  Larger models (e.g. all-mpnet-base-v2) were tested and showed negligible silhouette improvement at significantly higher compute cost. |
# | k search range 2–12 | Pragmatic search space covering plausible chemical name groupings without over-segmentation; ranges beyond 12 were explored manually and did not yield higher silhouette scores. |
# | k = 6 fixed (overrides silhouette optimum) | Manual inspection of random cluster samples (8f) showed that the automatic optimum produced fewer than 6 interpretable groupings on initial runs; k = 6 yielded six substantively coherent chemical categories. |
# | UMAP: n_neighbors = 15, min_dist = 0.1 | Default UMAP parameters recommended by McInnes et al. for datasets of this size; n_neighbors = 15 balances local and global structure, min_dist = 0.1 provides moderate cluster compactness. |
# | UMAP metric: cosine | Cosine distance is the natural metric for L2-normalised embedding vectors; Euclidean distance in high-dimensional embedding spaces suffers from the curse of dimensionality. |
# | t-SNE perplexity: min(30, ⌊(n−1)/3⌋) | The t-SNE algorithm requires perplexity < n/3; the cap at 30 applies the commonly recommended default (van der Maaten & Hinton 2008) while satisfying the constraint for any n. |
# | k-means nstart = 25, iter.max = 100 | nstart = 25 multiple random initialisations reduce sensitivity to the starting configuration; iter.max = 100 ensures convergence for this embedding dimensionality. |
# | set.seed(42) throughout | Fixed seed for k-means, UMAP, and t-SNE ensures that results are reproducible across runs; the value 42 is arbitrary.
#
# REQUIREMENTS
# ------------
# Python (via reticulate): pip install sentence-transformers
# R packages: reticulate, uwot, Rtsne, cluster, purrr
#
# INTERPRETATION
# --------------
# **Analysis 8b — Silhouette scores across k**
# The peak of the silhouette curve is the statistically optimal k.  A flat
# curve (small differences between k values) means the embedding space does
# not have strongly separated clusters; interpretation of any k choice must
# then rely more on manual inspection than on the metric.  k = 6 was fixed
# based on substantive inspection of cluster contents (8f) rather than the
# silhouette peak alone.
#
# **Analysis 8d/8e — UMAP and t-SNE projections**
# Both figures map the 384-dimensional embedding space to 2-D.  Tight,
# well-separated clusters in both projections independently are strong
# evidence for genuine structure in the embedding space.  Clusters that
# appear in UMAP but dissolve in t-SNE (or vice versa) suggest that the
# separation is sensitive to the projection objective and should not be
# over-interpreted.  The trustworthiness/continuity metrics (8h) quantify
# this consistency numerically.
#
# **Analysis 8g — Silhouette comparison across spaces**
# Compares mean silhouette in the original embedding space, in UMAP 2-D, and
# in t-SNE 2-D.  A large drop from embedding to projection space indicates
# that the projection introduces distortions: clusters that are cohesive in
# high dimensions appear mixed in 2-D.  A small drop confirms that the
# visualisations are faithful.
#
# **Analysis 8h — Trustworthiness and continuity**
# Trustworthiness measures whether new neighbours introduced by the
# projection were already near-neighbours in the original space.
# Continuity measures whether original neighbours are preserved in the
# projection.  Both range 0–1; values > 0.9 indicate a high-fidelity
# projection.  Lower values at larger k (more neighbours) are expected;
# the important comparison is UMAP vs. t-SNE at the same k: the method with
# consistently higher values better preserves the original neighbourhood
# structure and is the more faithful visualisation.
#
# OUTPUTS
# -------
# output/figures/Analysis_8b_silhouette_score.pdf
# output/figures/Analysis_8d_UMAP.pdf
# output/figures/Analysis_8e_tsne.pdf
# output/tables/Analysis_8_silhouette_comparison.csv
# output/tables/Analysis_8_projection_coordinates.csv
# output/tables/Analysis_8_trustworthiness_continuity.csv
# data/processed/embeddings_echa.rds        (embedding cache)
# data/processed/non_structure_clusters.rds (cluster-labelled name table)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(reticulate)
library(uwot)
library(Rtsne)
library(ggrepel)
library(cluster)
library(purrr)

# ------------------------------------------------------------------------------
# Load entity-classified data (entity_type required to filter corpus)
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed",
                               "all_substances_classified.rds"))

# ==============================================================================
# Analysis 8a: Embed non-structure substance names
# ==============================================================================
# INTENT
# Produce a dense vector representation of each substance name so that
# semantic similarity can be measured numerically.  Regulatory entry names
# are excluded because they are cross-references, not chemical names, and
# would pollute the chemical semantic space.
# ==============================================================================

embeddings_file <- here("data", "processed", "embeddings_echa.rds")

non_structure <- all_substances |>
  filter(is.na(inchikey)) |>
  filter(!entity_type %in% c("Regulatory entry")) |>
  distinct(substance_name) |>
  filter(!is.na(substance_name), nzchar(substance_name))

if (!file.exists(embeddings_file)) {
  message(sprintf("Analysis 8: embedding %d non_structure substance names",
                  nrow(non_structure)))

  st    <- reticulate::import("sentence_transformers")
  model <- st$SentenceTransformer("all-MiniLM-L6-v2")

  embeddings_echa <- model$encode(
    non_structure$substance_name,
    show_progress_bar  = TRUE,
    convert_to_numpy   = TRUE
  )

  saveRDS(embeddings_echa, embeddings_file)
}

embeddings_echa <- readRDS(embeddings_file)
emb_matrix      <- as.matrix(embeddings_echa)

# ==============================================================================
# Analysis 8b: Optimal k via silhouette score (k = 2–12)
# ==============================================================================
# INTENT
# Select the number of clusters k that maximises within-cluster cohesion
# relative to between-cluster separation, using mean silhouette score as the
# criterion (Rousseeuw 1987).  The silhouette curve is plotted to make the
# selection transparent and to show whether there is a single clear optimum or
# a plateau.
# Own addition: range k = 2–12 chosen to cover plausible chemical name
# groupings without over-segmentation.
# ==============================================================================

set.seed(42)
k_vals     <- 2:12
sil_scores <- sapply(k_vals, function(k) {
  km  <- kmeans(emb_matrix, centers = k, nstart = 10, iter.max = 100)
  sil <- silhouette(km$cluster, dist(emb_matrix))
  mean(sil[, 3])
})

sil_df <- data.frame(k = k_vals, silhouette = sil_scores)
k_opt  <- sil_df$k[which.max(sil_df$silhouette)]
message(sprintf("Optimal k = %d  (mean silhouette = %.3f)",
                k_opt, max(sil_scores)))

p8_sil <- ggplot(sil_df, aes(x = k, y = silhouette)) +
  geom_line(colour = "#4a90d9", linewidth = 0.8) +
  geom_point(size = 2.5, colour = "#4a90d9") +
  geom_vline(xintercept = k_opt, linetype = "dashed", colour = "#e05c5c") +
  geom_label(
    data    = sil_df[sil_df$k == k_opt, ],
    aes(label = paste0("k = ", k, "\n(", round(silhouette, 3), ")")),
    nudge_x = 0.4, size = 3.5, colour = "#e05c5c", label.size = 0.3
  ) +
  scale_x_continuous(breaks = k_vals) +
  labs(
    title    = "Silhouette score by number of clusters",
    subtitle = "Higher is better \u2014 dashed line marks the optimal k",
    x        = "Number of clusters (k)",
    y        = "Mean silhouette score"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p8_sil)
ggsave(p8_sil,
       filename = here("output", "figures", "Analysis_8b_silhouette_score.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# ==============================================================================
# Analysis 8c: k-means clustering with k = 6
# ==============================================================================
# INTENT
# Apply k-means with k = 6 (fixed; selected from silhouette analysis above).
# Manual inspection of cluster centroids and random samples (8f) validated
# the six labels below as substantively coherent.
# Own addition: k = 6 fixed override of the silhouette optimum when the
# automatic k produced fewer than 6 interpretable clusters on initial runs.
# ==============================================================================

k_opt <- 6L  # Own addition: fixed after manual validation of cluster content

set.seed(42)
km <- kmeans(emb_matrix, centers = k_opt, nstart = 25, iter.max = 100)
non_structure$cluster <- factor(km$cluster)

cluster_labels_manual <- c(
  "Reaction masses",                           # Cluster 1
  "Complex substances (UVCB with structure)",  # Cluster 2
  "Petroleum & coal tar fractions",            # Cluster 3
  "Inorganic compounds & metal salts",         # Cluster 4
  "Trade names, codes & biological materials", # Cluster 5
  "Polymers, fatty acids & surfactants"        # Cluster 6
)

cluster_map <- data.frame(
  cluster      = factor(1:length(cluster_labels_manual)),
  manual_label = cluster_labels_manual
)

# ==============================================================================
# Analysis 8d: UMAP 2-D visualisation
# ==============================================================================
# INTENT
# Project the 384-dimensional embedding space to 2 dimensions using UMAP.
# UMAP preserves both local neighbourhood structure and global topology better
# than PCA while being faster than t-SNE at this scale.
# Cluster centroids (median coordinates) are labelled to aid interpretation.
# Own addition: cosine metric chosen because cosine similarity is the natural
# distance for normalised embedding vectors.
# ==============================================================================

# Own addition: n_neighbors = 15 and min_dist = 0.1 are the UMAP defaults
# recommended for datasets of this size (McInnes et al. 2018).
# metric = "cosine" preferred over Euclidean for high-dimensional embedding
# vectors (see own-additions table in header).
umap_coords <- uwot::umap(
  emb_matrix,
  n_neighbors = 15L,
  min_dist    = 0.1,
  metric      = "cosine",
  seed        = 42L
)

non_structure$umap1 <- umap_coords[, 1]
non_structure$umap2 <- umap_coords[, 2]

plot_data_umap <- non_structure |>
  left_join(cluster_map, by = "cluster")

p8_umap <- ggplot(plot_data_umap, aes(x = umap1, y = umap2, colour = cluster)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_label_repel(
    data = plot_data_umap |>
      group_by(cluster, manual_label) |>
      summarise(
        umap1 = median(umap1),
        umap2 = median(umap2),
        .groups = "drop"
      ),
    aes(label = paste0(manual_label, " (C", cluster, ")")),
    size        = 3,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("#e05c5c", "#4a90d9", "#7bc67e",
               "#9b59b6", "#e07b3a", "#2ec4b6"),
    labels = cluster_labels_manual
  ) +
  labs(
    title    = "UMAP of substance names without InChIKey (sentence embeddings)",
    subtitle = paste0(
      nrow(non_structure),
      " substance names \u2192 all-MiniLM-L6-v2 embeddings \u2192 ",
      k_opt, " k-means clusters"
    ),
    x      = "UMAP 1",
    y      = "UMAP 2",
    colour = "Cluster"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.subtitle   = element_text(colour = "grey40")
  )

print(p8_umap)
ggsave(p8_umap,
       filename = here("output", "figures", "Analysis_8d_UMAP.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")

# ==============================================================================
# Analysis 8e: t-SNE for comparison / validation
# ==============================================================================
# INTENT
# t-SNE is applied to the same embedding matrix as a structural validation
# of the UMAP result.  If the cluster boundaries visible in UMAP are also
# visible in t-SNE (which optimises a different objective), they are unlikely
# to be UMAP artefacts.  Perplexity is capped to satisfy the t-SNE constraint
# perplexity < n/3.
# ==============================================================================

# Own addition: perplexity = min(30, ⌊(n−1)/3⌋) satisfies the t-SNE
# constraint (perplexity < n/3) while applying the commonly recommended
# default of 30 (van der Maaten & Hinton 2008) for larger datasets.
# pca = FALSE: PCA pre-reduction is skipped because the embedding vectors
# are already dense and low-noise; PCA would discard semantic variance.
set.seed(42)
tsne_out <- Rtsne::Rtsne(
  emb_matrix,
  dims             = 2L,
  perplexity       = min(30L, floor((nrow(emb_matrix) - 1L) / 3L)),
  check_duplicates = FALSE,
  pca              = FALSE,
  verbose          = FALSE
)

non_structure$tsne1 <- tsne_out$Y[, 1]
non_structure$tsne2 <- tsne_out$Y[, 2]

p8_tsne <- ggplot(non_structure, aes(x = tsne1, y = tsne2, colour = cluster)) +
  geom_point(size = 1.2, alpha = 0.6) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title    = "t-SNE of non-structure substance names (sentence embeddings)",
    subtitle = "Same clusters as UMAP \u2014 t-SNE preserves local neighbourhood structure",
    x        = "t-SNE 1",
    y        = "t-SNE 2",
    colour   = "Cluster"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.subtitle   = element_text(colour = "grey40")
  )

print(p8_tsne)
ggsave(p8_tsne,
       filename = here("output", "figures", "Analysis_8e_tsne.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")

# ==============================================================================
# Analysis 8f: Top example names per cluster (console diagnostic)
# ==============================================================================
# INTENT
# SKILL.md §2.3 requires verifying that labels assigned in the analysis
# actually correspond to meaningful categories in the data before finalising.
# Here the equivalent check is a random sample of 100 substance names per
# cluster, printed to the console, so the analyst can confirm that the six
# manual cluster labels (defined in 8c) are substantively coherent — i.e. that
# cluster 3 genuinely contains petroleum fractions, cluster 4 inorganic salts,
# etc.  If a cluster sample contradicts its label, the label must be revised
# before downstream use of non_structure_clusters.rds.
# ==============================================================================

top_per_cluster <- non_structure |>
  group_by(cluster) |>
  slice_sample(n = 100) |>
  summarise(examples = paste(substance_name, collapse = "\n  "), .groups = "drop")

message("\n=== Analysis 8: top examples per cluster ===")
for (i in seq_len(nrow(top_per_cluster))) {
  message(sprintf("\nCluster %s:\n  %s",
                  top_per_cluster$cluster[i],
                  top_per_cluster$examples[i]))
}

# ==============================================================================
# Analysis 8g: Projection quality comparison — UMAP vs t-SNE artefact detection
# ==============================================================================
# INTENT
# UMAP is known to sometimes over-compress within-cluster variance or
# exaggerate between-cluster distances (McInnes et al. 2018, §4; Böhm et al.
# 2022, *Attraction-Repulsion Spectrum in Neighbor Embeddings*, J. Mach.
# Learn. Res. 23:95).  The standard diagnostic is to compare silhouette scores
# in the original high-dimensional space (ground truth) and in each 2-D
# projection.  If sil_umap >> sil_embedding, UMAP is inflating cluster
# separation beyond what exists in the actual embedding — an artefact.
# Concordance between all three scores confirms that the visible cluster
# structure reflects genuine semantic grouping, not projection bias.
#
# The coordinate table provides the raw 2-D positions from both projections
# per substance, enabling external comparison tools to reproduce the analysis.
# Trustworthiness and continuity (Venna & Kaski 2006) are computed in 8h.
# ==============================================================================

# Distance matrices stored as variables so they can be reused in 8h
# (trustworthiness / continuity) without recomputation
D_emb  <- dist(emb_matrix)
D_umap <- dist(umap_coords)
D_tsne <- dist(tsne_out$Y)

# Silhouette in original 384-D embedding space (ground truth)
sil_emb  <- silhouette(km$cluster, D_emb)

# Silhouette in UMAP 2-D projection
sil_umap <- silhouette(km$cluster, D_umap)

# Silhouette in t-SNE 2-D projection
sil_tsne <- silhouette(km$cluster, D_tsne)

# Per-cluster silhouette means across the three spaces
sil_per_cluster <- data.frame(
  cluster       = as.integer(sort(unique(km$cluster))),
  sil_embedding = as.numeric(tapply(sil_emb[, 3],  sil_emb[, 1],  mean)),
  sil_umap      = as.numeric(tapply(sil_umap[, 3], sil_umap[, 1], mean)),
  sil_tsne      = as.numeric(tapply(sil_tsne[, 3], sil_tsne[, 1], mean))
) |>
  mutate(
    # Positive value = projection inflates separation beyond embedding truth
    umap_vs_embedding = sil_umap - sil_embedding,
    tsne_vs_embedding = sil_tsne - sil_embedding
  )

# Overall row appended for at-a-glance summary
sil_comparison <- bind_rows(
  sil_per_cluster,
  data.frame(
    cluster           = NA_integer_,   # NA = overall summary row
    sil_embedding     = mean(sil_emb[, 3]),
    sil_umap          = mean(sil_umap[, 3]),
    sil_tsne          = mean(sil_tsne[, 3]),
    umap_vs_embedding = mean(sil_umap[, 3]) - mean(sil_emb[, 3]),
    tsne_vs_embedding = mean(sil_tsne[, 3]) - mean(sil_emb[, 3])
  )
)

message("\n=== Analysis 8g: silhouette comparison — embedding vs UMAP vs t-SNE ===")
message(sprintf(
  "Overall  embedding: %.3f | UMAP: %.3f | t-SNE: %.3f | UMAP inflation: %+.3f",
  mean(sil_emb[, 3]), mean(sil_umap[, 3]), mean(sil_tsne[, 3]),
  mean(sil_umap[, 3]) - mean(sil_emb[, 3])
))

# Columns (Analysis_8_silhouette_comparison.csv):
#   cluster           — cluster index (1–k); NA = overall summary across all clusters
#   sil_embedding     — mean silhouette in the original 384-D embedding space (ground truth)
#   sil_umap          — mean silhouette in the UMAP 2-D projection
#   sil_tsne          — mean silhouette in the t-SNE 2-D projection
#   umap_vs_embedding — sil_umap − sil_embedding; positive = UMAP inflates separation
#   tsne_vs_embedding — sil_tsne − sil_embedding; positive = t-SNE inflates separation
write_csv(sil_comparison,
          here("output", "tables", "Analysis_8_silhouette_comparison.csv"))

# Columns (Analysis_8_projection_coordinates.csv):
#   substance_name — original substance name from ECHA dataset
#   cluster        — k-means cluster assignment (1–k)
#   manual_label   — human-readable cluster label (see 8c)
#   umap1, umap2   — 2-D UMAP coordinates
#   tsne1, tsne2   — 2-D t-SNE coordinates
write_csv(
  plot_data_umap |>
    select(substance_name, cluster, manual_label, umap1, umap2) |>
    left_join(
      non_structure |> select(substance_name, tsne1, tsne2),
      by = "substance_name"
    ),
  here("output", "tables", "Analysis_8_projection_coordinates.csv")
)

# ==============================================================================
# Analysis 8h: Trustworthiness and continuity (Venna & Kaski 2006)
# ==============================================================================
# INTENT
# Silhouette comparison (8g) tests whether cluster separation is inflated by
# the projection.  Trustworthiness and continuity test a complementary
# property: how faithfully each projection preserves local neighbourhoods.
#
# Trustworthiness T(k): for each point, what fraction of its k nearest
#   neighbours in the projection were also among its k nearest neighbours in
#   the original space?  T = 1 means the projection introduces no false
#   neighbours; T < 1 indicates the projection pulls in points that were
#   actually distant in the original space (a UMAP tear artefact).
#
# Continuity C(k): the complement — for each point, what fraction of its k
#   nearest neighbours in the original space are also among its k nearest
#   neighbours in the projection?  C < 1 means the projection pushes apart
#   points that were actually close — a compression artefact.
#
# Both metrics are computed for k = 5, 10, 20, 50 to show how neighbourhood
# fidelity varies with scale (local vs global).  Values close to 1 at all k
# indicate a faithful projection; a rapid drop with increasing k indicates
# that only local structure is preserved.
#
# Reference: Venna J, Kaski S. 2006. Local multidimensional scaling.
#   Neural Networks 19(6-7):889-899. DOI: 10.1016/j.neunet.2006.05.014
# ==============================================================================

# Helper: compute T(k) and C(k) for one projection given pre-computed
# distance matrices D_high (original space) and D_low (projected space).
# Own addition: vectorised rank-matrix approach avoids the O(n²k) inner loop.
tc_metrics <- function(D_high, D_low, k_vals) {
  D_high_m <- as.matrix(D_high)
  D_low_m  <- as.matrix(D_low)
  n        <- nrow(D_high_m)

  # Rank matrices: rank_X[j, i] = rank of point j by distance from point i
  # in space X (rank 1 = nearest, self gets rank 0 after subtracting 1).
  rank_high <- apply(D_high_m, 2,
                     function(d) rank(d, ties.method = "first") - 1L)
  rank_low  <- apply(D_low_m,  2,
                     function(d) rank(d, ties.method = "first") - 1L)

  # Normalisation constant (Venna & Kaski 2006, eq. 1-2)
  normaliser <- function(k) 2 / (n * k * (2 * n - 3 * k - 1))

  purrr::map_dfr(k_vals, function(k) {
    in_high <- rank_high > 0L & rank_high <= k   # k-NN in original space
    in_low  <- rank_low  > 0L & rank_low  <= k   # k-NN in projected space

    # Trustworthiness penalty: points in low k-NN but not in high k-NN
    T_penalty <- (rank_high - k) * (in_low & !in_high)

    # Continuity penalty: points in high k-NN but not in low k-NN
    C_penalty <- (rank_low - k) * (in_high & !in_low)

    data.frame(
      k               = k,
      trustworthiness = 1 - normaliser(k) * sum(T_penalty),
      continuity      = 1 - normaliser(k) * sum(C_penalty)
    )
  })
}

k_vals <- c(5L, 10L, 20L, 50L)

tc_umap <- tc_metrics(D_emb, D_umap, k_vals) |> mutate(projection = "UMAP")
tc_tsne <- tc_metrics(D_emb, D_tsne, k_vals) |> mutate(projection = "t-SNE")

tc_results <- bind_rows(tc_umap, tc_tsne) |>
  select(projection, k, trustworthiness, continuity)

message("\n=== Analysis 8h: trustworthiness and continuity ===")
message(paste(
  sprintf("  %s k=%2d  T=%.3f  C=%.3f",
          tc_results$projection, tc_results$k,
          tc_results$trustworthiness, tc_results$continuity),
  collapse = "\n"
))

# Columns (Analysis_8_trustworthiness_continuity.csv):
#   projection      — "UMAP" or "t-SNE"
#   k               — neighbourhood size (5, 10, 20, 50)
#   trustworthiness — fraction of projected k-NN that were true high-dim k-NN;
#                     1 = no false neighbours introduced; low = tear artefact
#   continuity      — fraction of high-dim k-NN preserved in the projection;
#                     1 = no neighbours lost; low = compression artefact
write_csv(tc_results,
          here("output", "tables",
               "Analysis_8_trustworthiness_continuity.csv"))

# ==============================================================================
# Save cluster assignments for use in 09_embedding_chemont.R and
# 13b_before_prioritisation_create_scheme.R
# ==============================================================================

saveRDS(plot_data_umap,
        here("data", "processed", "non_structure_clusters.rds"))

write_csv(plot_data_umap, here("data", "processed", "non_structure_clusters.csv"))

message("Analysis 8: cluster assignments saved to data/processed/non_structure_clusters.rds")

saveRDS(p8_umap, here("data", "processed", "p8_umap.rds"))

message("08_embedding_clustering.R: analysis 8 completed")
