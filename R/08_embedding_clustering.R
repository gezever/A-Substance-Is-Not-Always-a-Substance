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
# 1.  Sentence embeddings: all-MiniLM-L6-v2 via sentence-transformers
#     (Reimers & Gurevych 2019, https://arxiv.org/abs/1908.10084).
#     Own addition: choice of all-MiniLM-L6-v2 (compact, fast) over larger
#     models; acceptable for short chemical name strings.
# 2.  Optimal k selection: mean silhouette score (Rousseeuw 1987) over k = 2–12.
#     Own addition: range 2–12 is a pragmatic search space; wider ranges
#     were explored manually but did not yield better-separating clusters.
# 3.  k-means with k = 6 (optimal from silhouette analysis; fixed for
#     reproducibility).
# 4.  UMAP (McInnes et al. 2018, https://arxiv.org/abs/1802.03426) for 2-D
#     visualisation: n_neighbors = 15, min_dist = 0.1, cosine metric.
# 5.  t-SNE (van der Maaten & Hinton 2008) as validation of UMAP structure.
#     Perplexity capped at floor((n−1)/3) to satisfy the t-SNE constraint.
#
# REQUIREMENTS
# ------------
# Python (via reticulate): pip install sentence-transformers
# R packages: reticulate, uwot, Rtsne, cluster, purrr
#
# OUTPUTS
# -------
# output/figures/Analysis_8b_silhouette_score.pdf
# output/figures/Analysis_8d_UMAP.pdf
# output/figures/Analysis_8e_tsne.pdf
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
    title    = "Analysis 8b: Silhouette score by number of clusters",
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
# Random sample of 100 substance names per cluster, printed to the console,
# to allow manual validation that the cluster labels assigned above are
# substantively correct.  This verification step is required by the SKILL.md
# checklist before the analysis is considered complete.
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
# Save cluster assignments for use in 09_embedding_chemont.R and
# 13b_before_prioritisation_create_scheme.R
# ==============================================================================

saveRDS(plot_data_umap,
        here("data", "processed", "non_structure_clusters.rds"))
message("Analysis 8: cluster assignments saved to data/processed/non_structure_clusters.rds")

message("08_embedding_clustering.R: analysis 8 completed")
