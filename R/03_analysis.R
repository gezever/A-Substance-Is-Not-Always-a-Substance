# ==============================================================================
# 02_analysis.R
# Core analysis: "Regulatory substances are not interoperable entities"
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(igraph)
library(ggraph)
library(UpSetR)
library(scales)
library(reticulate)
library(uwot)
library(Rtsne)
library(ggrepel)
library(stringr)

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 1: What is a "substance"?
# Structure-defined vs non-structure entities
# ==============================================================================

df1 <- all_substances |>
  mutate(has_inchikey = !is.na(inchikey)) |>
  count(has_inchikey) |>
  mutate(
    label = if_else(has_inchikey, "Structure-defined\n(InChIKey present)", "Non-structure entity\n(InChIKey absent)"),
    pct   = n / sum(n)
  )

p1 <- ggplot(df1, aes(x = label, y = n, fill = has_inchikey)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(n, "\n(", percent(pct, accuracy = 0.1), ")")),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("FALSE" = "#e05c5c", "TRUE" = "#4a90d9")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Analysis 1: What is a 'substance'?",
    subtitle = "A large proportion of regulatory entries are not chemically identifiable substances",
    x        = NULL,
    y        = "Number of records"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p1)


# ==============================================================================
# Analysis 2: Per source — interoperability differs by regulation
# ==============================================================================

df2 <- all_substances |>
  group_by(source) |>
  summarise(
    with_inchikey    = sum(!is.na(inchikey)),
    without_inchikey = sum(is.na(inchikey)),
    total            = n(),
    .groups = "drop"
  ) |>
  mutate(pct_with = with_inchikey / total) |>
  arrange(desc(pct_with)) |>
  pivot_longer(c(with_inchikey, without_inchikey),
               names_to = "type", values_to = "n") |>
  mutate(
    type   = recode(type,
                    with_inchikey    = "InChIKey present",
                    without_inchikey = "InChIKey absent"),
    source = factor(source, levels = unique(source[type == "InChIKey present"]))
  )

p2 <- ggplot(df2, aes(x = source, y = n, fill = type)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c("InChIKey present" = "#4a90d9", "InChIKey absent" = "#e05c5c")) +
  coord_flip() +
  labs(
    title    = "Analysis 2: Interoperability by source",
    subtitle = "Regulation determines whether substances are chemically identifiable",
    x        = NULL,
    y        = "Share of records",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    plot.subtitle    = element_text(colour = "grey40")
  )

print(p2)


# ==============================================================================
# Analysis 3: CAS ↔ InChIKey inconsistencies
# ==============================================================================

# One CAS → multiple InChIKeys
cas_to_inchi <- all_substances |>
  filter(!is.na(cas_number), !is.na(inchikey)) |>
  group_by(cas_number) |>
  summarise(n_inchikeys = n_distinct(inchikey), .groups = "drop")

p3a <- ggplot(cas_to_inchi, aes(x = n_inchikeys)) +
  geom_histogram(binwidth = 1, fill = "#4a90d9", colour = "white") +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Analysis 3a: CAS → InChIKey (one CAS, how many InChIKeys?)",
    subtitle = "A CAS number should uniquely refer to a single structure",
    x        = "Number of distinct InChIKeys per CAS",
    y        = "Number of CAS numbers"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p3a)

# One InChIKey → multiple CAS numbers
inchi_to_cas <- all_substances |>
  filter(!is.na(cas_number), !is.na(inchikey)) |>
  group_by(inchikey) |>
  summarise(n_cas = n_distinct(cas_number), .groups = "drop")

p3b <- ggplot(inchi_to_cas, aes(x = n_cas)) +
  geom_histogram(binwidth = 1, fill = "#e07b3a", colour = "white") +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Analysis 3b: InChIKey → CAS (one structure, how many CAS numbers?)",
    subtitle = "The same chemical structure can carry multiple CAS numbers",
    x        = "Number of distinct CAS numbers per InChIKey",
    y        = "Number of InChIKeys"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p3b)

message(sprintf(
  "Analysis 3: %d CAS numbers map to >1 InChIKey; %d InChIKeys have >1 CAS",
  sum(cas_to_inchi$n_inchikeys > 1),
  sum(inchi_to_cas$n_cas > 1)
))


# ==============================================================================
# Analysis 4: Entity type classification
# Distinguishing molecules, substance groups, analytical parameters,
# mixtures, and regulatory/administrative entries
# ==============================================================================

all_substances <- all_substances |>
  mutate(
    entity_type = case_when(
      
      # 0. Individual molecule — structure-defined
      !is.na(inchikey) ~ "Molecule",
      
      # 1. Regulatory / administrative entry (most specific — check first)
      grepl(
        paste(
          "^entry\\s+\\d+$",
          "see group members",
          "substances which (are|is) classified",
          "regulation \\(ec\\)",
          "directive \\d{4}",
          "annex [ivx]+",
          "appendix \\d",
          sep = "|"
        ),
        substance_name, ignore.case = TRUE
      ) ~ "Regulatory entry",
      
      # 2. Analytical parameter / sum of analytes
      grepl("^reaction (mass|product) of|\\bC\\d+[-/]C\\d+\\b",
            substance_name, ignore.case = TRUE) |
        grepl("^[^;]+;[^;]+", substance_name) ~ "Analytical parameter",
      
      # 3. Mixture (trade names, complex compositions, petroleum products)
      grepl(
        paste(
          "trade name",
          "^a? ?mixture (of|composed)",
          "\\[a complex combination",
          "creosote",
          "petrolatum",
          "\\bnaphtha\\b",
          "slack wax",
          "kerosine",
          sep = "|"
        ),
        substance_name, ignore.case = TRUE
      ) ~ "Mixture",
      
      # 4. Substance group (chemical class, element + compounds, acronym groups)
      grepl(
        paste(
          "and its (compounds?|salts?|esters?|isomers?)",
          "\\(group(s)?\\)",
          "compounds?$",
          "fibres?$",
          "\\b(PAH|PCB|PFAS|PFC|PCDD|PBDE|PBB|PCT|HBCDD)\\b",
          sep = "|"
        ),
        substance_name, ignore.case = TRUE
      ) ~ "Substance group",
      
      # 5. Remainder: has a CAS but no resolvable structure
      !is.na(cas_number) ~ "CAS without structure",
      
      # 6. No identifier at all
      TRUE ~ "Unclassified"
    )
  )

# ------------------------------------------------------------------------------
# 4a: Overall distribution of entity types
# ------------------------------------------------------------------------------

df4 <- all_substances |>
  count(entity_type) |>
  mutate(pct = n / sum(n)) |>
  arrange(desc(n))

entity_colours <- c(
  "Molecule"              = "#4a90d9",
  "Substance group"       = "#7bc67e",
  "Analytical parameter"  = "#f5a623",
  "Mixture"               = "#9b59b6",
  "Regulatory entry"      = "#e05c5c",
  "CAS without structure" = "#aaaaaa",
  "Unclassified"          = "#dddddd"
)

p4a <- ggplot(df4, aes(x = reorder(entity_type, n), y = n, fill = entity_type)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(n, " (", percent(pct, accuracy = 0.1), ")")),
            hjust = -0.05, size = 3.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)), labels = comma) +
  scale_fill_manual(values = entity_colours) +
  coord_flip() +
  labs(
    title    = "Analysis 4a: What kind of entity is a regulatory 'substance'?",
    subtitle = "Most non-molecule entries are groups, analytical parameters, mixtures, or administrative placeholders",
    x        = NULL,
    y        = "Number of records"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p4a)

# ------------------------------------------------------------------------------
# 4b: Entity type breakdown per source
# ------------------------------------------------------------------------------

df4b <- all_substances |>
  count(source, entity_type) |>
  group_by(source) |>
  mutate(pct = n / sum(n)) |>
  ungroup() |>
  mutate(
    entity_type = factor(entity_type, levels = names(entity_colours)),
    source      = reorder(source, -n, sum)
  )

p4b <- ggplot(df4b, aes(x = source, y = pct, fill = entity_type)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = entity_colours) +
  coord_flip() +
  labs(
    title    = "Analysis 4b: Entity type composition per source",
    subtitle = "Each regulatory list has a different mix of molecules, groups, and administrative entries",
    x        = NULL,
    y        = "Share of records",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.subtitle   = element_text(colour = "grey40")
  ) +
  guides(fill = guide_legend(nrow = 2))

print(p4b)


# ==============================================================================
# Analysis 5: Overlap between lists (UpSet plot)
# ==============================================================================

# Use InChIKey as linking key; restrict to sources with at least 10 unique InChIKeys
inchi_per_source <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(source, inchikey)

sources_sufficient <- inchi_per_source |>
  count(source) |>
  filter(n >= 10) |>
  pull(source)

upset_input <- inchi_per_source |>
  filter(source %in% sources_sufficient) |>
  mutate(present = 1L) |>
  pivot_wider(names_from = source, values_from = present, values_fill = 0L) |>
  select(-inchikey) |>
  as.data.frame()

upset(
  upset_input,
  nsets       = ncol(upset_input),
  nintersects = 30,
  order.by    = "freq",
  mb.ratio    = c(0.55, 0.45),
  main.bar.color = "#4a90d9",
  sets.bar.color = "#e07b3a",
  text.scale  = c(1.2, 1.1, 1, 1, 1.1, 1),
  mainbar.y.label = "Overlap (number of shared InChIKeys)",
  sets.x.label    = "InChIKeys per list"
)
title("Analysis 5: Overlap of structure-defined substances across lists",
      line = 3, cex.main = 1)


# ==============================================================================
# Analysis 6: Coverage of structure-based linking
# ==============================================================================

df6 <- all_substances |>
  group_by(source) |>
  summarise(
    n_records          = n(),
    n_unique_inchikeys = n_distinct(inchikey[!is.na(inchikey)]),
    .groups = "drop"
  )

p6 <- ggplot(df6, aes(x = n_records, y = n_unique_inchikeys, label = source)) +
  geom_point(colour = "#4a90d9", size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_text(vjust = -0.7, hjust = 0.5, size = 3, colour = "grey30") +
  scale_x_continuous(labels = comma, trans = "log10") +
  scale_y_continuous(labels = comma, trans = "log10") +
  labs(
    title    = "Analysis 6: Coverage of structure-based linking per source",
    subtitle = "Points below the diagonal = more records than unique structures (duplication / ambiguity)",
    x        = "Number of records (log)",
    y        = "Number of unique InChIKeys (log)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p6)


# ==============================================================================
# Analysis 7: Network visualisation — bipartite graph (substance ↔ list)
# ==============================================================================

# Restrict to substances present in at least 2 lists (for readability)
inchi_multi <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(inchikey, source) |>
  group_by(inchikey) |>
  filter(n() >= 2) |>
  ungroup()

# Build edge list
edges <- inchi_multi |>
  select(from = inchikey, to = source)

# Build node list with type
node_substances <- data.frame(
  name = unique(edges$from),
  type = "substance",
  stringsAsFactors = FALSE
)
node_sources <- data.frame(
  name = unique(edges$to),
  type = "list",
  stringsAsFactors = FALSE
)
nodes <- bind_rows(node_substances, node_sources)

g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

set.seed(42)
p7 <- ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.15, colour = "grey60") +
  geom_node_point(aes(colour = type, size = type)) +
  geom_node_text(
    aes(label = if_else(type == "list", name, NA_character_)),
    repel = TRUE, size = 3.5, colour = "#333333"
  ) +
  scale_colour_manual(values = c("substance" = "#4a90d9", "list" = "#e05c5c")) +
  scale_size_manual(values  = c("substance" = 1.5,        "list" = 5)) +
  labs(
    title    = "Analysis 7: Bipartite network substance ↔ regulatory list",
    subtitle = paste0(
      "Only substances present in \u2265 2 lists shown (n = ",
      nrow(node_substances), " substances, ", nrow(node_sources), " lists)"
    ),
    colour   = NULL,
    size     = NULL
  ) +
  theme_graph(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.subtitle   = element_text(colour = "grey40")
  )

print(p7)


# ==============================================================================
# Analysis 8: Embedding + clustering of "Unclassified" substance names
# Uses sentence-transformers (Python) via reticulate → UMAP + k-means
# Requirements:
#   pip install sentence-transformers
#   install.packages(c("reticulate", "uwot", "Rtsne", "ggrepel"))
# ==============================================================================



# ------------------------------------------------------------------------------
# 8a: Embed substance names with sentence-transformers
# ------------------------------------------------------------------------------

embeddings_file <- here("data", "processed", "embeddings.rds")

if (!file.exists(embeddings_file)) {
  
  unclassified <- all_substances |>
    filter(entity_type == "Unclassified") |>
    distinct(substance_name) |>
    filter(!is.na(substance_name), nzchar(substance_name))
  
  message(sprintf("Analysis 8: embedding %d unclassified substance names", nrow(unclassified)))
  
  st <- reticulate::import("sentence_transformers")
  model <- st$SentenceTransformer("all-MiniLM-L6-v2")
  
  embeddings <- model$encode(
    unclassified$substance_name,
    show_progress_bar = TRUE,
    convert_to_numpy   = TRUE
  )
  
  saveRDS(embeddings, embeddings_file)
  
  
} else { 
  # load
  embeddings <- readRDS(embeddings_file) 
}

# embeddings is an (n × 384) numpy matrix
emb_matrix <- as.matrix(embeddings)

# ------------------------------------------------------------------------------
# 8b: k-means clustering on the embedding space
# ------------------------------------------------------------------------------

set.seed(42)
k <- 8L   # adjust after inspecting the UMAP plot

km <- kmeans(emb_matrix, centers = k, nstart = 25, iter.max = 100)
unclassified$cluster <- factor(km$cluster)

# ------------------------------------------------------------------------------
# 8c: UMAP for 2-D visualisation
# ------------------------------------------------------------------------------

umap_coords <- uwot::umap(
  emb_matrix,
  n_neighbors = 15L,
  min_dist    = 0.1,
  metric      = "cosine",
  seed        = 42L
)

unclassified$umap1 <- umap_coords[, 1]
unclassified$umap2 <- umap_coords[, 2]

# Per-cluster label: the substance name closest to each cluster centroid
cluster_labels <- unclassified |>
  group_by(cluster) |>
  slice_sample(n = 1) |>          # replace with centroid-nearest if preferred
  ungroup() |>
  select(cluster, label = substance_name)

plot_data_umap <- unclassified |>
  left_join(cluster_labels, by = "cluster")

p8_umap <- ggplot(plot_data_umap, aes(x = umap1, y = umap2, colour = cluster)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_label_repel(
    data          = plot_data_umap |> group_by(cluster) |>
      summarise(umap1 = median(umap1), umap2 = median(umap2),
                label = first(label), .groups = "drop"),
    aes(label = paste0("C", cluster, ": ", str_trunc(label, 40))),
    size          = 3,
    label.padding = unit(0.2, "lines"),
    show.legend   = FALSE
  ) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title    = "Analysis 8a: UMAP of unclassified substance names (sentence embeddings)",
    subtitle = paste0(
      nrow(unclassified), " substance names \u2192 all-MiniLM-L6-v2 embeddings \u2192 ",
      k, " k-means clusters"
    ),
    x        = "UMAP 1",
    y        = "UMAP 2",
    colour   = "Cluster"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.subtitle   = element_text(colour = "grey40")
  )

print(p8_umap)

# ------------------------------------------------------------------------------
# 8d: t-SNE for comparison / validation
# ------------------------------------------------------------------------------

set.seed(42)
tsne_out <- Rtsne::Rtsne(
  emb_matrix,
  dims         = 2L,
  perplexity   = min(30L, floor((nrow(emb_matrix) - 1L) / 3L)),
  check_duplicates = FALSE,
  pca          = FALSE,
  verbose      = FALSE
)

unclassified$tsne1 <- tsne_out$Y[, 1]
unclassified$tsne2 <- tsne_out$Y[, 2]

p8_tsne <- ggplot(unclassified, aes(x = tsne1, y = tsne2, colour = cluster)) +
  geom_point(size = 1.2, alpha = 0.6) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title    = "Analysis 8b: t-SNE of unclassified substance names (sentence embeddings)",
    subtitle = "Same clusters as UMAP — t-SNE preserves local neighbourhood structure",
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

# ------------------------------------------------------------------------------
# 8e: Top terms per cluster (most central substance names)
# ------------------------------------------------------------------------------

top_per_cluster <- unclassified |>
  group_by(cluster) |>
  slice_sample(n = 10) |>
  summarise(examples = paste(substance_name, collapse = "\n  "), .groups = "drop")

message("\n=== Analysis 8: top examples per cluster ===")
for (i in seq_len(nrow(top_per_cluster))) {
  message(sprintf("\nCluster %s:\n  %s",
                  top_per_cluster$cluster[i],
                  top_per_cluster$examples[i]))
}

message("02_analysis.R completed")
