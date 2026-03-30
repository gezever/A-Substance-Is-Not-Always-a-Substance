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
library(cluster)
library(purrr)
library(httr)
library(jsonlite)

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
ggsave(p1, 
       filename = here("output", "figures", "Analysis_1_What_is_a_substance.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

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
    title    = "Interoperability by source",
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
ggsave(p2, 
       filename = here("output", "figures", "Analysis_2_interoperability_differs_by_regulation.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")


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
  scale_y_log10(labels = comma) +
  labs(
    title    = "CAS → InChIKey (one CAS, how many InChIKeys?)",
    subtitle = "A CAS number should uniquely refer to a single structure",
    x        = "Number of distinct InChIKeys per CAS",
    y        = "Number of CAS numbers (log\u2081\u2080)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p3a)
ggsave(p3a, 
       filename = here("output", "figures", "Analysis_3a_CAS↔InChIKey_inconsistencies_multiple_INCHIKEYS.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# One InChIKey → multiple CAS numbers
inchi_to_cas <- all_substances |>
  filter(!is.na(cas_number), !is.na(inchikey)) |>
  group_by(inchikey) |>
  summarise(n_cas = n_distinct(cas_number), .groups = "drop")

p3b <- ggplot(inchi_to_cas, aes(x = n_cas)) +
  geom_histogram(binwidth = 1, fill = "#e07b3a", colour = "white") +
  scale_y_log10(labels = comma) +
  labs(
    title    = "InChIKey → CAS (one structure, how many CAS numbers?)",
    subtitle = "The same chemical structure can carry multiple CAS numbers",
    x        = "Number of distinct CAS numbers per InChIKey",
    y        = "Number of InChIKeys (log\u2081\u2080)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p3b)
ggsave(p3b, 
       filename = here("output", "figures", "Analysis_3b_CAS↔InChIKey_inconsistencies_multiple_CAS_numbers.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

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
    title    = "What kind of entity is a regulatory 'substance'?",
    subtitle = "Most non-molecule entries are groups, analytical parameters, mixtures, or administrative placeholders",
    x        = NULL,
    y        = "Number of records"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p4a)
ggsave(p4a, 
       filename = here("output", "figures", "Analysis_4a_Overall_distribution_of_entity_types.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

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
    title    = "Entity type composition per source",
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
ggsave(p4b, 
       filename = here("output", "figures", "Analysis_4b_Entity_type_breakdown_per_source.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

# ==============================================================================
# Analysis 4c: Parent compound extraction for Substance groups
# Goal: for "X and its compounds" entries, extract the parent compound and
#       look up an InChIKey via PubChem → linkable via ClassyFire
# ==============================================================================

extract_parent <- function(name) {
  m <- regmatches(name, regexpr(
    "^(.+?)\\s+(and its\\b|,\\s*(?:compounds?|salts?|esters?|isomers?)\\s*$)",
    name, perl = TRUE, ignore.case = TRUE
  ))
  if (length(m) == 0) return(NA_character_)
  trimws(gsub(
    "\\s+(and its\\b|,\\s*(?:compounds?|salts?|esters?|isomers?)\\s*)$",
    "", m, perl = TRUE, ignore.case = TRUE
  ))
}

lookup_parent_inchikey <- function(name) {
  if (is.na(name) || !nzchar(trimws(name))) return(NA_character_)

  cache_file <- file.path(here("data", "cache", "inchikey"), paste0(name, ".json"))

  if (file.exists(cache_file)) {
    json <- fromJSON(paste(readLines(cache_file, warn = FALSE), collapse = "\n"))
    return(json$PropertyTable$Properties$InChIKey)
  }

  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
    URLencode(name, reserved = TRUE),
    "/property/InChIKey/JSON"
  )

  res <- GET(url)
  Sys.sleep(0.2)

  if (status_code(res) == 200) {
    raw <- content(res, "text", encoding = "UTF-8")
    writeLines(raw, cache_file)
    json <- fromJSON(raw)
    return(json$PropertyTable$Properties$InChIKey)
  } else {
    return(NA_character_)
  }
}

substance_groups <- all_substances |>
  filter(entity_type == "Substance group") |>
  distinct(substance_name) |>
  mutate(
    base_compound   = sapply(substance_name, extract_parent),
    base_resolvable = !is.na(base_compound)
  )

substance_groups <- substance_groups |>
  mutate(
    base_inchikey = sapply(base_compound, lookup_parent_inchikey)
  )

# ------------------------------------------------------------------------------
# 4d: Linkability taxonomy — general overview of the complete ECHA dataset
# ------------------------------------------------------------------------------

linkability <- all_substances |>
  distinct(substance_name, entity_type, inchikey) |>
  left_join(
    substance_groups |> select(substance_name, base_inchikey),
    by = "substance_name"
  ) |>
  mutate(
    linkability = case_when(
      !is.na(inchikey)                                        ~ "Structure (InChIKey)",
      entity_type == "Substance group" & !is.na(base_inchikey) ~ "Group — base compound resolvable",
      entity_type == "Substance group"                        ~ "Group — base compound not resolvable",
      entity_type == "CAS without structure"                    ~ "CAS without structure",
      entity_type %in% c("Analytical parameter", "Mixture")    ~ "UVCB / Mixture",
      entity_type == "Regulatory entry"                        ~ "Regulatory entry",
      TRUE                                                      ~ "Unidentified"
    ),
    linkability = factor(linkability, levels = c(
      "Structure (InChIKey)",
      "Group — base compound resolvable",
      "Group — base compound not resolvable",
      "CAS without structure",
      "UVCB / Mixture",
      "Regulatory entry",
      "Unidentified"
    ))
  )

saveRDS(linkability, here("data", "processed", "linkability_taxonomy.rds"))

df4d <- linkability |>
  count(linkability) |>
  mutate(pct = round(n / sum(n) * 100, 1))

print(df4d)

linkability_colours <- c(
  "Structure (InChIKey)"                  = "#4a90d9",
  "Group — base compound resolvable"      = "#7bc67e",
  "Group — base compound not resolvable"  = "#b5e5a0",
  "CAS without structure"             = "#aaaaaa",
  "UVCB / Mixture"                    = "#9b59b6",
  "Regulatory entry"                  = "#e05c5c",
  "Unidentified"                      = "#dddddd"
)

p4d <- ggplot(df4d, aes(x = reorder(linkability, n), y = n, fill = linkability)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(n, " (", pct, "%)")),
            hjust = -0.05, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)), labels = comma) +
  scale_fill_manual(values = linkability_colours) +
  coord_flip() +
  labs(
    title    = "Linkability taxonomy of ECHA substances",
    subtitle = "What proportion can be unambiguously linked to a chemical structure?",
    x        = NULL,
    y        = "Number of unique substance names"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p4d)

ggsave(p4d, 
       filename = here("output", "figures", "Analysis_4d_Linkable-taxonomie—general_overview_of_the_complete_ECHA-dataset.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

saveRDS(linkability, here("data", "processed", "linkability_taxonomy.rds"))

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

pdf(
  here("output", "figures", "Analysis_5_Overlap_between_lists_UpSet.pdf"),
  width = 14, height = 7
)
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
grid::grid.text(
  "Analysis 5: Overlap of structure-defined substances across lists",
  x = 0.65, y = 0.97,
  gp = grid::gpar(fontsize = 12, fontface = "bold")
)
dev.off()


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
    title    = "Coverage of structure-based linking per source",
    subtitle = "Points below the diagonal = more records than unique structures (duplication / ambiguity)",
    x        = "Number of records (log)",
    y        = "Number of unique InChIKeys (log)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p6)

ggsave(p6, 
       filename = here("output", "figures", "Analysis_6_Coverage_of_structure-based_linking.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

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
    title    = "Bipartite network substance ↔ regulatory list",
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

ggsave(p7, 
       filename = here("output", "figures", "Analysis_7_Network_visualisation—bipartite_graph.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")


# ==============================================================================
# Analysis 8: Embedding + clustering of "non_structure" substance names
# Uses sentence-transformers (Python) via reticulate → UMAP + k-means
# Requirements:
#   pip install sentence-transformers
#   install.packages(c("reticulate", "uwot", "Rtsne", "ggrepel"))
# ==============================================================================



# ------------------------------------------------------------------------------
# 8a: Embed substance names with sentence-transformers
# ------------------------------------------------------------------------------

embeddings_file <- here("data", "processed", "embeddings_echa.rds")


  
non_structure <- all_substances |>
  filter(is.na(inchikey)) |>
  filter(!entity_type %in% c("Regulatory entry")) |>  # not chemical entities
  distinct(substance_name) |>
  filter(!is.na(substance_name), nzchar(substance_name))
  
if (!file.exists(embeddings_file)) {
  message(sprintf("Analysis 8: embedding %d non_structure substance names", nrow(non_structure)))
  
  st <- reticulate::import("sentence_transformers")
  model <- st$SentenceTransformer("all-MiniLM-L6-v2")
  
  embeddings_echa <- model$encode(
    non_structure$substance_name,
    show_progress_bar = TRUE,
    convert_to_numpy   = TRUE
  )
  
  saveRDS(embeddings_echa, embeddings_file)
  
  
} 

# load
embeddings_echa <- readRDS(embeddings_file) 


# embeddings is an (n × 384) numpy matrix
emb_matrix <- as.matrix(embeddings_echa)

# ------------------------------------------------------------------------------
# 8b: Optimal k via silhouette score (k = 2 … 12)
# ------------------------------------------------------------------------------



set.seed(42)
k_vals     <- 2:12
sil_scores <- sapply(k_vals, function(k) {
  km  <- kmeans(emb_matrix, centers = k, nstart = 10, iter.max = 100)
  sil <- silhouette(km$cluster, dist(emb_matrix))
  mean(sil[, 3])
})

sil_df <- data.frame(k = k_vals, silhouette = sil_scores)
k_opt  <- sil_df$k[which.max(sil_df$silhouette)]
message(sprintf("Optimal k = %d  (mean silhouette = %.3f)", k_opt, max(sil_scores)))

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
    subtitle = "Higher is better — dashed line marks the optimal k",
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

# ------------------------------------------------------------------------------
# 8c: k-means clustering with optimal k
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 6 clusters
# ------------------------------------------------------------------------------
k_opt = 6L


set.seed(42)
k  <- k_opt

km <- kmeans(emb_matrix, centers = k, nstart = 25, iter.max = 100)
non_structure$cluster <- factor(km$cluster)

# ------------------------------------------------------------------------------
# 8d: UMAP for 2-D visualisation
# ------------------------------------------------------------------------------


cluster_labels_manual <- c(
  # Cluster 1: All "Reaction mass of..." — defined mixtures of isomers/intermediates
  "Reaction masses",
  # Cluster 2: Complex individual substances — dyes, metal complexes, C8-18 surfactants
  "Complex substances (UVCB with structure)",
  # Cluster 3: Petroleum & coal tar fractions — UVCB hydrocarbons with Cn-Cm ranges
  "Petroleum & coal tar fractions",
  # Cluster 4: Inorganic compounds, metal salts, spinels, slags, nanomaterials
  "Inorganic compounds & metal salts",
  # Cluster 5: Trade names, alphanumeric codes, plant extracts, micro-organisms
  "Trade names, codes & biological materials",
  # Cluster 6: Polymers, fatty acid derivatives, ethoxylated surfactants
  "Polymers, fatty acids & surfactants"
)



umap_coords <- uwot::umap(
  emb_matrix,
  n_neighbors = 15L,
  min_dist    = 0.1,
  metric      = "cosine",
  seed        = 42L
)

non_structure$umap1 <- umap_coords[, 1]
non_structure$umap2 <- umap_coords[, 2]


# Create lookup table
cluster_map <- data.frame(
  cluster = factor(1:length(cluster_labels_manual)),
  manual_label = cluster_labels_manual
)

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
    size = 3,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("#e05c5c", "#4a90d9", "#7bc67e", "#9b59b6", "#e07b3a", "#2ec4b6"),
    labels = cluster_labels_manual
  ) +
  labs(
    title    = "UMAP of substance names without InChIKey (sentence embeddings)",
    subtitle = paste0(
      nrow(non_structure), " substance names \u2192 all-MiniLM-L6-v2 embeddings \u2192 ",
      k, " k-means clusters"
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

# ------------------------------------------------------------------------------
# 8e: t-SNE for comparison / validation
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

non_structure$tsne1 <- tsne_out$Y[, 1]
non_structure$tsne2 <- tsne_out$Y[, 2]

p8_tsne <- ggplot(non_structure, aes(x = tsne1, y = tsne2, colour = cluster)) +
  geom_point(size = 1.2, alpha = 0.6) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title    = "t-SNE of non_structure substance names (sentence embeddings)",
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
ggsave(p8_tsne, 
       filename = here("output", "figures", "Analysis_8e_tsne.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")
  



all_substances2 <- merge(all_substances,plot_data_umap,by="substance_name", all.x = TRUE) 

# ------------------------------------------------------------------------------
# 8f: Top terms per cluster (most central substance names)
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# 8g: Save cluster assignments for use in 06_create_scheme.R
# ------------------------------------------------------------------------------

saveRDS(plot_data_umap, here("data", "processed", "non_structure_clusters.rds"))
message("Analysis 8g: cluster assignments saved to data/processed/non_structure_clusters.rds")



# ==============================================================================
# Analysis 9: Embedding of chemont
# ==============================================================================

chemont_rds <- here("data", "source", "ChemOnt_2_1.rds")

chemont <- readRDS(chemont_rds)

# chemont$substance_name <- paste(
#   chemont$label,
#   chemont$definition,
#   chemont$altLabel
# )
chemont$substance_name <- chemont$label
chemont <- chemont |> select('substance', 'substance_name')


# ------------------------------------------------------------------------------
# 9a: Embed substance names with sentence-transformers
# ------------------------------------------------------------------------------

chemont_embeddings_file <- here("data", "processed", "embeddings_chemont_only_label.rds")

if (!file.exists(chemont_embeddings_file)) {
  
  message(sprintf("Analysis 9: embedding %d chemont substance names", nrow(chemont)))
  
  st <- reticulate::import("sentence_transformers")
  model_chemont <- st$SentenceTransformer("all-MiniLM-L6-v2")
  
  embeddings_chemont <- model_chemont$encode(
    chemont$substance_name,
    show_progress_bar = TRUE,
    convert_to_numpy   = TRUE
  )
  
  saveRDS(embeddings_chemont, chemont_embeddings_file)
  
  
} 

# load
embeddings_chemont <- readRDS(chemont_embeddings_file) 


# embeddings is an (n × 384) numpy matrix
emb_chemont_matrix <- as.matrix(embeddings_chemont)


# ------------------------------------------------------------------------------
# 9b: Cosine similarity
# ------------------------------------------------------------------------------

# normalise
normalize <- function(x) {
  x / sqrt(rowSums(x^2))
}

chemont_norm <- normalize(embeddings_chemont)
echa_norm    <- normalize(embeddings_echa)

# cosine similarity = matrix multiplication
sim_matrix <- echa_norm %*% t(chemont_norm)

# ------------------------------------------------------------------------------
# 9c: Beste match per substance
# ------------------------------------------------------------------------------

best_idx   <- max.col(sim_matrix)
best_score <- sim_matrix[cbind(1:nrow(sim_matrix), best_idx)]

matches <- tibble(
  substance_name     = non_structure$substance_name,
  chemont_label = chemont$substance_name[best_idx],
  score         = best_score
)

# ------------------------------------------------------------------------------
# 9d: filter
# ------------------------------------------------------------------------------

matches <- matches %>%
  mutate(
    matchable = !str_detect(
      substance_name,
      regex("reaction mass|petroleum|distillate|UVCB", ignore_case = TRUE)
    )
  )

# ------------------------------------------------------------------------------
# 9e:  Score distribution
# ------------------------------------------------------------------------------

p9e <- ggplot(matches, aes(score)) +
  geom_histogram(bins = 60) +
  theme_minimal() +
  labs(title = "Similarity score distribution")
print(p9e)

ggsave(p9e, 
       filename = here("output", "figures", "Analysis_9e_distribution_only_label.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")


# ------------------------------------------------------------------------------
# 9f: Threshold tuning (coverage)
# ------------------------------------------------------------------------------
thresholds <- seq(0.3, 0.9, by = 0.02)


coverage <- map_df(thresholds, function(t) { 
   tibble( 
     threshold = t, 
     n_matches = sum(matches$score >= t), 
     coverage = mean(matches$score >= t) 
     ) 
  }
)


p9f <- ggplot(coverage, aes(threshold, coverage)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Coverage vs threshold",
    x = "Similarity threshold",
    y = "Fraction matched"
  )
print(p9f)

ggsave(p9f, 
       filename = here("output", "figures", "Analysis_9f_coverage_threshold.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")

# ------------------------------------------------------------------------------
# 9g: Precision estimation
# ------------------------------------------------------------------------------

p9g <- ggplot(matches, aes(x = score)) +
  stat_ecdf() +
  theme_minimal() +
  labs(title = "ECDF of similarity scores")
print(p9g)
ggsave(p9f, 
       filename = here("output", "figures", "Analysis_9g_precision_estimation.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")

# ------------------------------------------------------------------------------
# 9h: Final output with threshold
# ------------------------------------------------------------------------------

threshold <- 0.7

final_matches <- matches %>%
  filter(score >= threshold)#, matchable)

head(final_matches, 20)

# ------------------------------------------------------------------------------
# 9i: Top-3 matches
# ------------------------------------------------------------------------------

top_k <- 3

top_idx <- t(apply(sim_matrix, 1, function(x) {
  order(x, decreasing = TRUE)[1:top_k]
}))

top_scores <- t(apply(sim_matrix, 1, function(x) {
  sort(x, decreasing = TRUE)[1:top_k]
}))

top_matches <- map_dfr(1:nrow(top_idx), function(i) {
  tibble(
    substance_name = non_structure$substance_name[i],
    chemont_label = chemont$substance_name[top_idx[i, ]],
    score = top_scores[i, ]
  )
})

threshold <- 0.8

final_top_matches <- top_matches %>%
  filter(score >= threshold)

best_matches <- top_matches %>%
  group_by(substance_name) %>%
  slice_max(score, n = 1) %>%
  ungroup()


# ------------------------------------------------------------------------------
# 9j: Ranking per substance
# ------------------------------------------------------------------------------


top_ranked <- top_matches %>%
  group_by(substance_name) %>%
  arrange(desc(score)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

top_wide <- top_ranked %>%
  select(substance_name, rank, score, chemont_label) %>%
  pivot_wider(
    names_from = rank,
    values_from = c(score, chemont_label),
    names_sep = "_"
  )
# ------------------------------------------------------------------------------
# 9k: Consistenty metrics
# ------------------------------------------------------------------------------
top_wide <- top_wide %>%
  mutate(
    gap_12 = score_1 - score_2,
    gap_23 = score_2 - score_3,
    mean_top3 = (score_1 + score_2 + score_3) / 3
  )
# ------------------------------------------------------------------------------
# 9l: Confidence score
# ------------------------------------------------------------------------------

top_wide <- top_wide %>%
  mutate(
    confidence =
      score_1 * 0.6 +      # absolute strength
      gap_12 * 0.3 +       # separation
      gap_23 * 0.1         # stability
  )

# ------------------------------------------------------------------------------
# 9m: Filters
# ------------------------------------------------------------------------------
filtered <- top_wide %>%
  filter(
    score_1 > 0.8,     # minimale kwaliteit
    gap_12 > 0.05      # duidelijk verschil
  )

# ------------------------------------------------------------------------------
# 9n: Visualisation
# ------------------------------------------------------------------------------

p9n <- ggplot(top_wide, aes(x = gap_12, y = score_1)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Match quality landscape",
    x = "Gap between top 1 and 2",
    y = "Top score"
  )
print(p9n)
ggsave(p9n, 
       filename = here("output", "figures", "Analysis_9n_quality_landscape.pdf"),
       device = "pdf",
       height = 25, width = 50, units = "cm")

# ------------------------------------------------------------------------------
# 9o: final matches
# ------------------------------------------------------------------------------
final_matches <- filtered %>%
  select(
    substance_name,
    best_match = chemont_label_1,
    score = score_1,
    confidence
  )

write.csv(
  final_matches,
  file = here("output", "tables", "final_matches.csv"),
  row.names = FALSE
)


message("03_analysis.R: analyses 1–9 completed")


# ==============================================================================
# Analysis 10: Workload — pairwise group-relation matching in sommatie_stoffen
# How many person-days are needed to determine all pairwise set relations?
# ==============================================================================

n_groups          <- all_substances |>
  filter(source == "sommatie_stoffen") |>
  distinct(substance_name) |>
  nrow()

relations_per_day <- 100L  # assumed throughput: group relations assessed per day

days_required <- function(n, y) ceiling(n * (n - 1L) / 2L / y)
marginal_days <- function(n, y) if (n <= 1L) 0L else ceiling((n - 1L) / y)

n_pairs     <- n_groups * (n_groups - 1L) / 2L
total_days  <- days_required(n_groups, relations_per_day)

message(sprintf(
  "Analysis 10: %d groups → %d unique pairs → %d person-days (@ %d relations/day)",
  n_groups, n_pairs, total_days, relations_per_day
))

# 10a: Cumulative days required as a function of n
p10a_data <- data.frame(
  n    = seq_len(n_groups),
  days = sapply(seq_len(n_groups), days_required, y = relations_per_day)
)

p10a <- ggplot(p10a_data, aes(x = n, y = days)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_vline(xintercept = n_groups, linetype = "dashed", colour = "#e05c5c") +
  annotate("label",
           x = n_groups, y = max(p10a_data$days) * 0.5,
           label = paste0("n = ", n_groups, "\n", total_days, " days"),
           colour = "#e05c5c", size = 3.5, label.size = 0.3) +
  labs(
    title    = paste("Analysis 10a: Days required at", relations_per_day, "relations per day"),
    subtitle = "Quadratic growth: n\u00d7(n\u22121)/2 pairs",
    x        = "Number of groups (n)",
    y        = "Days required"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p10a)
ggsave(p10a,
       filename = here("output", "figures", "Analysis_10a_workload_curve.pdf"),
       device = "pdf", height = 5, width = 10, units = "in")

# 10b: Marginal days per additional group
p10b_data <- data.frame(
  n         = 2L:n_groups,
  marg_days = sapply(2L:n_groups, marginal_days, y = relations_per_day)
)

p10b <- ggplot(p10b_data, aes(x = n, y = marg_days)) +
  geom_line(colour = "firebrick", linewidth = 1) +
  labs(
    title    = paste("Analysis 10b: Marginal days per additional group at", relations_per_day, "relations/day"),
    subtitle = paste0(
      "Adding 1 group to the existing ", n_groups,
      " groups costs ", marginal_days(n_groups, relations_per_day), " person-day(s)"
    ),
    x        = "Number of groups (n)",
    y        = "Marginal days"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p10b)
ggsave(p10b,
       filename = here("output", "figures", "Analysis_10b_workload_marginal.pdf"),
       device = "pdf", height = 5, width = 10, units = "in")


# ==============================================================================
# Analysis 11: Ambition vs. FTE — how many FTE are needed to meet the 2030 target?
# ==============================================================================

end_date             <- as.Date("2030-01-01")
start_date           <- Sys.Date()
remaining_years      <- as.numeric(difftime(end_date, start_date, units = "days")) / 365.25
working_days_to_2030 <- ceiling(remaining_years * 200)  # 200 working days per year

fte_required <- ceiling(total_days / working_days_to_2030)

message(sprintf(
  "Analysis 11: %d working days until 2030 | %d days of work | %d FTE required",
  working_days_to_2030, total_days, fte_required
))


# ==============================================================================
# Analysis 12: Pairwise source overlap — Jaccard heatmap + UpSet plot
# ==============================================================================

inchi_src <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(source, inchikey)

sources_vec <- sort(unique(inchi_src$source))

# Build pairwise Jaccard matrix
jaccard_mat <- outer(sources_vec, sources_vec, FUN = Vectorize(function(a, b) {
  set_a <- inchi_src$inchikey[inchi_src$source == a]
  set_b <- inchi_src$inchikey[inchi_src$source == b]
  n_int <- length(intersect(set_a, set_b))
  n_uni <- length(union(set_a, set_b))
  if (n_uni == 0L) 0 else n_int / n_uni
}))
rownames(jaccard_mat) <- sources_vec
colnames(jaccard_mat) <- sources_vec

# Hierarchical clustering on dissimilarity (1 - Jaccard)
hc        <- hclust(as.dist(1 - jaccard_mat), method = "average")
src_order <- sources_vec[hc$order]

jaccard_long <- as.data.frame(jaccard_mat) |>
  tibble::rownames_to_column("source_a") |>
  tidyr::pivot_longer(-source_a, names_to = "source_b", values_to = "jaccard") |>
  mutate(
    source_a = factor(source_a, levels = src_order),
    source_b = factor(source_b, levels = src_order)
  )

p12a <- ggplot(jaccard_long, aes(x = source_a, y = source_b, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(jaccard > 0.02, sprintf("%.2f", jaccard), "")),
            size = 2.5, colour = "grey20") +
  scale_fill_gradient(low = "white", high = "#084594", limits = c(0, 1),
                      name = "Jaccard") +
  coord_fixed() +
  labs(
    title    = "Analysis 12a: Pairwise Jaccard index between sources (on InChIKey)",
    subtitle = "Hierarchically clustered; 1 = complete overlap, 0 = no overlap",
    x        = NULL,
    y        = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.subtitle   = element_text(colour = "grey40"),
    legend.position = "right"
  )

print(p12a)
ggsave(p12a,
       filename = here("output", "figures", "Analysis_12a_Jaccard_heatmap.pdf"),
       device = "pdf", height = 10, width = 11, units = "in")

# 12b: UpSet plot restricted to obligation lists
obligation_sources <- c(
  "restriction_list", "candidate_list", "authorisation_list",
  "pops_list", "eu_positive_list", "harmonised_list", "svhc_identification"
)

upset_obl <- inchi_src |>
  filter(source %in% obligation_sources) |>
  mutate(present = 1L) |>
  tidyr::pivot_wider(names_from = source, values_from = present, values_fill = 0L) |>
  select(-inchikey) |>
  as.data.frame()

pdf(here("output", "figures", "Analysis_12b_UpSet_obligation_lists.pdf"),
    width = 14, height = 7)
UpSetR::upset(
  upset_obl,
  nsets           = ncol(upset_obl),
  nintersects     = 25,
  order.by        = "freq",
  mb.ratio        = c(0.55, 0.45),
  main.bar.color  = "#4a90d9",
  sets.bar.color  = "#e07b3a",
  text.scale      = c(1.2, 1.1, 1, 1, 1.1, 1),
  mainbar.y.label = "Overlap (shared InChIKeys)",
  sets.x.label    = "InChIKeys per list"
)
grid::grid.text(
  "Analysis 12b: Overlap of obligation lists (InChIKey)",
  x = 0.65, y = 0.97,
  gp = grid::gpar(fontsize = 12, fontface = "bold")
)
dev.off()


# ==============================================================================
# Analysis 13: Prioritisation via ClassyFire cache
# Which ChemOnt classes cover the most regulatory substances?
# ==============================================================================

# NULL-safe fallback operator
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Read direct_parent per InChIKey from the local ClassyFire cache
classyfire_cache_dir <- here("data", "cache", "classyfire")
cf_files <- list.files(classyfire_cache_dir, pattern = "\\.json$", full.names = TRUE)

cf_lookup <- purrr::map_dfr(cf_files, function(f) {
  parsed <- tryCatch(jsonlite::read_json(f), error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed[["inchikey"]]) || is.null(parsed[["direct_parent"]]))
    return(NULL)
  tibble::tibble(
    inchikey          = sub("^InChIKey=", "", parsed[["inchikey"]]),
    direct_parent_id  = parsed[["direct_parent"]][["chemont_id"]] %||% NA_character_,
    direct_parent_lbl = parsed[["direct_parent"]][["name"]]       %||% NA_character_
  )
})

# Join with all_substances to determine which sources contain each substance
substance_class <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(inchikey, source) |>
  inner_join(cf_lookup, by = "inchikey")

# Count unique substances and sources per ChemOnt class
class_coverage <- substance_class |>
  filter(!is.na(direct_parent_id)) |>
  group_by(direct_parent_id, direct_parent_lbl) |>
  summarise(
    n_inchikeys = n_distinct(inchikey),
    n_sources   = n_distinct(source),
    .groups     = "drop"
  ) |>
  arrange(desc(n_inchikeys))

saveRDS(class_coverage, here("output", "tables", "class_coverage.rds"))

# 13a: Top 30 ChemOnt classes by regulatory substance coverage
top30_classes <- class_coverage |> slice_head(n = 30)

p13a <- ggplot(top30_classes,
               aes(x = reorder(direct_parent_lbl, n_inchikeys),
                   y = n_inchikeys, fill = n_sources)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_inchikeys), hjust = -0.1, size = 3) +
  scale_fill_gradient(low = "#b3cde3", high = "#084594", name = "Number of\nsources") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title    = "Analysis 13a: Top 30 ChemOnt classes by regulatory substance coverage",
    subtitle = "Colour = number of sources; bar length = number of unique InChIKeys",
    x        = NULL,
    y        = "Number of unique substances (InChIKey)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p13a)
ggsave(p13a,
       filename = here("output", "figures", "Analysis_13a_ChemOnt_prioritisation_top30.pdf"),
       device = "pdf", height = 8, width = 12, units = "in")

# 13b: Histogram — distribution of substance coverage across ChemOnt classes
p13b <- ggplot(class_coverage, aes(x = n_inchikeys)) +
  geom_histogram(binwidth = 5, fill = "#4a90d9", colour = "white") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Analysis 13b: Distribution of substance coverage across ChemOnt classes",
    subtitle = "How many classes cover 1, 5, 10 \u2026 substances?",
    x        = "Number of unique substances per ChemOnt class",
    y        = "Number of ChemOnt classes"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

print(p13b)
ggsave(p13b,
       filename = here("output", "figures", "Analysis_13b_ChemOnt_coverage_histogram.pdf"),
       device = "pdf", height = 5, width = 10, units = "in")

message(sprintf(
  "Analysis 13: %d ChemOnt classes with coverage; top class '%s' covers %d substances",
  nrow(class_coverage),
  class_coverage$direct_parent_lbl[1],
  class_coverage$n_inchikeys[1]
))

message("03_analysis.R completed")
