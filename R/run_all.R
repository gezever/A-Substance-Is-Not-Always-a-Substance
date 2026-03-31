# ==============================================================================
# run_all.R
# Executes the full analysis pipeline in dependency order
# ==============================================================================

library(here)

scripts <- c(
  "01_download.R",
  "02_import.R",
  "03_analysis.R",                           # Analyses 1-3: substance identity and identifier consistency
  "04_entity_classification.R",              # Analysis 4:  entity type classification and linkability taxonomy
  "05_overlap_lists.R",                      # Analysis 5:  overlap between regulatory lists (UpSet)
  "06_coverage_linking.R",                   # Analysis 6:  coverage of structure-based linking per source
  "07_network_visualisation.R",              # Analysis 7:  bipartite substance ↔ list network
  "08_embedding_clustering.R",               # Analysis 8:  sentence embedding and clustering of non-structure names
  "09_embedding_chemont.R",                  # Analysis 9:  cosine similarity matching to ChemOnt classes
  "10_workload.R",                           # Analysis 10: pairwise group-relation workload estimation
  "11_ambition_fte.R",                       # Analysis 11: FTE required to meet 2030 target
  "12_pairwise_overlap.R",                   # Analysis 12: pairwise Jaccard heatmap + obligation UpSet
  "13a_classyfire_coverage.R",               # Analysis 13: ChemOnt class coverage via ClassyFire cache
  "13b_before_prioritisation_create_scheme.R", # Schema creation (prerequisite for 14)
  "14_prioritization.R"                      # Analysis 14: composite priority scoring and visualisations
)

for (script in scripts) {
  message("\n=== ", script, " ===")
  source(here("R", script))
}

message("\nDone. Figures: output/figures/  |  Tables: output/tables/")
