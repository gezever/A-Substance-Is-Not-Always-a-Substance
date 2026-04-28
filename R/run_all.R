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
  "09b_chemont_model_comparison.R",          # creates scibert/biobert embedding caches (required by 8i/8j)
  "08_embedding_clustering.R",               # Analysis 8:  sentence embedding and clustering of non-structure names
  "08b_cluster_regex.R",                     # Analysis 8b: regex reverse-engineering of embedding clusters
  "09_embedding_chemont.R",                  # Analysis 9:  cosine similarity matching to ChemOnt classes
  "10_workload.R",                           # Analysis 10: pairwise group-relation workload estimation
  "11_ambition_fte.R",                       # Analysis 11: FTE required to meet 2030 target
  "12_pairwise_overlap.R",                   # Analysis 12: pairwise Jaccard heatmap + obligation UpSet
  "13a_classyfire_coverage.R",               # Analysis 13: ChemOnt class coverage via ClassyFire cache
  "13b_before_prioritisation_create_scheme.R", # Schema creation (prerequisite for 14); creates substances_taxonomy_levels.ttl
  "09c_chemont_model_validation.R",          # Analysis 9c: requires substances_taxonomy_levels.ttl from 13b
  "09d_chemont_llm_experiment.R",            # Analysis 9d: requires Analysis_9c_model_validation.csv from 09c
  "09e_llm_nonstructure.R",                  # Analysis 9e: requires final_matches.csv from 09_embedding_chemont
  "09f_llm_group_equivalence.R",             # Analysis 9f: requires final_matches.csv from 09_embedding_chemont
  "09g_chemont_skos_mapping.R",              # Analysis 9g: requires Analysis_9f_group_equivalence.csv from 09f
  "14_prioritization.R",                      # Analysis 14: composite priority scoring and visualisations

)

for (script in scripts) {
  message("\n=== ", script, " ===")
  source(here("R", script))
}

message("\nDone. Figures: output/figures/  |  Tables: output/tables/")
