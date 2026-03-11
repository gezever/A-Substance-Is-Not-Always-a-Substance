# ==============================================================================
# run_all.R
# Voert de volledige analysepijplijn uit in volgorde
# ==============================================================================

library(here)

scripts <- c(
  "00_download.R",
  "01_import.R",
  "02_clean.R",
  "03_analysis.R",
  "04_visualize.R",
  "05_tables.R"
)

for (script in scripts) {
  message("\n=== ", script, " ===")
  source(here("R", script))
}

message("\nKlaar. Figuren: output/figures/  |  Tabellen: output/tables/")
