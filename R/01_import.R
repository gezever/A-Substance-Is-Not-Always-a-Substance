# ==============================================================================
# 01_import.R
# Load raw source data: ECHA/EU chemical substance lists
# ==============================================================================

library(readxl)
library(here)

src <- here("data", "source")

# ------------------------------------------------------------------------------
# Obligation lists (chem.echa.europa.eu/obligation-lists)
# ------------------------------------------------------------------------------

restriction_list    <- read_excel(file.path(src, "restriction_list_full.xlsx"),    sheet = "List_results")
candidate_list      <- read_excel(file.path(src, "candidate_list_full.xlsx"),      sheet = "List_results")
authorisation_list  <- read_excel(file.path(src, "authorisation_list_full.xlsx"),  sheet = "List_results")
pops_list           <- read_excel(file.path(src, "pops_list_full.xlsx"),           sheet = "List_results")
eu_positive_list    <- read_excel(file.path(src, "eu_positive_list_full.xlsx"),    sheet = "List_results")
harmonised_list     <- read_excel(file.path(src, "Harmonised_List.xlsx"),          sheet = "List_results")

# ------------------------------------------------------------------------------
# Activity lists (chem.echa.europa.eu/activity-lists)
# ------------------------------------------------------------------------------

restriction_process    <- read_excel(file.path(src, "restriction_process_full.xlsx"),    sheet = "List_results")
svhc_identification    <- read_excel(file.path(src, "svhc_identification_full.xlsx"),    sheet = "List_results")
authorisation_process  <- read_excel(file.path(src, "authorisation_process_full.xlsx"),  sheet = "List_results")
dossier_evaluation     <- read_excel(file.path(src, "dossier_evaluation_full.xlsx"),     sheet = "List_results")
clh_process            <- read_excel(file.path(src, "clh_process_full.xlsx"),            sheet = "List_results")
substance_evaluation   <- read_excel(file.path(src, "substance_evaluation_full.xlsx"),   sheet = "List_results")
pops_process           <- read_excel(file.path(src, "pops_process_full.xlsx"),           sheet = "List_results")
pbt_assessment         <- read_excel(file.path(src, "pbt_assessment.xlsx"),              sheet = "List_results")
ed_assessment          <- read_excel(file.path(src, "ed_assessment.xlsx"),               sheet = "List_results")

# ------------------------------------------------------------------------------
# REACH registrations (chem.echa.europa.eu)
# ------------------------------------------------------------------------------

reach_registrations <- read_excel(file.path(src, "reach_registrations.xlsx"), sheet = "Substances list")

# ------------------------------------------------------------------------------
# EU Pesticides Database — Active Substances
# Row 1-2: title/metadata; row 3: column headers; data from row 4 onward
# ------------------------------------------------------------------------------

pesticides <- read_excel(
  file.path(src, "Pesticides_ActiveSubstanceExport.xlsx"),
  sheet = "Active Substance Search export",
  skip  = 2
)

message("01_import.R completed")
