# ==============================================================================
# 01_import.R
# Load raw source data: ECHA/EU chemical substance lists
# ==============================================================================

library(readxl)
library(dplyr)
library(here)
library(httr)
library(jsonlite)

src <- here("data", "source")

# ------------------------------------------------------------------------------
# Obligation lists (chem.echa.europa.eu/obligation-lists)
# ------------------------------------------------------------------------------

restriction_list    <- read_excel(file.path(src, "restriction_list_full.xlsx"),    sheet = "List_results") |> mutate(source = "restriction_list")
candidate_list      <- read_excel(file.path(src, "candidate_list_full.xlsx"),      sheet = "List_results") |> mutate(source = "candidate_list")
authorisation_list  <- read_excel(file.path(src, "authorisation_list_full.xlsx"),  sheet = "List_results") |> mutate(source = "authorisation_list")
pops_list           <- read_excel(file.path(src, "pops_list_full.xlsx"),           sheet = "List_results") |> mutate(source = "pops_list")
eu_positive_list    <- read_excel(file.path(src, "eu_positive_list_full.xlsx"),    sheet = "List_results") |> mutate(source = "eu_positive_list")
harmonised_list     <- read_excel(file.path(src, "Harmonised_List.xlsx"),          sheet = "List_results") |> mutate(source = "harmonised_list")

# ------------------------------------------------------------------------------
# Activity lists (chem.echa.europa.eu/activity-lists)
# ------------------------------------------------------------------------------

restriction_process    <- read_excel(file.path(src, "restriction_process_full.xlsx"),    sheet = "List_results") |> mutate(source = "restriction_process")
svhc_identification    <- read_excel(file.path(src, "svhc_identification_full.xlsx"),    sheet = "List_results") |> mutate(source = "svhc_identification")
authorisation_process  <- read_excel(file.path(src, "authorisation_process_full.xlsx"),  sheet = "List_results") |> mutate(source = "authorisation_process")
dossier_evaluation     <- read_excel(file.path(src, "dossier_evaluation_full.xlsx"),     sheet = "List_results") |> mutate(source = "dossier_evaluation")
clh_process            <- read_excel(file.path(src, "clh_process_full.xlsx"),            sheet = "List_results") |> mutate(source = "clh_process")
substance_evaluation   <- read_excel(file.path(src, "substance_evaluation_full.xlsx"),   sheet = "List_results") |> mutate(source = "substance_evaluation")
pops_process           <- read_excel(file.path(src, "pops_process_full.xlsx"),           sheet = "List_results") |> mutate(source = "pops_process")
pbt_assessment         <- read_excel(file.path(src, "pbt_assessment.xlsx"),              sheet = "List_results") |> mutate(source = "pbt_assessment")
ed_assessment          <- read_excel(file.path(src, "ed_assessment.xlsx"),               sheet = "List_results") |> mutate(source = "ed_assessment")

# ------------------------------------------------------------------------------
# REACH registrations (chem.echa.europa.eu)
# ------------------------------------------------------------------------------

reach_registrations <- read_excel(file.path(src, "reach_registrations.xlsx"), sheet = "Substances list") |> mutate(source = "reach_registrations")

# ------------------------------------------------------------------------------
# EU Pesticides Database — Active Substances
# Row 1-2: title/metadata; row 3: column headers; data from row 4 onward
# ------------------------------------------------------------------------------

pesticides <- read_excel(
  file.path(src, "Pesticides_ActiveSubstanceExport.xlsx"),
  sheet = "Active Substance Search export",
  skip  = 2
) |> mutate(source = "pesticides")

# ------------------------------------------------------------------------------
# Gecombineerde dataframe (alle bronnen samengevoegd)
# Kolommen die niet in een bron voorkomen worden gevuld met NA
# ------------------------------------------------------------------------------

all_substances <- bind_rows(
  restriction_list,
  candidate_list,
  authorisation_list,
  pops_list,
  eu_positive_list,
  harmonised_list,
  restriction_process,
  svhc_identification,
  authorisation_process,
  dossier_evaluation,
  clh_process,
  substance_evaluation,
  pops_process,
  pbt_assessment,
  ed_assessment,
  reach_registrations,
  pesticides
)

# ------------------------------------------------------------------------------
# Unieke stoffen op basis van CAS- en EC-nummer
# Bronnen gebruiken wisselende kolomnamen; coalesce tot één waarde per rij
# ------------------------------------------------------------------------------

unique_substances <- all_substances |>
  transmute(
    substance_name = coalesce(`Substance name`, Name, Substance),
    ec_number      = coalesce(`EC number`, `EC Number`),
    cas_number     = coalesce(`CAS number`, `CAS Number`)
  ) |>
  distinct(ec_number, cas_number, .keep_all = TRUE) |>
  arrange(substance_name)

# ------------------------------------------------------------------------------
# InChIKey opzoeken via PubChem op basis van CAS-nummer
# JSON-responses worden gecached in data/cache/inchikey/ om herhaalde requests
# te voorkomen. Ongeldige CAS-nummers ("-" of NA) worden overgeslagen.
# PubChem rate limit: max 5 req/s → Sys.sleep(0.2) tussen requests.
# ------------------------------------------------------------------------------

inchikey_cache_dir <- here("data", "cache", "inchikey")
dir.create(inchikey_cache_dir, recursive = TRUE, showWarnings = FALSE)

get_inchikey <- function(cas) {
  if (is.na(cas) || cas == "-" || !nzchar(trimws(cas))) return(NA_character_)

  cache_file <- file.path(inchikey_cache_dir, paste0(cas, ".json"))

  if (file.exists(cache_file)) {
    json <- fromJSON(cache_file)
    return(json$PropertyTable$Properties$InChIKey[1])
  }

  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
    URLencode(cas, reserved = TRUE),
    "/property/InChIKey/JSON"
  )

  res <- GET(url)
  Sys.sleep(0.2)

  if (status_code(res) == 200) {
    raw <- content(res, "text", encoding = "UTF-8")
    writeLines(raw, cache_file)
    json <- fromJSON(raw)
    return(json$PropertyTable$Properties$InChIKey[1])
  } else {
    return(NA_character_)
  }
}

unique_substances$inchikey <- sapply(unique_substances$cas_number, get_inchikey)

message("01_import.R completed")
