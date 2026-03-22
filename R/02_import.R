
# ==============================================================================
# 01_import.R
# Load raw source data: ECHA/EU chemical substance lists
# ==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(here)
library(httr)
library(jsonlite)
library(stringr)

src <- here("data", "source")


# ------------------------------------------------------------------------------
# Helper function: multiple numbers in one cell are split
# ------------------------------------------------------------------------------

expand_df_on_whitespace <- function(df) {
  for (col in c("ec_number", "cas_number")) {
    df <- df %>%
      separate_rows(all_of(col), sep = "\\s+") %>%
      distinct()
  }
  return(df)
}

# ------------------------------------------------------------------------------
# Helper function: read list and standardise column names
# ------------------------------------------------------------------------------

read_list <- function(file, sheet = "List_results", source_name, skip = 0) {
  
  read_excel(file.path(src, file), sheet = sheet, skip = skip) |>
    rename(
      substance_name = any_of(c("Substance name", "Name", "Substance")),
      cas_number     = any_of(c("CAS number", "CAS Number")),
      ec_number      = any_of(c("EC number", "EC Number"))
    ) |>
    mutate(source = source_name)
}
# ------------------------------------------------------------------------------
# VMM lists
# ------------------------------------------------------------------------------

sommatie_stoffen <- readRDS(here("data", "source", "sommatie_stoffen.rds"))  |>
  rename(
    substance_name = any_of(c("altLabel_en")),
    cas_number     = any_of(c("casNumber")),
    ec_number      = any_of(c("ecNumber"))
  ) |>
  mutate(source = 'sommatie_stoffen')|>
  subset(inScheme == 'https://data.omgeving.vlaanderen.be/id/conceptscheme/sommatie_stoffen')
# ------------------------------------------------------------------------------
# Obligation lists
# ------------------------------------------------------------------------------

restriction_list    <- read_list("restriction_list_full.xlsx",    source_name = "restriction_list")
candidate_list      <- read_list("candidate_list_full.xlsx",      source_name = "candidate_list")
authorisation_list  <- read_list("authorisation_list_full.xlsx",  source_name = "authorisation_list")
pops_list           <- read_list("pops_list_full.xlsx",           source_name = "pops_list")
eu_positive_list    <- read_list("eu_positive_list_full.xlsx",    source_name = "eu_positive_list")
harmonised_list     <- read_list("Harmonised_List.xlsx",          source_name = "harmonised_list")

# ------------------------------------------------------------------------------
# Activity lists
# ------------------------------------------------------------------------------

restriction_process    <- read_list("restriction_process_full.xlsx",    source_name = "restriction_process")
svhc_identification    <- read_list("svhc_identification_full.xlsx",    source_name = "svhc_identification")
authorisation_process  <- read_list("authorisation_process_full.xlsx",  source_name = "authorisation_process")
dossier_evaluation     <- read_list("dossier_evaluation_full.xlsx",     source_name = "dossier_evaluation")
clh_process            <- read_list("clh_process_full.xlsx",            source_name = "clh_process")
substance_evaluation   <- read_list("substance_evaluation_full.xlsx",   source_name = "substance_evaluation")
pops_process           <- read_list("pops_process_full.xlsx",           source_name = "pops_process")
pbt_assessment         <- read_list("pbt_assessment.xlsx",              source_name = "pbt_assessment")
ed_assessment          <- read_list("ed_assessment.xlsx",               source_name = "ed_assessment")

# ------------------------------------------------------------------------------
# REACH registrations
# ------------------------------------------------------------------------------

reach_registrations <- read_list(
  "reach_registrations.xlsx",
  sheet = "Substances list",
  source_name = "reach_registrations"
)

# ------------------------------------------------------------------------------
# EU Pesticides Database
# ------------------------------------------------------------------------------

pesticides <- read_list(
  "Pesticides_ActiveSubstanceExport.xlsx",
  sheet = "Active Substance Search export",
  skip  = 2,
  source_name = "pesticides"
)

# ------------------------------------------------------------------------------
# Combine all sources
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
  pesticides,
  sommatie_stoffen
)|>
  select(substance_name, ec_number, cas_number, source) |> 
  distinct()|>
  mutate(
    cas_number = str_trim(cas_number),
    ec_number = str_trim(ec_number)
  )

# ------------------------------------------------------------------------------
# Cleaning
# ------------------------------------------------------------------------------

pattern <- "^[0-9]+-[0-9]+-[0-9]+.*$"

all_substances$ec_number[!grepl(pattern, all_substances$ec_number)] <- NA
all_substances$cas_number[!grepl(pattern, all_substances$cas_number)] <- NA

all_substances <- expand_df_on_whitespace(all_substances)

all_substances$ec_number[!grepl(pattern, all_substances$ec_number)] <- NA
all_substances$cas_number[!grepl(pattern, all_substances$cas_number)] <- NA


# ------------------------------------------------------------------------------
# Unique substances (CAS + EC)
# ------------------------------------------------------------------------------

unique_substances <- all_substances |>
  select(cas_number) |> distinct()

#unique_substances <- all_substances |>
#  distinct(ec_number, cas_number, .keep_all = TRUE) |>
#  arrange(substance_name)

# ------------------------------------------------------------------------------
# PubChem InChIKey lookup with caching
# ------------------------------------------------------------------------------

inchikey_cache_dir <- here("data", "cache", "inchikey")
dir.create(inchikey_cache_dir, recursive = TRUE, showWarnings = FALSE)

get_inchikey <- function(cas) {

  if (is.na(cas) || cas == "-" || !nzchar(trimws(cas)))
    return(NA_character_)

  cache_file <- file.path(inchikey_cache_dir, paste0(cas, ".json"))

  if (file.exists(cache_file)) {
    json <- fromJSON(paste(readLines(cache_file, warn = FALSE), collapse = "\n"))
    return(json$PropertyTable$Properties$InChIKey)
  } #else {    return(NA_character_)  }

  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
    URLencode(cas, reserved = TRUE),
    "/property/InChIKey/JSON"
  )

  message(paste('handling ',cas))

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

unique_substances <- unique_substances |>
  mutate(inchikey = lapply(cas_number, get_inchikey)) |>
  unnest(inchikey)



# ------------------------------------------------------------------------------
# Add InChIKey back to all substances
# ------------------------------------------------------------------------------

lookup <- unique_substances |>
  distinct(cas_number, inchikey)

all_substances$inchikey <- lookup$inchikey[
  match(all_substances$cas_number, lookup$cas_number)
]

all_substances <- all_substances |>
  select(substance_name, ec_number, cas_number, inchikey, source) |> distinct()

# ------------------------------------------------------------------------------
# Write all substances 
# ------------------------------------------------------------------------------

saveRDS(all_substances, here("data", "processed", "all_substances.rds"))
write.csv(all_substances, here("data", "processed", "all_substances.csv"))
#colnames(all_substances)

message("01_import.R completed")
