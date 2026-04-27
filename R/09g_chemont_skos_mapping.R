# ==============================================================================
# 09g_chemont_skos_mapping.R
# SKOS mapping triples for ChemOnt equivalences from experiment 09f
#
# PURPOSE
# -------
# Converts the 09f equivalence results (is_group, equivalent_node,
# relationship_type) into formal SKOS mapping triples against the ChemOnt
# SKOS file, using the predicate direction:
#
#   subject = substance (legislative entry)
#   object  = ChemOnt concept
#
#   relationship_type = "exact"   → skos:exactMatch
#   relationship_type = "subset"  → skos:broadMatch  (ChemOnt is broader)
#   relationship_type = "broader" → skos:narrowMatch  (ChemOnt is narrower)
#
# URI conventions (consistent with existing RDF pipeline):
#   Substance : MD5 hash of substance_name
#     https://data.omgeving.vlaanderen.be/id/concept/chemical_substance/{MD5}
#   ChemOnt   : numeric CHEMONTID from ChemOnt_2_1_skos.ttl
#     https://data.omgeving.vlaanderen.be/id/concept/chemical_substance/{CHEMONTID}
#
# DATA PROVENANCE
# ---------------
# Input 1: output/tables/Analysis_9f_group_equivalence.csv
# Input 2: data/processed/chemont_label_uri.csv  (built by chemont_label_uri.py)
#
# OUTPUTS
# -------
# data/processed/rdf/chemont_nonstructure_mapping.ttl
# output/tables/Analysis_9g_skos_mapping.csv
# ==============================================================================

library(digest)
library(dplyr)
library(here)
library(readr)
library(stringr)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

RESULTS_CSV   <- here("output", "tables", "Analysis_9f_group_equivalence.csv")
LABEL_URI_CSV <- here("data", "processed", "chemont_label_uri.csv")
OUT_TTL       <- here("data", "processed", "rdf", "chemont_nonstructure_mapping.ttl")
OUT_CSV       <- here("output", "tables", "Analysis_9g_skos_mapping.csv")

CONCEPT_BASE  <- "https://data.omgeving.vlaanderen.be/id/concept/chemical_substance/"

SKOS_PREDICATE <- c(
  exact   = "skos:exactMatch",
  subset  = "skos:broadMatch",
  broader = "skos:narrowMatch"
)

dir.create(here("data", "processed", "rdf"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "tables"),          showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Step 1: Build label → URI lookup
# ------------------------------------------------------------------------------

if (!file.exists(LABEL_URI_CSV)) {
  stop("Label → URI lookup not found: ", LABEL_URI_CSV,
       "\nRun: python3 python/chemont_label_uri.py ",
       "data/source/chemont/ChemOnt_2_1_skos.ttl ",
       "data/processed/chemont_label_uri.csv")
}

label_uri <- read_csv(LABEL_URI_CSV, show_col_types = FALSE)
message(sprintf("Label → URI lookup: %d entries", nrow(label_uri)))

# ------------------------------------------------------------------------------
# Step 2: Load 09f results and compute substance URIs
# ------------------------------------------------------------------------------

results <- read_csv(RESULTS_CSV, show_col_types = FALSE)
message(sprintf("09f results: %d rows", nrow(results)))

results <- results |>
  mutate(
    substance_uri = paste0(
      CONCEPT_BASE,
      sapply(substance_name, digest, algo = "md5", USE.NAMES = FALSE)
    )
  )

# ------------------------------------------------------------------------------
# Step 3: Filter to mappable rows and join ChemOnt URIs
# ------------------------------------------------------------------------------

mappings <- results |>
  filter(
    is_group == TRUE,
    !is.na(equivalent_node),
    !is.na(relationship_type),
    relationship_type %in% names(SKOS_PREDICATE)
  ) |>
  left_join(label_uri, by = c("equivalent_node" = "label")) |>
  rename(chemont_uri = uri) |>
  filter(!is.na(chemont_uri)) |>
  mutate(predicate = SKOS_PREDICATE[relationship_type]) |> 
  subset(confidence == "high") # only high confidence

message(sprintf(
  "Mappable rows: %d (exact: %d, subset: %d, broader: %d)",
  nrow(mappings),
  sum(mappings$relationship_type == "exact"),
  sum(mappings$relationship_type == "subset"),
  sum(mappings$relationship_type == "broader")
))

n_unresolved <- results |>
  filter(is_group == TRUE, !is.na(equivalent_node)) |>
  anti_join(mappings, by = "substance_name") |>
  nrow()
if (n_unresolved > 0L)
  message(sprintf("Warning: %d is_group rows with equivalent_node not found in label lookup", n_unresolved))

# ------------------------------------------------------------------------------
# Step 4: Write CSV summary
# ------------------------------------------------------------------------------

mappings |>
  select(
    substance_name, substance_uri, predicate,
    equivalent_node, chemont_uri,
    relationship_type, equivalent_level, confidence,
    entity_type
  ) |> 
  write_csv(OUT_CSV)
message(sprintf("CSV written → %s", OUT_CSV))

# ------------------------------------------------------------------------------
# Step 5: Write Turtle
# ------------------------------------------------------------------------------

escape_literal <- function(x) {
  x |>
    str_replace_all("\\\\", "\\\\\\\\") |>
    str_replace_all('"',    '\\\\"')    |>
    str_replace_all("[\r\n\t]", " ")
}

prefixes <- c(
  "@prefix concept: <https://data.omgeving.vlaanderen.be/id/concept/chemical_substance/> .",
  "@prefix skos:    <http://www.w3.org/2004/02/skos/core#> .",
  "@prefix rdfs:    <http://www.w3.org/2000/01/rdf-schema#> .",
  ""
)

triples <- mappings |>
  mutate(
    block = sprintf(
      paste0(
        "<%s>\n",
        "    %s <%s> ;\n",
        "    rdfs:label \"%s\"@nl ;\n",
        "    skos:note \"confidence: %s\" .\n"
      ),
      substance_uri,
      predicate,
      chemont_uri,
      escape_literal(substance_name),
      confidence
    )
  ) |>
  pull(block)

writeLines(c(prefixes, triples), OUT_TTL)
message(sprintf(
  "TTL written → %s  (%d triples)",
  OUT_TTL, nrow(mappings)
))

# ------------------------------------------------------------------------------
# Step 6: Validate with Apache Jena riot (if available)
# ------------------------------------------------------------------------------

riot <- "/opt/apache-jena-5.6.0/bin/riot"
if (file.exists(riot)) {
  ret <- system2(riot, args = c("--validate", OUT_TTL), stdout = "", stderr = "")
  if (ret == 0L) {
    message("riot validation: OK")
  } else {
    warning("riot validation failed — check TTL syntax")
  }
} else {
  message("riot not found, skipping validation")
}
