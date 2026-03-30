library(here)
library(tidyr)
library(dplyr)
library(jsonlite)
library(stringr)
library(digest)
library(httr)



# functie om dataframe om te zetten naar jsonld
to_jsonld <- function(dataframe) {
  # lees context
  context <- jsonlite::read_json(here("data", "source", "context", "context.json"))
  # jsonld constructie
  df_in_list <- list('@graph' = dataframe, '@context' = context)
  df_in_json <- toJSON(df_in_list, auto_unbox=TRUE)
  return(df_in_json)
}

get_classyfire_direct_parent <- function(inchikey) {
  cache_file <- here("data", "cache", "classyfire", paste0(inchikey, ".json"))
  if (!file.exists(cache_file)) {
    url <- paste0("https://cfb.fiehnlab.ucdavis.edu/entities/", inchikey)
    response <- tryCatch(GET(url), error = function(e) NULL)
    if (is.null(response) || http_error(response)) return(NA_character_)
    writeBin(content(response, "raw"), cache_file)
  }
  parsed <- tryCatch(
    jsonlite::read_json(cache_file),
    error = function(e) NULL
  )
  if (is.null(parsed)) return(NA_character_)
  chemont_id <- parsed$direct_parent$chemont_id
  if (is.null(chemont_id)) NA_character_ else chemont_id
}


all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

all_substances <- all_substances %>%
  mutate(hash = ifelse(
    !is.na(inchikey) & inchikey != "",
    sapply(inchikey, function(x) digest::digest(x, algo = "md5")),
    sapply(substance_name, function(x) digest::digest(x, algo = "md5"))
  )) %>%
  mutate(direct_parent_chemont_id = sapply(inchikey, function(x) {
    if (is.na(x) || x == "") NA_character_
    else get_classyfire_direct_parent(x)
  })) %>% 
  mutate(type = ifelse(!is.na(inchikey) & inchikey != "", "dbo:ChemicalSubstance", NA_character_))

# merge all classyfire cache JSON files into a single JSON-LD document
classyfire_files <- list.files(here("data", "cache", "classyfire"), pattern = "\\.json$", full.names = TRUE)
classyfire_graph <- lapply(classyfire_files, function(f) {
  parsed <- tryCatch(jsonlite::read_json(f), error = function(e) NULL)
  if (!is.null(parsed) && !is.null(parsed[["inchikey"]])) {
    inchikey <- sub("^InChIKey=", "", parsed[["inchikey"]])
    parsed[["id"]] <- digest::digest(inchikey, algo = "md5")
  }
  parsed
})
classyfire_graph <- Filter(Negate(is.null), classyfire_graph)
context <- jsonlite::read_json(here("data", "source", "context", "context.json"))
classyfire_jsonld <- toJSON(list(`@graph` = classyfire_graph, `@context` = context), auto_unbox = TRUE)
write(classyfire_jsonld, here("data", "processed", "rdf", "classyfire.jsonld"))


# ------------------------------------------------------------------------------
# Linkability taxonomy → JSON-LD
# ------------------------------------------------------------------------------

linkability_rds <- here("data", "processed", "linkability_taxonomy.rds")

if (file.exists(linkability_rds)) {

  linkability <- readRDS(linkability_rds)

  cosc_prefix    <- "https://data.omgeving.vlaanderen.be/id/collection/chemical_substance/"
  concept_prefix <- "https://data.omgeving.vlaanderen.be/id/concept/chemical_substance/"

  # Compute hash URI for each substance (same logic as all_substances$hash)
  linkability <- linkability |>
    mutate(
      linkability = as.character(linkability),
      hash = ifelse(
        !is.na(inchikey) & inchikey != "",
        sapply(inchikey,       digest::digest, algo = "md5"),
        sapply(substance_name, digest::digest, algo = "md5")
      ),
      base_hash = ifelse(
        !is.na(base_inchikey) & base_inchikey != "",
        sapply(base_inchikey, digest::digest, algo = "md5"),
        NA_character_
      )
    )

  # Fixed slugs per linkability level
  slug_map <- c(
    "Structure (InChIKey)"                  = "linkability_structure_inchikey",
    "Group \u2014 base compound resolvable"      = "linkability_group_resolvable",
    "Group \u2014 base compound not resolvable"  = "linkability_group_not_resolvable",
    "CAS without structure"                 = "linkability_cas_without_structure",
    "UVCB / Mixture"                        = "linkability_uvcb_mixture",
    "Regulatory entry"                      = "linkability_regulatory_entry",
    "Unidentified"                          = "linkability_unidentified"
  )

  # One skos:Collection block per linkability category
  collection_blocks <- linkability |>
    group_by(linkability) |>
    summarise(hashes = list(unique(hash)), .groups = "drop") |>
    purrr::pmap_chr(function(linkability, hashes) {
      slug    <- slug_map[[linkability]]
      uri     <- paste0("<", cosc_prefix, slug, ">")
      members <- paste0(
        "    skos:member <", concept_prefix, hashes, "> ;",
        collapse = "\n"
      )
      paste0(
        uri, "\n",
        "    a skos:Collection ;\n",
        "    skos:prefLabel \"", linkability, "\"@en ;\n",
        members, "\n",
        "    .\n"
      )
    })

  # skos:narrower triples: group substance → base compound
  narrower_blocks <- linkability |>
    filter(!is.na(base_hash)) |>
    distinct(hash, base_hash) |>
    purrr::pmap_chr(function(hash, base_hash) {
      paste0(
        "<", concept_prefix, hash, ">\n",
        "    skos:narrower <", concept_prefix, base_hash, "> .\n"
      )
    })

  linkability_ttl_path <- here("data", "processed", "rdf", "linkability_taxonomy.ttl")
  header <- paste(
    "PREFIX skos:    <http://www.w3.org/2004/02/skos/core#>",
    paste0("PREFIX concept: <", concept_prefix, ">"),
    paste0("PREFIX cosc:    <", cosc_prefix, ">"),
    "",
    sep = "\n"
  )

  write(
    paste(c(header, collection_blocks, narrower_blocks), collapse = "\n\n"),
    linkability_ttl_path
  )

  message(sprintf(
    "06_create_scheme.R: %d linkability collections + %d narrower triples written to %s",
    length(collection_blocks), length(narrower_blocks), linkability_ttl_path
  ))

} else {
  message("06_create_scheme.R: linkability_taxonomy.rds not found — skipping. Run 03_analysis.R first.")
}

write(to_jsonld(all_substances), here("data", "processed", "rdf", "substances.jsonld"))
system("bash bash/jsonld-to-ttl.bash")


saveRDS(all_substances, here("data", "processed", "all_substances_inchikey.rds"))
all_substances <- readRDS(here("data", "processed", "all_substances_inchikey.rds"))



# ==============================================================================
# Add k-means clusters from Analysis 8 as skos:Collection to the TTL
# ==============================================================================

clusters_rds <- here("data", "processed", "non_structure_clusters.rds")

if (file.exists(clusters_rds)) {

  non_structure_clusters <- readRDS(clusters_rds)

  # Build URI lookup: substance_name → hash (same logic as all_substances$hash)
  hash_lookup <- all_substances |>
    select(substance_name, hash) |>
    distinct()

  # Prefixes already defined in substances_taxonomy.ttl
  cosc_prefix    <- "https://data.omgeving.vlaanderen.be/id/collection/chemical_substance/"
  concept_prefix <- "https://data.omgeving.vlaanderen.be/id/concept/chemical_substance/"

  # Build one skos:Collection block per cluster
  cluster_turtle <- non_structure_clusters |>
    inner_join(hash_lookup, by = "substance_name") |>
    filter(!is.na(hash)) |>
    group_by(cluster, manual_label) |>
    summarise(
      hashes = list(hash),
      .groups = "drop"
    ) |>
    purrr::pmap_chr(function(cluster, manual_label, hashes) {
      slug    <- gsub("[^a-z0-9]+", "_", tolower(manual_label))
      uri     <- paste0("<", cosc_prefix, "cluster_", slug, ">")
      members <- paste0(
        "    skos:member <", concept_prefix, hashes, "> ;",
        collapse = "\n"
      )
      paste0(
        uri, "\n",
        "    a skos:Collection ;\n",
        "    skos:prefLabel \"", manual_label, "\"@en ;\n",
        members, "\n",
        "    .\n"
      )
    })

  # Write as a separate Turtle file
  cluster_ttl_path <- here("data", "processed", "rdf", "clusters.ttl")
  header <- "PREFIX skos: <http://www.w3.org/2004/02/skos/core#>\n"
  write(paste(c(header, cluster_turtle), collapse = "\n\n"), cluster_ttl_path)

  message(sprintf(
    "06_create_scheme.R: %d cluster collections written to %s",
    length(cluster_turtle), cluster_ttl_path
  ))

} else {
  message("06_create_scheme.R: non_structure_clusters.rds not found — skipping cluster TTL. Run 03_analysis.R first.")
}



# ==============================================================================
# OA Annotations: substance → source tagging
# ==============================================================================

annotation_prefix <- "https://data.omgeving.vlaanderen.be/id/annotation/chemical_substance/"
concept_prefix_ann <- "https://data.omgeving.vlaanderen.be/id/concept/chemical_substance/"

annotation_blocks <- all_substances |>
  distinct(hash, source) |>
  filter(!is.na(hash), !is.na(source)) |>
  purrr::pmap_chr(function(hash, source) {
    source_slug <- gsub("[^a-zA-Z0-9_-]", "_", source)
    ann_uri <- paste0("<", annotation_prefix, hash, "_", source_slug, ">")
    paste0(
      ann_uri, "\n",
      "    rdf:type        oa:Annotation ;\n",
      "    dc:source       \"", source, "\" ;\n",
      "    skos:prefLabel  \"", source, "\" ;\n",
      "    oa:hasBody      <", concept_prefix_ann, source_slug, "> ;\n",
      "    oa:hasTarget    <", concept_prefix_ann, hash, "> ;\n",
      "    oa:motivatedBy  oa:tagging .\n"
    )
  })

source_concept_blocks <- all_substances |>
  distinct(source) |>
  filter(!is.na(source)) |>
  purrr::pmap_chr(function(source) {
    source_slug <- gsub("[^a-zA-Z0-9_-]", "_", source)
    paste0(
      "<", concept_prefix_ann, source_slug, ">\n",
      "    rdf:type        skos:Concept ;\n",
      "    skos:prefLabel  \"", source, "\" .\n"
    )
  })

annotation_ttl_path <- here("data", "processed", "rdf", "annotations.ttl")
ann_header <- paste(
  "PREFIX rdf:        <http://www.w3.org/1999/02/22-rdf-syntax-ns#>",
  "PREFIX dc:         <http://purl.org/dc/elements/1.1/>",
  "PREFIX skos:       <http://www.w3.org/2004/02/skos/core#>",
  "PREFIX oa:         <http://www.w3.org/ns/oa#>",
  "PREFIX annotation: <https://data.omgeving.vlaanderen.be/id/annotation/chemical_substance/>",
  "PREFIX csc:        <https://data.omgeving.vlaanderen.be/id/concept/chemical_substance/>",
  "",
  sep = "\n"
)

write(
  paste(c(ann_header, source_concept_blocks, annotation_blocks), collapse = "\n\n"),
  annotation_ttl_path
)

message(sprintf(
  "06_create_scheme.R: %d source concepts + %d annotations written to %s",
  length(source_concept_blocks), length(annotation_blocks), annotation_ttl_path
))


system("bash bash/merge.bash")

# ==============================================================================
# add chebi
# ==============================================================================

system("bash bash/chebi.sh")

chebi <- read.csv(here("data", "processed", "chebi.csv"), stringsAsFactors = FALSE)


