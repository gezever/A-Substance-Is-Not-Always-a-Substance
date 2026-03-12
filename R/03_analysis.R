# ==============================================================================
# 03_analysis.R
# Kernanalyse: voorkomen van chemische stoffen in requesten en referentielijsten
# ==============================================================================

library(dplyr)
library(tidyr)
library(here)

# ------------------------------------------------------------------------------
# Inladen schone data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))


# ------------------------------------------------------------------------------
# Analyse 1: Frequentie van stoffen in requesten
# ------------------------------------------------------------------------------

substance_freq <- all_substances |>
  count(cas_number, source, sort = TRUE) |>
  rename(n_requests = n)

# ------------------------------------------------------------------------------
# Analyse 2: Aanwezigheid op referentielijsten per stof
# ------------------------------------------------------------------------------

substance_list_presence <- substances_annotated |>
  filter(!is.na(list_name)) |>
  group_by(cas_number, name) |>
  summarise(
    n_lists     = n_distinct(list_name),
    lists       = paste(sort(unique(list_name)), collapse = "; "),
    .groups     = "drop"
  )

# ------------------------------------------------------------------------------
# Analyse 3: Stoffen die zowel vaak voorkomen in requesten ALS op lijsten staan
# ==============================================================================
# Kernvraag: "Is een stof altijd dezelfde stof?" — stoffen kunnen onder
# verschillende namen/CAS-nummers worden aangevraagd.

combined <- substance_freq |>
  left_join(substance_list_presence, by = c("cas_number", "name")) |>
  mutate(
    on_reference_list = !is.na(lists),
    n_lists           = replace_na(n_lists, 0L)
  ) |>
  arrange(desc(n_requests), desc(n_lists))

# ------------------------------------------------------------------------------
# Analyse 4: Overlap tussen requesten per tijdsperiode (optioneel)
# ------------------------------------------------------------------------------

requests_over_time <- requests |>
  mutate(year = format(date, "%Y")) |>
  count(year, name = "n_requests")

# ------------------------------------------------------------------------------
# Resultaten opslaan
# ------------------------------------------------------------------------------

saveRDS(substance_freq,          here("data", "processed", "substance_freq.rds"))
saveRDS(substance_list_presence, here("data", "processed", "substance_list_presence.rds"))
saveRDS(combined,                here("data", "processed", "combined.rds"))
saveRDS(requests_over_time,      here("data", "processed", "requests_over_time.rds"))

message("03_analysis.R voltooid")
