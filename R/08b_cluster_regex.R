# ==============================================================================
# 08b_cluster_regex.R
# Analysis 8b: Reverse-engineering regular expressions for embedding clusters
#
# PURPOSE
# -------
# The six embedding clusters produced by 08_embedding_clustering.R were
# labelled manually based on substantive inspection of cluster contents.
# This script asks whether those labels can be operationalised as explicit
# regular expressions — i.e., whether the lexical surface of substance names
# contains enough systematic signal to discriminate clusters without embeddings.
#
# The question is relevant for two practical reasons:
#   1. Interpretability: a regex makes the cluster criterion explicit and
#      falsifiable, whereas an embedding-based cluster is a black box.
#   2. Prospective classification: new substance names can be tested against
#      the regexes without recomputing embeddings.
#
# DATA PROVENANCE
# ---------------
# Input:  `data/processed/non_structure_clusters.csv`
#         Produced by `08_embedding_clustering.R`.  Contains `substance_name`,
#         `cluster` (integer 1–6), `umap1`, `umap2`, and `manual_label`.
#
# METHODOLOGY
# -----------
# For each cluster, a candidate regex was hand-crafted after inspecting random
# samples (n = 50) and systematic frequency counts of candidate patterns.
# Candidate patterns were chosen to target the most distinctive lexical
# features of each cluster:
#
#   Cluster 1  Reaction masses
#              Metal-element group names (Lead, Nickel, etc.), "and its
#              salts/compounds", PFAS abbreviations, ceramic fibres.
#
#   Cluster 2  Complex substances (UVCB with structure)
#              "(petroleum)" in the name, ECHA long-form descriptions starting
#              with "[A complex combination", boiling-point specifications.
#              This cluster is by far the most lexically distinctive.
#
#   Cluster 3  Petroleum & coal tar fractions
#              IUPAC locant prefixes (^[0-9]+,[0-9] / ^[0-9]+-[a-z]),
#              stereo-descriptor parentheses (^([0-9] / ^([±]),
#              "covering any … isomers", polyhalogenated aromatics (diphenyl
#              ethers, naphthalenes, cyclododecane), polychlorinated dibenzo
#              compounds, organotins, diisocyanates.
#
#   Cluster 4  Inorganic compounds & metal salts
#              Anthracene oil/paste/fraction, nonylphenol, hexachloro-
#              cyclohexane, C-range ethoxylated alcohols, formaldehyde
#              releasers, azo colorants, endosulfan, pentachlorophenol,
#              chlorinated paraffins, tributyltin, specific benzenedicarboxylic
#              acid esters.
#
#   Cluster 5  Trade names, codes & biological materials
#              "Reaction mass of:" (with colon — more distinctive than without),
#              trade code patterns (3–15 upper-case letters + digit suffix),
#              short all-caps codes (^[A-Z]{2,5}-[0-9]{3,}), pigment/dye codes.
#
#   Cluster 6  Polymers, fatty acids & surfactants
#              Microbial biocontrol agents (Bacillus, Beauveria, Aspergillus …),
#              natural waxes and gums (beeswax, carnauba, xanthan gum …),
#              structural proteins (gelatin, casein, albumin),
#              strain identifiers (DSM, ATCC, FMCH, FZB …).
#
# Performance is evaluated via precision, recall, and F1 against the cluster
# labels.  Because the clusters were formed by semantic embedding similarity
# rather than explicit lexical rules, recall is inherently limited — substances
# within a cluster share chemical meaning, not necessarily surface-level text
# patterns.  Precision is therefore the more informative metric: a high-
# precision regex with modest recall is still useful for prospective
# classification and for communicating what is *most typical* of each cluster.
#
# LIMITATIONS
# -----------
# - Recall is low for all clusters except Cluster 2 (~12 % for clusters 1/3/5,
#   ~2–6 % for clusters 4/6).  A single regex cannot exhaustively represent an
#   embedding-based grouping.
# - Clusters 4 and 5 show the lowest F1 scores because they function as
#   catch-all categories for heterogeneous substance names that did not fit the
#   other clusters.
# - The regexes are designed for human interpretability, not for maximum
#   performance; a machine-learning classifier would outperform them but would
#   not be self-explanatory.
#
# OUTPUTS
# -------
# output/figures/Analysis_8b_cluster_regex_precision_recall.pdf
# output/tables/Analysis_8b_cluster_regex_performance.csv
# (cluster_regexes: named character vector, available in R session)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(here)

# ------------------------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------------------------

non_structure_clusters <- read.csv(
  here("data", "processed", "non_structure_clusters.csv"),
  stringsAsFactors = FALSE
)

# ==============================================================================
# Analysis 8b-i: Define representative regular expressions per cluster
# ==============================================================================
# INTENT
# Hand-crafted regexes that target the most characteristic and diagnostic lexical
# patterns in each cluster.  All patterns use (?i) for case-insensitive matching.
# Patterns are combined with | within each cluster regex.
# ==============================================================================

cluster_regexes <- c(

  # --------------------------------------------------------------------------
  # Cluster 1 — Reaction masses
  # Targets: metal-element group names at the start of the name; the phrase
  # "and its salts/compounds"; PFAS abbreviations; refractory ceramic fibres;
  # medium-chain chlorinated paraffin groups; names ending in "fibres".
  # --------------------------------------------------------------------------
  "1" = paste0(
    "(?i)(",
    "^(Lead|Mercury|Arsenic|Nickel|Cadmium|Chromium|Cobalt|Antimony|Bismuth|",
    "Barium|Thallium|Copper|Selenium)\\b",
    "|\\band its (salts?|compounds?)( and|$|,)",
    "|\\bPF(AS|BS|DA|HxA|HxS|NA|UnDA)\\b",
    "|Refractory Ceramic Fibres?",
    "|Medium.chain chlorinated",
    "|\\bfibres?$",
    "|Polychlorin.*(paraffin|terphenyl|biphenyl)",
    ")"
  ),

  # --------------------------------------------------------------------------
  # Cluster 2 — Complex substances (UVCB with structure)
  # Targets: "(petroleum)" anywhere in the name; ECHA long-form description
  # marker "[A complex combination"; boiling-point or boiling-range
  # specifications.  These three patterns together cover 58 % of the cluster
  # with 85 % precision — by far the most lexically distinctive cluster.
  # --------------------------------------------------------------------------
  "2" = paste0(
    "(?i)(",
    "\\(petroleum\\)",
    "|\\[A complex combination",
    "|boiling (point|range)",
    ")"
  ),

  # --------------------------------------------------------------------------
  # Cluster 3 — Petroleum & coal tar fractions
  # Targets: IUPAC locant prefixes (^2,3- / ^4-amino etc.); stereo-descriptor
  # parentheses (^([0-9] or ^([±])); "covering any … isomers"; poly-halogenated
  # aromatics (diphenyl ethers, naphthalenes, cyclododecane); polychlorinated
  # dibenzo-dioxins/furans/naphthalenes; organotins; diisocyanates;
  # chlorinated alkane entries starting with "Alkanes, Cxx-yy, chloro".
  # --------------------------------------------------------------------------
  "3" = paste0(
    "(?i)(",
    "^[0-9]+,[0-9]",
    "|^[0-9]+-[a-z]",
    "|^\\([0-9]",
    "|^\\([±]",
    "|covering any.*isomers",
    "|^(Poly|Hexa|Penta|Tetra|Octa)(brom|chlor).*(ether|diphenyl|naphthal|cyclododec)",
    "|Polychlorin.*(dibenzo|dioxin|naphthal)",
    "|\\bdiisocyanate\\b",
    "|\\borganotin\\b",
    "|\\bstannan\\b",
    "|^Alkanes?, C[0-9]+-[0-9]+, chloro",
    "|^Polychlorinated naphth",
    ")"
  ),

  # --------------------------------------------------------------------------
  # Cluster 4 — Inorganic compounds & metal salts
  # Targets: anthracene oil/paste/fraction products; nonylphenol variants;
  # hexachlorocyclohexane (lindane group); C-range ethoxylated or branched
  # alcohols; formaldehyde releasers; azo colorants; endosulfan;
  # pentachlorophenol; chlorinated paraffins; tributyltin compounds;
  # 1,2-benzenedicarboxylic acid di-C[n] esters.
  # --------------------------------------------------------------------------
  "4" = paste0(
    "(?i)(",
    "^Anthracene (oil|paste|fraction)",
    "|\\bnonylphenol\\b",
    "|hexachlorocyclohexane",
    "|Alcohols?, C[0-9]+-[0-9]+.*(ethoxylated|branched)",
    "|formaldehyde.*releas",
    "|Azocolour",
    "|endosulfan",
    "|pentachlorophenol",
    "|chlorinated paraffins",
    "|^4-Nonylphenol",
    "|tributyltin compounds",
    "|^1,2-Benzenedicarboxylic acid.*di-C[0-9]",
    ")"
  ),

  # --------------------------------------------------------------------------
  # Cluster 5 — Trade names, codes & biological materials
  # Targets: "Reaction mass of:" (with colon, more distinctive than without);
  # trade code pattern: 3–15 upper-case letters followed by a dash/space and
  # digits+letters (e.g. UC-706, HYAFF 11, POLYFOX PF-3320); short codes
  # (^[A-Z]{2-5}-[0-9]{3+}); pigment and dye codes (PIGMENT YELLOW …,
  # DISPERSE BLUE …, REACTIVE RED …).
  # --------------------------------------------------------------------------
  "5" = paste0(
    "(?i)(",
    "^[Rr]eaction mass of:",
    "|^[A-Z]{3,15}[- ][0-9][A-Z0-9.-]{1,12}$",
    "|^[A-Z]{2,5}-[0-9]{3,}$",
    "|^PIGMENT ",
    "|^DISPERSE ",
    "|^REACTIVE (BLUE|RED|YELLOW|BLACK)",
    "|^[A-Z]{4,}[_ ][A-Z]{2,}[- _][0-9A-Z]{2,}$",
    ")"
  ),

  # --------------------------------------------------------------------------
  # Cluster 6 — Polymers, fatty acids & surfactants
  # Targets: microbial biocontrol agents (Bacillus, Beauveria, Aspergillus, …)
  # identified by genus name at the start; natural waxes (beeswax, carnauba,
  # montan wax, ceresin) and gums/resins (xanthan gum, guar gum, tragacanth,
  # gum arabic, dammar, rosin); structural proteins (gelatin, casein, albumin);
  # strain identifiers (DSM, ATCC, FMCH, FZB, …) co-occurring with "strain".
  # --------------------------------------------------------------------------
  "6" = paste0(
    "(?i)(",
    "^(Bacillus|Beauveria|Aspergillus|Trichoderma|Metarhizium|Isaria|",
    "Paecilomyces|Clonostachys|Phlebiopsis|Aureobasidium|Adoxophyes|",
    "Ampelomyces|Steinernema|Pythium|Lecanicillium)\\b",
    "|\\b(beeswax|carnauba|montan wax|ceresin|gelatin|casein|albumin|rosin|",
    "xanthan gum|guar gum|tragacanth|gum arabic|dammar)\\b",
    "|\\bstrain(s)?\\b.*(DSM|ATCC|IAB|QST|FMCH|FZB|BV-|GC-|MBI|PTA-)",
    ")"
  )
)

# ==============================================================================
# Analysis 8b-ii: Evaluate precision, recall, and F1 per cluster
# ==============================================================================
# INTENT
# Quantify how well each regex discriminates its target cluster from all others.
# Precision = TP / (TP + FP): fraction of regex matches that are truly in the
#   cluster.  High precision means few false positives.
# Recall    = TP / N: fraction of cluster members that the regex captures.
#   Low recall is expected given the embedding-based cluster formation.
# F1        = harmonic mean of precision and recall.
# ==============================================================================

names_vec <- non_structure_clusters$substance_name
cluster_vec <- non_structure_clusters$cluster

regex_performance <- lapply(names(cluster_regexes), function(cl_id) {
  pat  <- cluster_regexes[cl_id]
  cl_n <- as.integer(cl_id)

  matches  <- grepl(pat, names_vec, perl = TRUE)
  tp       <- sum(matches & cluster_vec == cl_n)
  fp       <- sum(matches & cluster_vec != cl_n)
  fn       <- sum(!matches & cluster_vec == cl_n)
  N        <- tp + fn

  precision <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  recall    <- if (N > 0)       tp / N          else NA_real_
  f1        <- if (!is.na(precision) & !is.na(recall) & (precision + recall) > 0)
                 2 * precision * recall / (precision + recall) else NA_real_

  label <- unique(non_structure_clusters$manual_label[cluster_vec == cl_n])

  data.frame(
    cluster   = cl_n,
    label     = label,
    N         = N,
    TP        = tp,
    FP        = fp,
    recall    = recall,
    precision = precision,
    F1        = f1,
    stringsAsFactors = FALSE
  )
}) |> bind_rows()

print(regex_performance |>
  mutate(
    recall    = sprintf("%.0f%%", 100 * recall),
    precision = sprintf("%.0f%%", 100 * precision),
    F1        = round(F1, 2)
  ) |>
  select(cluster, label, N, TP, FP, recall, precision, F1))

write.csv(
  regex_performance,
  here("output", "tables", "Analysis_8b_cluster_regex_performance.csv"),
  row.names = FALSE
)

# ==============================================================================
# Analysis 8b-iii: Visualise precision–recall per cluster
# ==============================================================================
# INTENT
# Plot precision and recall as a paired dot-plot so that the trade-off is
# immediately visible.  Each cluster is one row; the two metrics are shown as
# separate dots connected by a line segment.  Cluster 2 should stand out
# clearly as the most discriminable cluster.
# ==============================================================================

pr_long <- regex_performance |>
  select(cluster, label, precision, recall) |>
  pivot_longer(cols = c(precision, recall),
               names_to  = "metric",
               values_to = "value") |>
  mutate(
    label  = factor(
      paste0("Cluster ", cluster, ": ", label),
      levels = paste0("Cluster ", 6:1, ": ",
                      regex_performance$label[order(regex_performance$cluster,
                                                    decreasing = TRUE)])
    ),
    metric = factor(metric, levels = c("recall", "precision"))
  )

p8b <- ggplot(pr_long, aes(x = value, y = label, colour = metric, shape = metric)) +
  geom_line(
    aes(group = label),
    colour = "grey70", linewidth = 0.6
  ) +
  geom_point(size = 3.5) +
  geom_text(
    aes(label = sprintf("%.0f%%", 100 * value)),
    nudge_y = 0.35,
    size     = 3,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(add = 0.05)
  ) +
  scale_colour_manual(
    values = c(recall = "#4a90d9", precision = "#e05c5c"),
    labels = c(recall = "Recall (TP / N)", precision = "Precision (TP / (TP+FP))")
  ) +
  scale_shape_manual(
    values = c(recall = 16, precision = 17),
    labels = c(recall = "Recall (TP / N)", precision = "Precision (TP / (TP+FP))")
  ) +
  labs(
    title    = "Regex performance per cluster",
    subtitle = paste0(
      "Cluster 2 (UVCB) is the most lexically distinctive. ",
      "Low recall in other clusters reflects embedding-based (not rule-based) cluster formation."
    ),
    x        = NULL,
    y        = NULL,
    colour   = "Metric",
    shape    = "Metric"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    panel.grid.minor = element_blank()
  )

print(p8b)
ggsave(
  p8b,
  filename = here("output", "figures",
                  "Analysis_8b_cluster_regex_precision_recall.pdf"),
  device   = "pdf",
  height   = 5, width = 9, units = "in"
)

# ==============================================================================
# Analysis 8b-iv: Cross-cluster heatmap of false-positive rates
# ==============================================================================
# INTENT
# Show which clusters are confused by each regex: cell (i, j) is the number of
# names from cluster j that are matched by cluster i's regex.  The diagonal
# gives TP; off-diagonal cells give FP broken down by source cluster.
# A sparse off-diagonal confirms that regexes have good specificity despite
# low recall.
# ==============================================================================

confusion_mat <- outer(
  names(cluster_regexes),
  sort(unique(cluster_vec)),
  FUN = Vectorize(function(regex_id, true_cl) {
    pat  <- cluster_regexes[regex_id]
    sum(grepl(pat, names_vec[cluster_vec == true_cl], perl = TRUE))
  })
)
rownames(confusion_mat) <- paste0("Regex cl", names(cluster_regexes))
colnames(confusion_mat) <- paste0("True cl", sort(unique(cluster_vec)))

confusion_long <- as.data.frame(as.table(confusion_mat)) |>
  rename(regex_cluster = Var1, true_cluster = Var2, count = Freq) |>
  mutate(
    diagonal = sub("Regex cl", "", regex_cluster) ==
               sub("True cl",  "", true_cluster)
  )

p8b_heat <- ggplot(
  confusion_long,
  aes(x = true_cluster, y = regex_cluster, fill = count)
) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(
    aes(label = count, colour = diagonal),
    size = 3.5, fontface = "bold"
  ) +
  scale_fill_gradient(low = "white", high = "#4a90d9",
                      name = "Matches (n)") +
  scale_colour_manual(values = c("FALSE" = "grey40", "TRUE" = "#c0392b"),
                      guide = "none") +
  labs(
    title    = "Regex cross-cluster match counts",
    subtitle = "Diagonal (red) = true positives; off-diagonal = false positives by source cluster",
    x        = "True cluster (name origin)",
    y        = "Applied regex"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x   = element_text(angle = 30, hjust = 1),
    plot.subtitle = element_text(colour = "grey40", size = 9),
    panel.grid    = element_blank()
  )

print(p8b_heat)
ggsave(
  p8b_heat,
  filename = here("output", "figures",
                  "Analysis_8b_cluster_regex_confusion.pdf"),
  device   = "pdf",
  height   = 5, width = 7, units = "in"
)
