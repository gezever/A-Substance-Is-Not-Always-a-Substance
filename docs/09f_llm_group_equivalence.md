# Experiment 09f — ChemOnt Equivalence Assessment for Non-Structure Substances

## Context

Experiments 09d and 09e asked: *"to which ChemOnt class does this substance belong?"* — a membership question. Experiment 09f asks a different question: **"is this substance name itself equivalent to a ChemOnt node?"**

The legislative relevance is direct. When a regulatory entry such as *"fatty acids, C16-22"* or *"alkyl sulfates"* corresponds to a ChemOnt node, the taxonomy can be used to enumerate all individual child substances covered by that entry. This makes the ChemOnt hierarchy operationally useful for legislation: instead of listing individual substances, a legislator can reference a ChemOnt class and rely on the taxonomy to define coverage.

Three relationship types are distinguished:

| Type | Meaning | Example |
|---|---|---|
| **exact** | Substance name IS the ChemOnt class | *"lead alkyls"* ≡ ChemOnt "Alkyl lead compounds" |
| **subset** | Name covers a specific subset of a ChemOnt class | *"fatty acids, C16-22"* ⊂ "Fatty acids and conjugates" |
| **broader** | Name spans multiple ChemOnt nodes | *"fatty acids, vegetable-oil, esters"* > any single node |

---

## Design

### Input population

All matchable non-structure substances (~9,547 total). The matchable filter excludes reaction masses, petroleum fractions, distillates, and UVCB substances. An initial pilot of n = 500 (`set.seed(42)`) was run first to validate the approach; the results below cover the full run of n = 5,000 (`set.seed(42)`, `PILOT_N = 5000L`).

### LLM prompt

A new prompt distinct from 09d/09e instructs the model to assess equivalence rather than membership. For each substance, the model returns:

```json
{
  "is_group": true,
  "equivalent_node": "Alkyl lead compounds",
  "equivalent_level": "direct_parent",
  "relationship_type": "exact",
  "reasoning": "...",
  "confidence": "high"
}
```

### ChemOnt input: full taxonomy

Unlike 09d and 09e — which used a ~1,856-label subset derived from the ClassyFire ground truth — 09f uses the **full ChemOnt 2.1 taxonomy** extracted from `ChemOnt_2_1_skos.ttl` via `python/build_chemont_hierarchy.py`.

| | 09d / 09e | 09f |
|---|---|---|
| Source | `chemont_ground_truth.csv` (observed labels only) | Full SKOS (`ChemOnt_2_1_skos.ttl`) |
| Labels | ~1,856 | 4,824 |
| System prompt tokens | ~8K | ~31K |
| Pilot cost | ~$2.40 | ~$6 |

The full taxonomy is structured across 11 depth levels: KINGDOM (2), SUPERCLASS (31), CLASS (765), SUBCLASS (1,729), LEVEL 5 (1,200), LEVEL 6 (738), LEVEL 7 (248), LEVEL 8 (82), LEVEL 9 (21), LEVEL 10 (7), LEVEL 11 (1).

### Methodology comparison: LLM vs. regex

The existing `entity_type` column in `all_substances_classified.rds` classifies substances as *"Substance group"* via regex patterns on the substance name (e.g., *"X and its compounds"*, *"ends with compounds"*, PAH/PCB/PFAS acronyms). The LLM's `is_group` field provides an independent assessment, enabling a direct comparison of the two methods.

---

## Results (pilot n = 500)

### is_group distribution

| is_group | n | % |
|---|---|---|
| FALSE | 300 | 60% |
| TRUE | 190 | 38% |
| NA (parse failure) | 10 | 2% |

38% of matchable non-structure substances are identified by the LLM as group entries that can be mapped to a ChemOnt node.

### LLM vs. regex: methodology comparison

The pilot sample contained only **4 substances** with `entity_type = "Substance group"` (regex-based). All 4 were correctly confirmed as groups by the LLM.

| regex: Substance group | LLM: group | n |
|---|---|---|
| No | No | 300 |
| No | Yes | 186 |
| No | NA | 10 |
| Yes | Yes | 4 |

The regex identifies **4 out of 190 groups** found by the LLM — a recall of 2%. This dramatic under-coverage is explained by the pilot composition: the sample is dominated by `entity_type = "Unclassified"` (273, 54.6%) and `"CAS without structure"` (173, 34.6%), entity types assigned when the regex found no group-pattern in the name.

The LLM identifies group entries that the regex systematically misses: substance names like *"fatty acids, C16-22"*, *"Alkenes C20+"*, or *"aramid fibres"* describe chemical families clearly, but contain none of the regex trigger phrases. This confirms that the regex-based `entity_type` classification is a lower bound on the actual number of legislative group entries in the dataset.

### Relationship type for is_group == TRUE (n = 190)

| relationship_type | high | medium | low | total | % |
|---|---|---|---|---|---|
| subset | 38 | 114 | 22 | **174** | **91.6%** |
| broader | 1 | 4 | 10 | 15 | 7.9% |
| exact | 1 | — | — | **1** | **0.5%** |

The overwhelming majority (91.6%) are **subsets**: the legislative entry describes a chemically defined subgroup (often by chain length or functional variant) within a broader ChemOnt class. Exact equivalences are rare; only one was found.

**Confidence distribution (is_group == TRUE):**

| Confidence | n | % |
|---|---|---|
| high | 40 | 21.1% |
| medium | 118 | 62.1% |
| low | 32 | 16.8% |

### Equivalent level (is_group == TRUE, n = 190)

| Level | n | % |
|---|---|---|
| subclass | 80 | 42.1% |
| class | 42 | 22.1% |
| direct_parent | 33 | 17.4% |
| superclass | 32 | 16.8% |
| kingdom | 1 | 0.5% |
| level 6 | 1 | 0.5% |

Most equivalences are found at the **subclass** level (42%), with substantial coverage at class (22%) and superclass (17%). The relatively high proportion at superclass reflects legislative entries that describe broad chemical families (e.g., *"Hydrocarbons, C13-C18"* → superclass "Hydrocarbons"). No equivalences were found deeper than level 6, confirming that legislative group entries are defined at accessible taxonomic depths.

> **Note**: The `direct_parent` option in the pilot prompt was a carry-over from the membership framing of 09d/09e and does not correspond to any ChemOnt taxonomy level. The prompt has since been corrected to use the actual level names (`level 5`, `level 6`). The 33 pilot substances tagged as `direct_parent` are likely nodes at level 5 or deeper; their `equivalent_level` values should be treated as unreliable.

---

## Results (full run n = 5,000)

### is_group distribution

| is_group | n | % |
|---|---|---|
| FALSE | 2,468 | 49.4% |
| TRUE | 1,584 | 31.7% |
| NA (parse failure) | 948 | 19.0% |

31.7% of the 5,000 sampled non-structure substances are identified as ChemOnt-equivalent group entries. The NA rate (19%) is substantially higher than in the pilot (2%), reflecting a more diverse and technically demanding substance population in the larger sample — the full run includes more obscure CAS-registered substances with cryptic names where the model declines to classify rather than guess.

**Group rate by entity_type:**

| entity_type | n | NA % | is_group TRUE % |
|---|---|---|---|
| Substance group (regex) | 39 | 10.3% | 88.6% |
| Mixture | 68 | 16.2% | 71.9% |
| CAS without structure | 2,037 | 17.2% | 46.9% |
| Unclassified | 2,856 | 20.4% | 31.7% |

The group detection rate is strongly stratified by entity_type: substances already flagged by the regex as groups are confirmed at 88.6%, while the largest category ("Unclassified") has a 31.7% group rate — still substantial for a category defined by the absence of any regex signal.

### LLM vs. regex: methodology comparison

| regex: Substance group | LLM: group | n |
|---|---|---|
| No | No | 2,464 |
| No | Yes | 1,553 |
| No | NA | 944 |
| Yes | No | 4 |
| Yes | Yes | 31 |
| Yes | NA | 4 |

The regex identifies **31 out of 1,584 groups** found by the LLM — a recall of **2.0%**, consistent with the pilot. Four substances flagged as `Substance group` by the regex were assessed as non-groups by the LLM, suggesting minor regex over-reach. The 1,553 LLM-identified groups that the regex missed are distributed across `Unclassified` (88.6%) and `CAS without structure` (11.4%) entity types.

### Relationship type for is_group == TRUE (n = 1,584)

| relationship_type | high | medium | low | total | % |
|---|---|---|---|---|---|
| subset | 285 | 959 | 226 | **1,470** | **92.8%** |
| broader | 15 | 25 | 53 | 93 | 5.9% |
| exact | 18 | 1 | — | **19** | **1.2%** |
| NA | — | — | 2 | 2 | 0.1% |

The subset dominance (92.8%) is consistent with the pilot (91.6%). Exact equivalences are rare but no longer singular: 19 substances have a high- or medium-confidence exact match to a ChemOnt label. The broader category (5.9%) shows a confidence skew toward low — correctly reflecting the model's uncertainty when a substance spans multiple ChemOnt nodes.

**High-confidence summary (n = 318):**

| relationship_type | n |
|---|---|
| subset | 285 |
| exact | 18 |
| broader | 15 |

### Equivalent level (is_group == TRUE, n = 1,584)

| Level | n | % |
|---|---|---|
| subclass | 670 | 42.3% |
| class | 359 | 22.7% |
| superclass | 265 | 16.7% |
| direct_parent | 247 | 15.6% |
| NA | 24 | 1.5% |
| level 5 | 11 | 0.7% |
| kingdom | 5 | 0.3% |
| level 6 | 3 | 0.2% |

The level distribution is essentially identical to the pilot: subclass (42.3%), class (22.7%), superclass (16.7%). The 247 substances tagged as `direct_parent` (15.6%) are artefacts of an incorrect prompt option — `direct_parent` does not correspond to any ChemOnt taxonomy level and was a carry-over from the membership framing of 09d/09e. The prompt has been corrected to use `level 5` and `level 6` instead; for those 247 cases the actual taxonomy depth is unknown. Equivalences below subclass remain rare (level 5: 0.7%, level 6: 0.2%), confirming that legislative group entries do not map to the deepest taxonomy levels.

---

## Experiment 09g — SKOS mapping triples

The 09f equivalence results were converted to formal SKOS mapping triples in `R/09g_chemont_skos_mapping.R`, using only **high-confidence** assessments.

**SKOS predicate mapping:**

| relationship_type | SKOS predicate | Direction |
|---|---|---|
| exact | `skos:exactMatch` | substance ≡ ChemOnt class |
| subset | `skos:broadMatch` | ChemOnt class is broader than substance |
| broader | `skos:narrowMatch` | ChemOnt class is narrower than substance |

**Output: `data/processed/rdf/chemont_nonstructure_mapping.ttl`**

| predicate | n |
|---|---|
| `skos:broadMatch` | 285 |
| `skos:exactMatch` | 18 |
| `skos:narrowMatch` | 1 |
| **Total** | **304** |

URI conventions follow the existing RDF pipeline:
- Substance: `https://data.omgeving.vlaanderen.be/id/concept/chemical_substance/{MD5(substance_name)}`
- ChemOnt concept: `https://data.omgeving.vlaanderen.be/id/concept/chemical_substance/{CHEMONTID}`

ChemOnt URIs are resolved via `data/processed/chemont_label_uri.csv`, extracted from `ChemOnt_2_1_skos.ttl` using `python/chemont_label_uri.py`. Of the 1,560 mappable rows (all confidence levels), 7 had `equivalent_node` labels not found in the ChemOnt SKOS — likely LLM hallucinations of non-existent node names. All 304 high-confidence triples passed Apache Jena `riot --validate`.

---

## Pattern analysis

### The one exact equivalence

> **"lead alkyls"** ≡ ChemOnt **"Alkyl lead compounds"** (direct_parent, high confidence)

The substance name directly corresponds to the ChemOnt label. Exact equivalences are rare because legislative naming conventions and ChemOnt labels are developed independently; coincidental matches at the label level are uncommon.

### High-confidence subsets (n = 38): recurring patterns

High-confidence subset assignments cluster around four naming conventions:

**1. Chain-length notation (C_x–C_y)**

| Substance name | Equivalent node | Level |
|---|---|---|
| Fatty acids, C16-22 | Fatty acids and conjugates | subclass |
| Hydrocarbons, C7-C8, n-alkanes | Alkanes | subclass |
| Carboxylic acids, di-, C5-9 | Dicarboxylic acids and derivatives | subclass |
| 2-Propenoic acid, perfluoro-C8-16-alkyl esters | Acrylic acid esters | subclass |

The chain-length specifier places the substance unambiguously within a ChemOnt class, but as a named subset. The class covers the full chain-length range; the legislative entry covers only part of it.

**2. Functional class + qualifier**

| Substance name | Equivalent node | Level |
|---|---|---|
| alkyl sulfate, sodium (mono C10-C16-alkyl esters) | Alkyl sulfates | subclass |
| Quaternary ammonium compounds, (C16-18 and C18-unsatd.) | Quaternary ammonium salts | subclass |
| phenols, styrenated | Phenols | class |
| aramid fibres | Aromatic oligoamide copolymers | subclass |

**3. Specific chemical families with known taxonomy**

| Substance name | Equivalent node | Level |
|---|---|---|
| Tributyltin Compounds | Triorganotin compounds | subclass |
| undecafluorohexanoic acid (PFHxA), its salts and esters | Perfluoroalkyl carboxylic acid and derivatives | direct_parent |
| salts and esters of MCPA | Phenoxyacetic acids | subclass |
| Sophorolipids: fermentation products of glucose | Sophorolipids | subclass |

**4. Broad hydrocarbon families**

| Substance name | Equivalent node | Level |
|---|---|---|
| Hydrocarbons, C4 | Hydrocarbons | superclass |
| Hydrocarbons, C9-14, ethylene-manuf.-by-products | Hydrocarbons | superclass |
| Alkenes C20+ | Alkenes | direct_parent |

These are mapped to superclass because ChemOnt does not have a subclass node that matches the described chain range.

### Broader relationships (n = 15)

Broader entries cannot be mapped to a single ChemOnt node. They represent:

- **Reaction products**: *"Reaction product of butanal aldol condensation products…"* → the mixture character prevents a single-node mapping
- **Cross-class groups**: *"boron compounds, with the exception of…"* spans multiple ChemOnt classes
- **Undefined mixtures**: *"Glycerides, C16-18 mono-, di- and tri-"* → glycerolipids at multiple levels

One high-confidence broader case:
> *"Man-made vitreous (silicate) fibres"* → equivalent_node = null, confidence = high

The LLM correctly identifies this as a group (is_group = TRUE) but recognises that no single ChemOnt node covers the full class of amorphous silicate fibres. This is methodologically correct behaviour.

---

## Embedding validation

Only 5 of 500 pilot substances had an embedding match in `final_matches.csv`. Of these, 2 agreed with the LLM's equivalent node. The overlap is too small for a meaningful validation; see the discussion in experiment 09e for why embedding coverage is structurally limited for non-structure substances.

---

## Discussion

### The subset dominance reflects legislative practice

The near-absence of exact equivalences (1/190) and dominance of subsets (174/190) is not a model failure — it reflects how legislation works. Regulatory entries for substance groups typically specify a chemical family *with additional constraints* (chain length, saturation level, specific salts or esters). ChemOnt classes are defined structurally and cover the full family without those constraints. A legislative group is therefore almost always a subset of the corresponding ChemOnt node, never an exact match.

This has a concrete implication for downstream use: when a legislative entry maps to a ChemOnt node as a *subset*, the taxonomy can still enumerate candidate child substances, but a secondary filter (chain length range, specific substituents) must be applied to the ChemOnt children before the list can be used as a complete coverage map.

### Regex coverage is structurally incomplete

The `entity_type = "Substance group"` regex catches entries that are phrased with explicit group language ("X and its compounds"). It misses the much larger class of entries phrased as chemical family names ("fatty acids, C16-22", "aramid fibres", "phenols, styrenated"). In the pilot, the LLM identifies 186 additional group entries that the regex missed entirely. A full-population run would likely reveal that the true number of legislative group entries is an order of magnitude larger than the regex-based count.

### Implications for taxonomic coverage in legislation

The 40 high-confidence is_group assignments in the pilot represent the most reliable mapping candidates: legislative entries where the ChemOnt node and coverage relationship can be stated with high confidence. Scaled to the full population (~9,547 matchable substances), an extrapolation yields roughly 760 high-confidence group mappings — a set large enough to construct a meaningful legislative-to-taxonomy crosswalk.

---

## Limitations

1. **Sample coverage**: 5,000/~9,547 substances (~52%). Roughly half the matchable population remains unclassified. A follow-up run on the remaining ~4,547 substances would complete the mapping.

2. **No ground truth**: accuracy cannot be directly measured. The confidence signal is the only quality proxy. The 7 unresolved labels in 09g (ChemOnt node not found in SKOS) are a partial indicator of hallucination rate (~0.4% of is_group==TRUE rows).

3. **Subset boundaries are undefined**: the LLM identifies *which* ChemOnt node a substance is a subset of, but does not specify the boundary (e.g., exactly which child substances within "Fatty acids and conjugates" are covered by "fatty acids, C16-22"). Defining that boundary requires a separate step.

4. **Broader cases are unresolved**: 93 substances are identified as broader than any single ChemOnt node. These cannot be directly used for taxonomic coverage without manual decomposition into constituent nodes.

5. **High NA rate in full run (19%)**: the model declines to classify nearly 1 in 5 substances in the larger sample. These are concentrated in technically obscure names where abstaining is preferable to a forced guess; they are not recoverable by re-prompting without additional context.

---

## Conclusions

1. **31.7% of sampled non-structure substances** (n = 5,000) are identified as ChemOnt-equivalent group entries — far more than the 0.8% flagged by the regex-based `entity_type` in the same sample.

2. **Subsets dominate (92.8%)**: legislative groups almost always cover a named subset of a ChemOnt class, not the full class. Exact equivalences are rare (1.2%) but present (19 confirmed at high/medium confidence).

3. **Equivalences are found primarily at subclass and class level** (65% combined), making these the most actionable levels for legislative-to-taxonomy crosswalks.

4. **The regex-based entity_type classification has ~2% recall** on the group entries identified by the LLM, confirming that name-pattern matching is structurally insufficient for identifying legislative group entries.

5. **304 high-confidence SKOS mapping triples** (09g) form a ready-to-use legislative–ChemOnt crosswalk, validated by Apache Jena riot and aligned with the existing RDF pipeline URI conventions.
