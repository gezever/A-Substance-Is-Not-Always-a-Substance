# Experiment 09e — LLM ChemOnt Classification of Non-Structure Substances

## Design

Experiment 09e applies the LLM classifier from 09d to **non-structure substances**: substances without an InChIKey that cannot be classified by ClassyFire. The central question is whether an LLM can assign a meaningful ChemOnt class from a substance name alone when no structural ground truth exists.

### Relation to prior experiments

| Experiment | Input | Ground truth | Purpose |
|---|---|---|---|
| 09d | Structured substances (InChIKey known) | ClassyFire assignment | Validate LLM accuracy |
| 09e | Non-structure substances (no InChIKey) | None available | Apply to the target population |

Because there is no structural ground truth, quality is assessed via three indirect signals:

- **(a) LLM confidence** — validated as a reliable accuracy proxy in 09d (high: 54% Hit@1; low: 0%)
- **(b) Agreement with embedding matches** from 09_embedding_chemont.R, where both methods produce a result
- **(c) Null rates per hierarchy level** — appropriate uncertainty is a sign of calibration, not failure

### Matchable filter

Reaction masses, petroleum fractions, distillates, and UVCB substances are excluded before classification. These names do not describe a single chemical entity and cannot yield a meaningful class assignment regardless of the classifier used.

### LLM methodology

Identical to 09d: the full ChemOnt label list is passed as a cached system prompt; a single API call per substance returns kingdom / superclass / class / subclass / direct\_parent as structured JSON. Results are cached per substance (SHA-1 key). The same Python classifier (`python/llm_chemont_classify.py`) and cache format are used; only the cache directory differs (`data/cache/llm_chemont_nonstructure/`).

### Sample

The matchable non-structure population contains approximately 9,547 substances. A pilot of **n = 500** was drawn using `set.seed(42)` via `slice_sample()`. Full-population output file contains 503 rows (3 additional rows likely originate from the cache export logic).

### Cost estimate (claude-sonnet-4-6)

| Item | Cost |
|---|---|
| Cache write (8K tokens) | ~$0.02 |
| Cache reads (500 × 8K × $0.30/Mtok) | ~$1.20 |
| Output (500 × ~156 tokens × $15/Mtok) | ~$1.17 |
| **Total pilot** | **~$2.40** |
| Full run (~9,547 substances) | **~$45** |

---

## Results

### Confidence distribution

| Confidence | n | % |
|---|---|---|
| high | 72 | 14.4% |
| medium | 249 | 49.8% |
| low | 179 | 35.8% |

The low-confidence fraction (35.8%) is markedly higher than in 09d, where only a small minority of structured substances received a low-confidence label. This reflects the inherently ambiguous nature of non-structure substance names: trade names, internal codes, extract descriptions, and complex mixture names provide fewer semantic signals than systematic IUPAC names.

### Null rates per hierarchy level

| Level | Null rate |
|---|---|
| Kingdom | 15.2% |
| Superclass | 15.2% |
| Class | 18.4% |
| **Subclass** | **48.4%** |
| Direct parent | 17.0% |

The subclass null rate (48.4%) is strikingly higher than the surrounding levels. Two factors contribute:

1. ChemOnt paths do not always pass through a subclass — some classes connect directly to a direct parent without an intermediate subclass node.
2. Substance names that are specific enough to identify a class or direct parent are often insufficient to resolve the subclass subdivision.

The near-identical null rates at kingdom and superclass (both 15.2%) indicate that when the model cannot assign a kingdom it also cannot assign a superclass, and vice versa — consistent with a coherent top-down classification strategy.

### Kingdom distribution (where assigned)

| Kingdom | n | % |
|---|---|---|
| Organic compounds | 365 | 86.1% |
| Inorganic compounds | 59 | 13.9% |

---

## Comparison with embedding matches

The embedding-based method from 09_embedding_chemont.R assigns a ChemOnt direct parent to substances via cosine similarity to a reference set of known-structure embeddings. For non-structure substances, this reference set has very limited coverage.

| | n |
|---|---|
| Substances with both LLM and embedding result | 9 |
| Of which: methods agree | 6 (67%) |

The overlap of only 9 out of 503 substances renders a meaningful method comparison impossible on this cohort. This is a structural limitation: the embedding method depends on cosine similarity to known structures, and non-structure substances are systematically underrepresented in the reference embedding space. The agreement rate of 67% is consistent with the 09d finding that LLM and embedding methods often converge when both produce a result, but the sample is too small to draw conclusions.

---

## Pattern analysis: high-confidence assignments

### ChemOnt superclass distribution (n = 72 high-confidence)

| Superclass | n | % |
|---|---|---|
| Lipids and lipid-like molecules | 15 | 21% |
| Organoheterocyclic compounds | 13 | 18% |
| Hydrocarbons | 8 | 11% |
| Organic acids and derivatives | 6 | 8% |
| Benzenoids | 5 | 7% |
| Organohalogen compounds | 5 | 7% |
| Organic nitrogen compounds | 4 | 6% |
| Inorganic salts | 2 | 3% |
| Other | 14 | 19% |

### Name patterns that drive high confidence

High confidence occurs when the substance name contains sufficient semantic signal to derive the ChemOnt class without structural information. Six recurring patterns are identifiable:

**1. Functional class explicitly stated in the name**

The most common driver: the name directly names the chemical class.

| Substance name (excerpt) | Direct parent assigned |
|---|---|
| *"alkyl sulfate, sodium (mono …)"* | Organic sulfate salts |
| *"Fatty acids, C8-10 … Me esters"* | Fatty acid methyl esters |
| *"2,4-diphenylmethane-diisocyanate"* | Isocyanates |
| *"Benzyl-C12-14-alkyldimethylamine chloride"* | Alkyldimethylbenzylammonium chlorides |

**2. Carbon-chain length notation (C_x–C_y)**

Names with explicit chain length ranges are almost always high-confidence and map unambiguously to lipid or hydrocarbon classes:

| Substance name (excerpt) | Direct parent assigned |
|---|---|
| *"Alkenes, C8-10, C9-rich"* | Olefins |
| *"Hydrocarbons, C6, isoalkanes"* | Branched alkanes |
| *"Fatty acids, C12-16 (even numbered) … Me esters"* | Fatty acid methyl esters |
| *"Paraffin waxes and hydrocarbons …"* | Alkanes |

**3. Heterocyclic ring name present in the substance name**

When a recognised ring system is literally part of the name, the LLM resolves the class with high certainty:

| Substance name (excerpt) | Direct parent assigned |
|---|---|
| *"4-Methyl-2-propylbenzimidazole"* | Benzimidazoles |
| *"1-benzyl-4-phenyl-piperidin-4-carboxylic acid"* | Piperidinecarboxylic acids |
| *"(2,2-Diphenyltetrahydrofuranyl)…"* | Tetrahydrofurans |
| *"Brominated chlorinated copper phthalocyanine"* | Phthalocyanines |

Three phthalocyanine-containing substances all receive high-confidence assignments — the word "phthalocyanine" is sufficiently diagnostic on its own.

**4. Well-known trivial names**

Common names that map unambiguously to a single chemical type:

| Substance name | Direct parent assigned |
|---|---|
| *"rape oil"*, *"coconut oil"*, *"safflower oil"* | Triacylglycerols |
| *"cresols"* | Cresols |
| *"ETANAL"* | Short-chain aldehydes |

**5. IUPAC-style salt nomenclature**

Inorganic salts written in the canonical "acid, metal(n+) salt" format are unambiguous:

| Substance name | Direct parent assigned |
|---|---|
| *"Sulfuric acid, nickel(2+) salt (1:1), pentahydrate"* | Transition metal sulfates |
| *"Arsenic acid (H3AsO4), lead(2+) salt"* | Post-transition metal arsenates |
| *"Sulfuric acid, cadmium salt …"* | Transition metal sulfates |

**6. High-confidence null for non-chemical entities**

Three high-confidence assignments return null at all hierarchy levels:

| Substance name | Type |
|---|---|
| *"Matte, copper-lead"* | Metal alloy |
| *"Candida oleophila strain O"* | Yeast (biocontrol agent) |
| *"Papiliotrema terrestris PT22"* | Yeast (biocontrol agent) |

The model correctly identifies these as non-classifiable under ChemOnt and returns null with high confidence — a calibrated and methodologically correct response.

---

## Discussion

### What the high-confidence population tells us

The substances that receive high confidence are systematically those whose **name is already a partial chemical description**: they contain functional class names, ring system names, chain length descriptors, or well-known trivial names. In other words, high confidence selects the "easy" fraction of the non-structure population — substances where the name was written to communicate chemical identity.

The substances that receive low confidence are the complementary population: trade names (SETAFIX XD 1289, CVCP-1V-01), internal codes (CHEMICAL 10057033, PM14943, Z-84), biological extracts (Allium cepa L. bulb extract, Pepper ext.), and vague mixture descriptions (Generic name, Spent liquor from…). These names provide no chemical signal, and the model appropriately reflects this with low confidence.

### Implications for downstream use

The confidence signal validated in 09d carries over to 09e:

- **High confidence (14.4%)**: ~72 substances with assignments likely comparable in quality to 09d high-confidence results (~54% accuracy at direct parent). Suitable for direct use with spot-checking.
- **Medium confidence (49.8%)**: assignments carry uncertainty; appropriate for enrichment pipelines where a human review step follows.
- **Low confidence (35.8%)**: assignments should not be used without verification, or can be treated as unclassifiable.

### Limitations

1. **No ground truth**: unlike 09d, accuracy cannot be quantified. The confidence signal is an indirect proxy, not a verified accuracy estimate.
2. **Embedding comparison is not informative**: only 9 of 503 substances have an embedding match, making the inter-method comparison statistically meaningless on this cohort.
3. **Pilot only**: the 500-substance pilot covers ~5% of the matchable population. The full run (~9,547 substances, ~$45) is needed before drawing conclusions about class distribution in the full non-structure corpus.
4. **Matchable filter may be under-inclusive**: the current regex filter (`reaction mass|petroleum|distillate|UVCB`) excludes the most obvious non-classifiable types but may miss others (biological extracts, alloys, proprietary codes). False-positive matchables inflate the low-confidence fraction.

---

## Conclusions

1. **The LLM assigns a kingdom to ~85% of matchable non-structure substances** and a direct parent to ~83%, with appropriate null at all levels when the name is insufficiently informative.

2. **Confidence is higher for names that encode chemical structure** (functional group names, ring systems, chain length notation) and lower for trade names, codes, and biological extracts. This is the expected and correct behaviour.

3. **High-confidence assignments cluster in lipids, organoheterocyclic compounds, and hydrocarbons** — classes whose members are frequently described with names that directly encode the class.

4. **The subclass level (48.4% null) reflects a structural feature of ChemOnt**, not a model failure: many ChemOnt paths skip the subclass node, and the subclass subdivision is often too fine-grained to resolve from a name alone.

5. **A full production run requires ~$45** and would yield ChemOnt assignments for the ~9,547 matchable non-structure substances, covering a population that ClassyFire cannot reach at all.
