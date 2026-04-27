# Experiment 09d — ChemOnt Classification via LLM

## Design

The experiment tests whether a large language model (Claude Sonnet 4.6) can classify chemical substances into the ChemOnt 2.1 taxonomy from substance names alone — without any structural information (InChIKey, SMILES).

### Baseline comparison

The design is identical to experiment 09c, which evaluates four embedding models on the same task: miniLM, mpnet, scibert, and biobert. The embedding models rank all ~4,825 ChemOnt classes by cosine similarity to the substance name; the LLM directly returns five nested classes.

### Ground truth

The ground truth consists of substances with a known structure (InChIKey) to which ClassyFire has assigned a ChemOnt class. The evaluation metric is `chemont_label`, the direct parent in the hierarchy (most specific level).

### LLM approach

The system prompt contains the full list of ChemOnt labels (kingdom → superclass → class → subclass / direct\_parent), cached after the first call. A single API call is made per substance, returning five levels as structured JSON:

```
kingdom → superclass → class → subclass → direct_parent
```

These five outputs form an implicit ranked list: `direct_parent` at rank 1, `kingdom` at rank 5. Hit@k and MRR are therefore directly comparable to the embedding results from 09c.

### Sample

A random sample of 500 was drawn from the full validation set (n ≈ 28,000 substance–class pairs) using `set.seed(42)`. The 95% confidence intervals apply to this sample and are calculated via the normal approximation for proportions:

```
SE = sqrt(p × (1 − p) / n)
CI = p ± 1.96 × SE
```

### Cost estimate

| Item | Cost |
|---|---|
| Cache write (8K tokens) | ~$0.02 |
| Cache reads (500 × 8K × $0.30/Mtok) | ~$1.20 |
| Output (500 × ~156 tokens × $15/Mtok) | ~$1.17 |
| **Total** | **~$2.40** |

Mean output was 156 tokens per call (median 154). The cache hit rate was 99.6%: after the first call, the ChemOnt labels were consistently read from the prompt cache.

---

## Results

### Hit@k and MRR

| Model | Hit@1 | 95% CI | Hit@3 | Hit@5 | MRR | n |
|---|---|---|---|---|---|---|
| miniLM | 5.75% | — | 9.86% | 12.34% | 9.33% | 28,204 |
| mpnet | 0.93% | — | 1.52% | 1.84% | 1.55% | 27,942 |
| scibert | 0.40% | — | 0.71% | 1.02% | 0.81% | 27,942 |
| biobert | 0.47% | — | 0.77% | 0.99% | 0.90% | 27,942 |
| **Claude Sonnet 4.6** | **39.2%** | [34.9–43.5%] | **41.2%** | **41.4%** | **40.2%** | 500 |

Claude outperforms the best embedding model (miniLM) by roughly seven-fold on Hit@1. The difference between Hit@1 and Hit@5 for Claude is small (39.2% vs. 41.4%), indicating that ranks 2–5 (higher levels in the hierarchy) yield very few additional hits: when the direct-parent level is wrong, the levels above it are typically wrong as well.

### Hit@1 per hierarchy level

| Level | Hit@1 | 95% CI |
|---|---|---|
| Kingdom | 96.5% | [94.9–98.0%] |
| Superclass | 52.3% | [47.9–56.6%] |
| Class | 40.7% | [36.4–45.1%] |
| Subclass | 27.6% | [23.7–31.4%] |
| Direct parent (= chemont\_label) | 39.2% | [34.9–43.5%] |

The pattern is monotonic: the broader the category, the better the performance. Kingdom (7 classes) is near-trivial; subclass (hundreds of classes) is the most difficult. Notably, *direct parent* scores slightly higher than subclass: the most specific class sometimes corresponds precisely to what the substance name suggests, even when the intermediate subclass level is already wrong.

### Performance per kingdom

| Kingdom | n | Hit@1 direct parent | Hit@1 kingdom |
|---|---|---|---|
| Organic compounds | 442 | 38.5% | 96.6% |
| Inorganic compounds | 40 | 52.5% | 97.5% |
| (no ground truth) | 18 | 27.8% | — |

Inorganic substances score *higher* at the direct-parent level than organic ones. This is likely because the inorganic ChemOnt hierarchy is shallower and less branched, making the step from kingdom to direct parent shorter.

---

## Pattern Analysis

### Confidence as a quality signal

The model reports a self-assessed confidence for each classification: `high`, `medium`, or `low`.

| Confidence | n | Hit@1 direct parent |
|---|---|---|
| high | 306 | 53.9% |
| medium | 172 | 18.0% |
| low | 22 | 0.0% |

The correlation is strong and practically useful: `high`-confidence calls are more than three times as accurate as `medium` ones, and `low` confidence does not yield a single correct direct parent. The confidence score is therefore a viable filter for downstream use.

### Errors at kingdom level

Kingdom errors are rare (n = 16 out of 482 with ground truth):

- 9 cases: the model returned `null` for kingdom — substances with ambiguous names (see below).
- 6 cases: an organic substance classified as inorganic — consistently for metal salts or metal chelates, where the organic component is overshadowed by the metal ion.
- 1 case: an inorganic substance incorrectly classified as organic.

### Direct-parent errors: structural patterns

The most common error classes for `direct_parent`:

1. **Hierarchical confusion**: the model selects a parent class rather than the correct child class. Examples: `Benzoic acid esters` → `p-Phthalate esters` (4×); `Fatty acid esters` → `Enoate esters` (3×). The model recognises the general type but misses the specific subtype.

2. **Salient functional group overshadows structural type**: for complex molecules with multiple functional groups, the model selects the most prominent group as the classification criterion, even when it is not the decisive one in the ChemOnt hierarchy. Examples:
   - An azo dye (ground truth: `2-naphthalene sulfonates`) is classified as `Azobenzenes` — the azo bond is visually dominant in the name, but the naphthalene sulphonate structure is what determines the ChemOnt position.
   - `1,1'-isopropylidenebis(p-phenyleneoxy)dipropan-2-ol` (ground truth: `Diphenylmethanes`) → `Phenol ethers`: the model identifies the phenoxy ether linkages but misses the diphenylmethane core.

3. **Subtype errors at subclass level**: the most frequent substitution is `[more specific class]` → `Carboxylic acid esters` or `Carboxylic acid derivatives` (repeatedly for benzoic acid derivatives, dicarboxylic acid derivatives, and acrylate derivatives). The model knows the supertype but not the finer subdivision.

4. **Null at subclass level** (18% null rate vs. ~2% at other levels): the model more often returns `null` at the subclass level than at higher levels. This reflects appropriate uncertainty: when the substance name provides insufficient information to determine the subclass, the model prefers `null` over a guess.

### JSON parse errors

Three substances produce the message `JSON parse error: Expecting value: line 1 column 1 (char 0)`:

- `Serricornin` — a pheromone with little name recognition
- `T001484` — an internal code with no chemical meaning
- `Tritosulfuron` — a sulphonylurea herbicide

In all three cases, the API likely returns an empty or invalid response. The cause is probably a name too complex or obscure for the model to produce a valid JSON output. These cases count as `null` at all levels and are missed in the evaluation.

---

## Discussion of Reasoning

The reasoning texts reveal how the model constructs its classification.

### Correct cases (high confidence)

For correctly classified substances, the reasoning follows a consistent pattern:

1. Identify the core structure or functional group from the name.
2. Link it to a ChemOnt class via an explicit chemical property.
3. Verify the hierarchical position.

Example (`Ethyl bromoacetate`, correct → `Alpha-halocarboxylic acid derivatives`):
> *"Ethyl bromoacetate (BrCH₂COOC₂H₅) is an ester of bromoacetic acid, where the bromine is on the alpha-carbon relative to the carboxyl group, making it an alpha-halocarboxylic acid derivative; it is also an organobromide."*

The reasoning is chemically correct and precise: the model identifies the alpha position of the bromine and links it to the appropriate ChemOnt class.

### Incorrect cases with explicit conflict

Particularly instructive are cases in which the model itself acknowledges a conflict but still makes the wrong choice:

Example (`Dipotassium [[N,N'-ethylenebis[N-(carboxylatomethyl)glycinato]](4-)…]zincate(2-)`, ground truth: `Tetracarboxylic acids and derivatives`, incorrect → `Amino acids, peptides, and analogues`):
> *"…the EDTA ligand contains amino acid-like carboxymethyl groups… It is best classified as an organic metal salt involving amino acid-derived chelating ligands, but the closest available ChemOnt label reflecting the organic acid/amino acid backbone of EDTA is 'Amino acids, peptides, and analogues'."*

The model recognises that the EDTA structure is in fact a tetracarboxylic acid, but nonetheless selects the amino acid category owing to structural similarity. This illustrates a fundamental limitation: the model optimises for structural analogy rather than for the exact ChemOnt definition.

### Medium confidence: focus problem

With medium-confidence errors, the model systematically selects a salient feature over the defining one:

Example (`imazapyr`, ground truth: `Alpha amino acids and derivatives`, incorrect → `Pyridinecarboxylic acids and derivatives`):
> *"Imazapyr is a herbicide containing a pyridine ring with a carboxylic acid group and an imidazolinone ring system; the pyridine carboxylic acid moiety is the most structurally defining feature in ChemOnt classification."*

The pyridine carboxylate group is indeed prominent in the name, but the alpha-amino acid structure (imidazolinone) is what determines the ChemOnt position. The model is unaware of ChemOnt's internal prioritisation rules.

---

## Conclusions

1. **LLM classification is viable** for an initial assignment of ChemOnt classes from substance names: ~39% accuracy at the most difficult level (direct parent), ~97% at kingdom.

2. **Confidence is a reliable quality signal**: high-confidence results (61% of cases) achieve 54% accuracy — sufficient for semi-automated enrichment with human verification of the remaining cases.

3. **Systematic error classes are identifiable** and partly avoidable: hierarchical confusion, salience bias for complex names, and null responses for ambiguous names. Fine-tuning on ChemOnt definition pairs or few-shot examples per class could reduce these errors.

4. **Unknown and code names** (T001484, pheromones, internal codes) cause parse failures. A more robust Python wrapper with fallback parsing and retry logic is advisable.

5. **The cost profile is favourable**: ~$2.40 for 500 substances, with 99.6% cache reuse. Scaling to the full validation set (~28,000 substances) would cost approximately $130 — a fraction of the cost of manual classification.
