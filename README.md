# Hydrogel formulation screening

Data and computational tools accompanying **“From target rheology to formulation: benchmarking data-driven models for shear modulus prediction and computational candidate design for multi-component biopolymer hydrogels”**, by Aleš Ručigaj, Petar Vračar and Tilen Kopač.

This release matches the available analyses behind `Manuscript_v02(1).docx`, including the later SVR and histogram-gradient-boosting comparison. It packages the source-aware Gaussian process retained in v6, composition-blocked validation, preparation Monte Carlo, exploratory hyperparameter MCMC and candidate search. No new experimental observations were generated.

**Scientific scope:** predictions concern the imported generalized-Maxwell-derived modulus G, conditional on gel formation. Candidate formulations are untested. Monte Carlo ranges describe specified input perturbations through a fixed model, not calibrated prediction intervals or synthesis-success probabilities.

## For reviewers: start here

**The key data and results can be reviewed without installing software.** Start with the five links below; the [reviewer guide](docs/REVIEWER_GUIDE.md) provides a clickable map to each manuscript table and figure.

| What would you like to inspect? | Open this |
|---|---|
| The experimental dataset | [All 205 retained records](data/curated/records.csv) · [Column definitions and units](data/data_dictionary.csv) |
| Prediction accuracy | [New-series model comparison](results/tables/table_4_extended_new_series_metrics.csv) · [Individual held-out predictions](revision_analysis/new_series_predictions.csv) |
| Proposed formulations | [Full-precision candidate recipes and Monte Carlo summaries](results/tables/table_6_candidates_full_precision.csv) |
| A particular manuscript table or figure | [Reviewer guide and clickable manuscript map](docs/REVIEWER_GUIDE.md) |
| The interactive formulation tool | **Code → Download ZIP → extract → open `index.html` locally.** No Python installation required. |

CSV files can be viewed as tables on GitHub or downloaded for Excel. Clicking `index.html` in the GitHub file list opens its file page; download it or the repository to run the offline application. The full collection of scripts, JSON model files and NPZ chains supports detailed reproduction; it is not necessary to open each file during an initial review.

**Scope of reproduction:** the historical Ca/no-ion GP and the new-series workflows are executable. The remaining original historical benchmark pipelines and raw rheological fitting inputs are not available. See [methods](docs/methods.md) and [reproducibility status](docs/reproducibility_status.md).

## Reproduce the available analyses

Tested with Python 3.12.14 and Node.js 24.22.1. Python dependencies are pinned in `requirements.txt`; JavaScript uses only Node built-ins. Node is needed for command-line candidate search, but not for opening the offline browser app.

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python reproduce.py --mode check
```

`check` reconstructs and verifies all 205 input records, checks training/test group exclusions, verifies reported new-series metrics, refits the final GP, checks known non-gels and reproduces the candidate preparation Monte Carlo. It also exercises the standalone application when Node is installed. Outputs are written to `results/` and `revision_analysis/reconstructed_data/`.

```bash
python reproduce.py --mode all
```

`all` reruns the historical Ca/no-ion GP, the v5 nested comparison, GP/SVR/histogram-boosting audit, retrospective historical-only test, archived-chain MCMC diagnostics and calibration, full candidate search, scenario sensitivity, figure generation and checks. It overwrites generated analysis outputs. The entry point sets numerical-library threads to one for reproducibility. Allow several minutes; runtime depends on the machine. This command **does not** rerun the unavailable remaining historical model families or fit rheological sweeps.

The standard reproduction uses the supplied MCMC chains. To run new chains rather than load them:

```bash
python reproduce.py --mode mcmc-fresh
```

This creates a separate `results/mcmc_fresh/` directory with no chain cache. Sampling uses the original seeds, warm-up and convergence-extension rules. Existing release models and chains are retained. Rename that output directory before repeating an uncached run.

Individual commands:

```bash
python analysis_source/model/hydrogel.py --c 0.5525 0.2975 0 --ion Ca --dose 0.033785
python analysis_source/monte_carlo.py --stocks 1.7 1.7 1.7 --fractions 0.65 0.35 0 --ion Ca --dose 0.06757 --target 15000 --n 5000
python reproduce.py --mode figures
python build_app.py
```

The forward CLI above uses **nominal final concentrations**. The Monte Carlo CLI uses **stock concentrations and stock crosslinker dose**. Their factor of two is a volume convention, independent of ion valence.

## Repository contents

| Path | Purpose |
|---|---|
| `index.html`, `app/template.html`, `build_app.py` | Standalone English interface and its reproducible build |
| `data/raw/` | Recovered historical source CSV, new-series spreadsheet and historical benchmark report |
| `data/curated/` | All 205 records, provenance, cleaning decisions and quarantined conflicts |
| `analysis_source/model/` | Deployed v6 GP, metadata/support geometry, Python and JavaScript engines, exploratory MCMC model |
| `analysis_source/mcmc/` | Full-data and 16 outer-fold posterior chains and diagnostics |
| `analysis_v4/` | Earlier source-aware GP comparison and frozen historical model snapshot |
| `analysis_v5/` | Nested GP/ridge/ensemble comparison and saved predictions |
| `analysis_historical_gp/` | Recovered original Ca/no-ion GP code and archived record-level predictions/folds |
| `revision_analysis/` | New comparators, fold manifest, calibration, searches and scenario analyses |
| `results/tables/` | CSV exports corresponding to Tables 1–6, plus outer-fold membership |
| `revision_figures/` | Regenerated figures from available outputs; see manuscript map for differences |
| `docs/` | Methods, data limits, manuscript mapping and provenance |
| `tests/`, `verify_release.py` | Numerical references and release verification |

## Data and validation

- Historical: 168 retained records, of which 164 have positive G and four are observed non-gels.
- New series (`KS`): 37 records, of which 34 have positive G and three are observed non-gels.
- Final model: 198 positive responses; seven failures are preserved as exact-match exclusions.
- Historical subsets: A = 116 Ca/no-ion records; B = 152 divalent/no-ion records. Both contain 24 polymer-composition groups. All historical non-calcium divalent records occupy one polymer composition.
- New-series evaluation: 16 outer composition folds; headline metrics use 30 positive non-control records in 15 groups. Four independent near-repeat/repeat controls remain separate. Other new-series compositions enter training; this is not an untouched prospective test.
- Headline GP: RMSE(log₁₀ G) 0.302701, median error factor 1.47244, 22/30 within a factor of two. Mean baseline: 0.294181; added SVR: 0.309321; histogram boosting: 0.278834. Histogram boosting is **not XGBoost**.
- All original source response precision is retained; this does not imply equivalent measurement accuracy.

## What is not available

The historical report contains the aggregate scores for the seven-family/three-family benchmarks. The Ca/no-ion GP implementation and its record-level predictions/folds were recovered and can be rerun with `python reproduce_historical_gp.py`. The other six subset-A families and all three subset-B original implementations, tuning details and record-level split/prediction exports were not recovered. Tables 2–3 can be reconstructed from the report, but the full historical model comparison cannot be rerun exactly from this deposit. The historical-only transfer analysis is a separate retrospective analysis and must not be substituted for those scores.

Raw rheometer sweeps, the original Maxwell fitting implementation and its parameter uncertainty, complete preparation/batch metadata, and prospective tests of the proposed recipes are not supplied. The final-concentration convention describes nominal coordinates under additive 1:1 volumes; it is not a measured concentration in a sprayed or externally bathed gel.

## Citation and licensing status

`CITATION.cff` contains the manuscript authors and release title. Add the actual publication DOI, repository URL and archived-release DOI when assigned; no DOI or publication status has been invented. A release checksum manifest is supplied in `MANIFEST.sha256`.

No public reuse licence was specified in the supplied materials. The authors should set the intended code/data licences before public release. This package does not assign third-party rights or silently apply a new licence to the research data.
