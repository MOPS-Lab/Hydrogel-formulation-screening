# Hydrogel formulation screening

Data, models and computational workflows accompanying **“From target rheology to formulation: benchmarking data-driven models for shear modulus prediction and computational candidate design for multi-component biopolymer hydrogels”**, by Aleš Ručigaj, Petar Vračar and Tilen Kopač.

This repository contains the curated formulation dataset, historical benchmark summaries and recovered Ca/no-ion Gaussian process (GP), composition-blocked new-series validation, preparation Monte Carlo, exploratory hyperparameter MCMC, and an offline formulation-screening application.

**Scope:** the response is an imported generalized-Maxwell-derived modulus G, conditional on gel formation. Candidate formulations have not been tested experimentally. Preparation Monte Carlo ranges describe specified input perturbations through a fixed model; they are not calibrated prediction intervals or synthesis-success probabilities.

## Start here

The main data and results can be inspected without installing software.

| To inspect | Open |
|---|---|
| Experimental records and units | [205 retained records](data/curated/records.csv) · [Data dictionary](data/data_dictionary.csv) |
| Predictive performance | [Model comparison](results/tables/table_4_extended_new_series_metrics.csv) · [Held-out predictions](revision_analysis/new_series_predictions.csv) |
| Candidate formulations | [Recipes and Monte Carlo summaries](results/tables/table_6_candidates_full_precision.csv) |
| A manuscript table or figure | [Reviewer guide](docs/REVIEWER_GUIDE.md) |
| Computational definitions and limitations | [Methods](docs/methods.md) · [Reproducibility scope](docs/reproducibility_status.md) |

CSV files can be viewed on GitHub or downloaded for spreadsheet software. The [reviewer guide](docs/REVIEWER_GUIDE.md) maps the manuscript results to their data and scripts.

## Use the offline application

Select **Code → Download ZIP**, extract the archive, and open `index.html` in a browser. The application runs offline and does not require Python or Node.js. Opening `index.html` in GitHub displays its source-file page.

The application supports forward prediction, preparation-scenario simulation and computational candidate search. Stock inputs and nominal final concentrations are distinguished in the interface. A support flag describes the location of a formulation relative to measured inputs; it does not certify prediction accuracy.

## Reproduce the analyses

The workflows were tested with Python 3.12.14 and Node.js 24.22.1. Python dependencies are pinned in [requirements.txt](requirements.txt); JavaScript uses Node built-ins. Node is required for command-line candidate search and application checks.

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python reproduce.py --mode check
```

`check` reconstructs all 205 input records, verifies group exclusions and saved new-series scores, refits the final GP, checks known non-gels and reproduces candidate preparation Monte Carlo. It also checks application logic when Node is installed. Outputs are written to `results/` and `revision_analysis/reconstructed_data/`.

```bash
python reproduce.py --mode all
```

`all` reruns the recovered historical Ca/no-ion GP, nested model comparisons, the GP/SVR/histogram-boosting audit, retrospective historical-only transfer, archived-chain diagnostics and calibration, candidate searches, scenario sensitivity, figure generation and checks. It overwrites generated outputs and sets numerical-library threads to one. Runtime depends on the machine. The unavailable historical pipelines and raw rheology fitting cannot be rerun from this deposit.

The standard workflow reuses the archived MCMC chains. Fresh sampling is available separately:

```bash
python reproduce.py --mode mcmc-fresh
```

This mode writes to `results/mcmc_fresh/` without replacing the archived release chains. It uses the original seeds, warm-up and convergence-extension rules. Rename an existing `results/mcmc_fresh/` directory before another uncached run.

Examples of individual commands:

```bash
python analysis_source/model/hydrogel.py --c 0.5525 0.2975 0 --ion Ca --dose 0.033785
python analysis_source/monte_carlo.py --stocks 1.7 1.7 1.7 --fractions 0.65 0.35 0 --ion Ca --dose 0.06757 --target 15000 --n 5000
python reproduce.py --mode figures
python build_app.py
```

The forward CLI uses **nominal final concentrations**; the Monte Carlo CLI uses **stock concentrations and stock crosslinker dose**. Under the specified 1:1 volume convention, nominal final concentrations are halved for both divalent and trivalent ions. This is a volume convention, not an ion-binding stoichiometry.

## Dataset and interpretation

- Historical series: 168 records, including 164 positive responses and four observed non-gels.
- New series (`KS`): 37 records, including 34 positive responses and three observed non-gels.
- Final GP: fitted to 198 positive responses; seven known failures are retained as exact-match exclusions.
- Historical subsets: A contains 116 Ca/no-ion records; B contains 152 divalent/no-ion records. Both have 24 polymer-composition groups. All historical non-calcium divalent records share one polymer composition.
- New-series validation: 16 outer composition folds cover 34 positives. Headline scores use 30 records in 15 groups; four independent repeat/near-repeat controls are reported separately. Other new-series compositions enter training, so this is composition-blocked evaluation, not an untouched prospective test.

The retained GP has RMSE(log₁₀ G) **0.302701**, median error factor **1.47244**, and **22/30** predictions within a factor of two. It does not have the lowest RMSE: the training-mean baseline scores 0.294181, RBF-SVR 0.309321, and histogram gradient boosting 0.278834. The [complete comparison](results/tables/table_4_extended_new_series_metrics.csv) includes the remaining models and R² values. Histogram gradient boosting is the scikit-learn implementation, distinct from the historical XGBoost benchmark.

Source response precision is preserved; the number of stored digits does not establish measurement accuracy. Nominal concentration coordinates assume additive 1:1 volumes and do not establish measured concentrations inside externally crosslinked gels.

## Repository map

| Path | Contents |
|---|---|
| `index.html`, `app/`, `build_app.py` | Offline application, template and build script |
| `data/raw/`, `data/curated/` | Source tables, retained records, cleaning decisions and provenance |
| `analysis_source/` | Deployed GP, Python/JavaScript engines, Monte Carlo and MCMC archives |
| `analysis_historical_gp/` | Recovered original Ca/no-ion GP and record-level predictions/folds |
| `analysis_v4/`, `analysis_v5/` | Earlier GP comparison and nested alternative-model comparison |
| `revision_analysis/` | Comparator audit, fold manifest, calibration, searches and sensitivity |
| `results/tables/` | Machine-readable tables and outer-fold membership |
| `figures/manuscript_as_supplied/`, `revision_figures/` | Manuscript-source images and separately regenerated figures |
| `docs/` | Methods, reviewer guide, manuscript mapping and reproducibility scope |
| `tests/`, `verify_release.py` | Numerical references and verification |

## Reproducibility boundaries

The historical Ca/no-ion GP and the new-series workflows are executable. The other six original subset-A families and all three original subset-B pipelines were not recovered; their aggregate reported scores are supplied. The retrospective historical-only transfer analysis is a separate analysis and does not replace those benchmarks.

Raw rheometer sweeps, the original Maxwell-fitting implementation/settings and parameter uncertainty, complete preparation/batch metadata, and prospective tests of the candidate formulations are not supplied. See the [analysis-by-analysis scope](docs/reproducibility_status.md).

## Citation and reuse

Use the metadata in [CITATION.cff](CITATION.cff) when citing this computational companion. The repository is available at <https://github.com/MOPS-Lab/Hydrogel-formulation-screening>. A publication DOI or archival DOI is not recorded in this version. For reproducible attribution, identify the tag or commit actually used.

Reuse terms are described in [LICENSING.md](LICENSING.md). [MANIFEST.sha256](MANIFEST.sha256) records checksums of the packaged files, excluding the manifest itself. The [data availability statement](docs/data_availability.txt) describes the supplied materials and their limits.
