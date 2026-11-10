# Reviewer guide

The main data and results can be inspected without installing software. CSV files display as tables on GitHub and can also be downloaded for Excel. JSON model files and NPZ posterior archives are mainly for computational reproduction.

## The five-minute route

1. **Dataset:** [all 205 retained records](../data/curated/records.csv), with a [column and unit dictionary](../data/data_dictionary.csv). There are 198 positive responses and seven recorded non-gels. G is an imported Maxwell-derived descriptor, not G′ measured at one frequency.
2. **Prediction accuracy:** [new-series model comparison](../results/tables/table_4_extended_new_series_metrics.csv) and [individual held-out predictions](../revision_analysis/new_series_predictions.csv). Headline scoring covers 30 records in 15 held-out composition groups; other new-series records enter training. GP RMSE(log10 G) is 0.302701, median error factor 1.47244, and 22/30 predictions are within a factor of two. The GP does not have the lowest RMSE. The expanded table includes RBF-SVR and scikit-learn histogram gradient boosting as additional comparators.
3. **Recipe candidates:** [full-precision recipes and Monte Carlo summaries](../results/tables/table_6_candidates_full_precision.csv). These are untested computational candidates. Stock and nominal final concentrations are separate columns. Scenario target fractions are not synthesis-success probabilities.
4. **Preparation sensitivity and uncertainty:** [scenario sensitivity results](../revision_analysis/mc_sensitivity.json), [MCMC calibration](../revision_analysis/mcmc_calibration.json) and [regenerated figures](../revision_figures/). At the 15 kPa target, 81.38% of preparation perturbations leave the measured support; nominal MCMC 80% intervals cover 18/30 observations.
5. **Interactive tool:** download the repository using **Code → Download ZIP**, extract the ZIP, then open `index.html` in a browser. The application works offline; Python and Node.js are unnecessary for this route. Opening `index.html` on GitHub displays the file page rather than launching the application.

## Links by manuscript item

| Manuscript item | First file to open | Further detail |
|---|---|---|
| Dataset / Section 2.4 | [records.csv](../data/curated/records.csv) | [Raw input tables](../data/raw/) and [cleaning audit](../data/curated/historical_cleaning_audit.csv) |
| Table 1 | [Historical ion counts](../results/tables/table_1_historical_B_ion_counts.csv) | [Original benchmark report](../data/raw/historical_benchmark_report.xlsx) |
| Table 2 | [Historical subset A scores](../results/tables/table_2_historical_reported_metrics.csv) | [Recovered Ca/no-ion GP predictions and folds](../analysis_historical_gp/archived_results/oof_calcium_blocked_GP.csv) |
| Table 3 | [Historical subset B scores](../results/tables/table_3_historical_reported_metrics.csv) | Aggregate report; original subset-B fitting pipelines are unavailable |
| Table 4 and additional comparators | [Extended new-series comparison](../results/tables/table_4_extended_new_series_metrics.csv) | [Per-record predictions](../revision_analysis/new_series_predictions.csv), [outer fold membership](../results/tables/outer_fold_membership.csv), [complete inner/outer fold manifest](../revision_analysis/fold_manifest.json) |
| Table 5 | [Independent repeat controls](../results/tables/table_5_independent_repeat_controls.csv) | Includes exact doses; Ba and Al are near-repeats |
| Table 6 | [Candidate recipes](../results/tables/table_6_candidates_full_precision.csv) | [Complete restricted search](../revision_analysis/search_results.json) |
| Figures 1–8 | [Manuscript-source images](../figures/manuscript_as_supplied/) | [Regenerated figures](../revision_figures/); see the [figure provenance](reproducibility_status.md) |
| Methods and limitations | [Computational methods](methods.md) | [Reproducibility status](reproducibility_status.md) |
| Verification performed for this release | [Verification summary](../results/verification.json) | [Historical GP reproduction](../results/historical_GP_reproduction.json), [audit summary](../revision_analysis/audit_summary.json) |

The [complete manuscript map](manuscript_map.csv) also lists the relevant scripts. Original manuscript images and regenerated later figures are kept separately and are not claimed to be pixel-identical.

## If you want to rerun the work

Follow the installation commands in the [README](../README.md). Then:

```bash
python reproduce.py --mode check
python reproduce.py --mode all
```

`check` verifies input reconstruction, fold exclusions, saved scores, the final GP refit and Monte Carlo consistency. `all` reruns the available analyses and regenerates outputs. Standard MCMC reproduction uses archived chains; `--mode mcmc-fresh` starts separate fresh sampling. The README explains each mode and its dependencies.

## Boundaries of the evidence

The code for the original historical Ca/no-ion GP was recovered. The other six subset-A model families and all three original subset-B implementations were not recovered; their aggregate scores are supplied. Raw rheological sweeps and the original Maxwell-fitting implementation are also unavailable. The release therefore supports the available computational analyses without claiming complete reproduction of every historical result or experimental validation of proposed recipes.
