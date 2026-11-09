# Reproducibility and manuscript alignment

This is a source/data release prepared against `Manuscript_v02(1).docx`. The manuscript itself and its private editorial comments are not included.

| Analysis | What is provided | Status |
|---|---|---|
| Dataset and cleaning | Source tables; all retained rows; exclusions/conflicts; corrections; executable reconstruction | Reconstructed and verified against all 205 retained records |
| Historical A/B benchmarks, Tables 2–3 | Original unrounded aggregate report, inputs, CSV exports, chart-generation code; recovered Ca/no-ion GP implementation and per-record folds/predictions | Ca/no-ion GP is reproducible; the other six A families and all three original B pipelines remain missing |
| Source-aware GP and v5 alternatives, Table 4 | Training grids, nested validation, per-record predictions, final fitted models | Refit and score reproduction available |
| Added SVR and histogram boosting in abstract | Exact later audit code, fold manifest, per-record predictions, scores | Reproduction available; current manuscript body/Table 4 do not yet describe/include both added models |
| Four controls, Table 5 | Original pairs, doses and response ratios | Recalculation available; four pairs do not define an experimental error floor |
| Candidate recipes, Table 6 | Search engine, complete candidate lists, full-precision stock recipes, seeds, perturbation outputs | Reproduction available; recipes are computational candidates, not experiments |
| MCMC | All 17 posterior archives, sampler, diagnostics and calibration code | Cached-chain reproduction available; explicit fresh-sampling mode included |
| Rheological target extraction | Imported G values only | Raw sweeps, original Maxwell fitting code/settings and fitting uncertainty are not supplied |

## Concrete discrepancies in the current manuscript

0. **Section 2.5 generalises inner-CV hyperparameter selection to every benchmark model.** The recovered historical Ca/no-ion GP fits kernel parameters by log marginal likelihood on each training set; its inner grouped loop calibrates intervals. Describe this separately from the source-aware GP grid selection. Its recovered original polymer coordinates are pre-crosslinker values; uniform halving before standardisation would give the same standardised features.
1. **Section 2.4 says Fe/Al occur in the new series only.** The retained historical positives include 12 Fe/Al records. They are excluded from the historical divalent benchmark B but included in the full 198-positive GP fit. Correct this sentence.
2. **The abstract adds SVR and histogram boosting, while Table 4 and the corresponding methods omit them.** The recovered audit reproduces RMSE 0.3093206632 and 0.2788343621. Use `results/tables/table_4_extended_new_series_metrics.csv` and the declared grids in `audit_revision.py` to align the methods, table and discussion. Histogram boosting must not be called XGBoost.
3. **The body says historical-only transfer was not evaluated.** The later supporting-data package contains a retrospective historical-only test (RMSE 0.3612048702, median factor 1.7392353, 19/30 within ×2). It is reproduced here as an additional analysis. It is not an untouched prospective experiment and is not interchangeable with new-series LOCO or the historical benchmarks.
4. **Experimental mixing description versus numerical convention.** The manuscript describes weight ratios and spraying/bath crosslinking; model inputs use liquid-volume fractions and nominal additive-volume dilution. Retain the explicit conditional interpretation, and verify the original experimental meaning before calling these measured gel concentrations.
5. **Figures are not all the same revision.** `figures/manuscript_as_supplied/` contains exact embedded images from the current Word file. `revision_figures/` contains regenerable later analysis figures. Figure 7 in particular uses a revised display of conditional spread. These are labelled separately so that a regenerated figure is not represented as a pixel-identical source figure.
6. **Availability wording must match the deposit.** Do not claim that every historical benchmark fit is exactly reproducible while its original code/folds are missing. The supplied aggregates and full new-series workflow are available. An accurate suggested paragraph is in `docs/data_availability.txt`.

Original v6 notes are retained as provenance in `docs/original_v6_notes.txt`. Their older wording does not override the limitations documented here.

## Remaining author-supplied materials

For a complete historical-benchmark deposit, add the original Python/R scripts, package versions, hyperparameter grids, preprocessing/selection procedures, exact outer/inner record assignments, and out-of-fold predictions for the six unrecovered subset-A families and all three subset-B models. The recovered original Ca/no-ion GP is in `analysis_historical_gp/`. Its nested inner loop calibrates intervals; the GP kernel is fitted by marginal likelihood, unlike the source-aware GP's inner grid selection. For upstream target reproducibility, add rheometer exports and the Maxwell fitting implementation/settings. These gaps cannot be recovered from aggregate scores alone.

Before public release, set an author-agreed licence, assign a real repository URL and use a versioned release or commit for the manuscript citation. No repository was published and no licence was granted by assembling this package.
