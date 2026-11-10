# Computational conventions

## Inputs and target

The target is the imported positive Maxwell-derived modulus G in Pa. The model fits log10(G/Pa), not directly measured G′ at a specified frequency. Zero responses remain labelled non-gels and are excluded from logarithmic regression. No missing target is imputed.

The historical source CSV stores ion amounts as mol per 100 mL of crosslinker solution. Multiplication by 10 gives stock mol/L. The equal-volume convention then divides stock ion dose by two. Polymer component coordinates are stock concentration (percentage points) × polymer-liquid fraction ÷ 2. This conversion applies equally to divalent and trivalent salts; valence does not change the dilution factor.

`logdose = log10(1 + nominal_final_ion_M / 0.005)` is equivalent to `log10(1 + stock_ion_M / 0.01)` under the same 1:1 convention. The constant 0.005 M is a numerical feature scale, not a chemical binding threshold or dissociation constant.

For the new spreadsheet the measured response is column **P, “tt”**. Column O and the `G_Pa_full_precision` provenance field are earlier model predictions, not the new measured response. KS_6 fractions 0.33/0.33/0.33 are normalised to thirds. KS_23 assigns the 0.43 fraction to TOCNF, following the confirmed correction preserved in v4. `prepare_data.py` checks these transformations against every retained record.

Historical cleaning removes 65 additive-containing rows, collapses 32 exact duplicates and quarantines 14 rows with conflicting responses at identical recorded inputs, leaving 168 rows from 279. The full decisions and conflicting records are retained. This policy should not be interpreted as proof that every identical-input record was experimentally duplicated; original preparation identifiers are unavailable.

Fixed metadata: pH 7; crosslinking 24 h at 5 °C; rheology 25 °C; TOCNF carboxyl content 1.5 mmol/g dry matter; alginate G fraction 0.70 and M fraction 0.30. The model cannot estimate effects of these fixed conditions. Spraying/external-bath preparation does not establish physical homogeneous twofold dilution; numerical coordinates retain the convention used in the supplied models.

## Recovered historical Ca/no-ion GP

The original subset-A GP implementation, archived predictions and fold assignments are in `analysis_historical_gp/`. `reproduce_historical_gp.py` reproduces RMSE(log10 G) 0.134773550357283. It uses five-fold composition-grouped validation; each GP fits kernel hyperparameters by log marginal likelihood, with two optimiser restarts. The nested four-fold loop calibrates empirical intervals, rather than selecting kernel parameters by a candidate grid. Original pre-crosslinker polymer concentrations are preserved in this recovered pipeline. Uniformly halving each concentration before per-feature standardisation would leave its standardised values unchanged. The six remaining subset-A algorithms and all original subset-B implementations are not recovered.

## Source-aware GP and validation

Continuous features are three nominal polymer concentrations and logdose. Standardisation and the log-target mean are estimated only on training rows in every split. The Matérn-5/2 kernel is multiplied by ion-category and source-category correlation factors. Kernel amplitude is fixed to 1; cross-ion correlation is fixed to 0.3. Noise is an additive **variance** in log-target units.

Outer splits leave out one new-series polymer composition, using coordinates rounded to eight decimals. Historical records matching any held-out new-series composition are also removed. The 16 outer folds cover 34 positive new records. Four controls with historical matches are excluded only from headline scoring, leaving 30 records/15 groups. The inner four-fold GroupKFold repeats the composition exclusion, with training-only scaling and tuning. `revision_analysis/fold_manifest.json` contains every training, inner-validation and outer-test record ID.

The GP grid is length ∈ {0.5, 1.5, 4}, noise variance ∈ {0.03, 0.3}, source correlation ∈ {0.3, 0.9}. Minimum inner mean squared log error selects parameters. The final refit uses all 198 positives and selects length 4, noise variance 0.3 and source correlation 0.9. These final settings were not fixed before validation and need not equal every outer-fold choice. The final deployed model retains this GP even though other approaches have smaller headline RMSE.

The alternative ridge/GP/ensemble grids and nested algorithm selector are declared in `analysis_v5/train_v5.py`. The later `audit_revision.py` adds RBF-SVR (12 candidates) and scikit-learn histogram gradient boosting (12 candidates), using the same outer and inner composition exclusions. SVR tunes C, gamma and epsilon. Histogram boosting tunes learning rate, leaf count and L2 regularisation; early stopping is disabled. This is not the historical XGBoost implementation. Algorithm-family comparisons remain exploratory after inspecting the same validation results.

Reported factor errors are 10 to the absolute log10 prediction error. Median factor is 10 to the median absolute log error; geometric mean factor is 10 to its mean. Group-weighted metrics and paired cluster-bootstrap comparisons are also provided. Bootstrap intervals condition on the existing out-of-fold predictions and omit refitting/tuning instability.

## Input support, search and preparation Monte Carlo

Support is a convex hull for each polymer system and ion, computed in the four input coordinates with SVD reduction to the observed affine subspace. Combined-source and new-series-only support are separate. Numerical tolerance is 1e-7. Being inside a hull does not establish accuracy or probability coverage.

Search uses equal active stock concentrations from 1 to 3 wt.% in steps of 0.1; component-liquid fractions in steps of 0.05; 15 logarithmic ion-dose points plus observed new-series doses, then local refinement. Ca ALG/TOCNF search examines 7,182 points; Ca across all systems examines 87,318. Candidate ranking includes nominal target tolerance and support, followed by scenario target fraction and 90th-percentile target error. It screens up to 24 candidates with 1,000 draws and evaluates up to eight with 5,000 draws. It does not guarantee a global optimum or cover every feasible preparation procedure.

Independent uniform scenario half-widths are ±0.1 percentage points for polymer stocks, ±2% for each dispensed liquid volume, and ±5% for crosslinker stock concentration. Perturbed concentrations are recomputed from dispensed volumes. Seed 20261014 and the same integer random generator are used in Python and JavaScript. Perturbations outside support are counted, not discarded. The target-hit fraction is conditional on the model and valid draws; Wilson numerical intervals only describe Monte Carlo sampling precision.

Table 6 recipes are untested. For the 15 kPa target the nominal selected point is about 13.650 kPa and 4,069/5,000 perturbations (81.38%) are outside the four-point ALG/TOCNF/Ca hull. A narrow conditional spread does not demonstrate a narrow experimental prediction interval. `mc_sensitivity.py` additionally changes scenario widths, seed, sample count and shared calibration assumptions.

## Exploratory MCMC

The retained point GP is unchanged. The separate posterior integrates log length, log noise variance and logit source correlation. Normal prior means are log(2), log(0.1), log(0.7/0.3), with standard deviations 1, 1.2, 1.5. Mean centring is empirical and uses each training set; amplitude 1 and ion correlation 0.3 remain fixed. There are four adaptive random-walk Metropolis chains, 1,000 warm-up iterations each, and initially 2,500 retained iterations each, extended when necessary. Adaptation occurs only in warm-up.

Validation uses 128 representative posterior draws; the browser uses 32. Predictions include observation-noise variance in the exploratory intervals. Archived chains exist for the full fit and 16 outer folds. The original rank-normalised split-Rhat and bulk-ESS implementation is included; no additional convergence statistic is implied. The audited maximum Rhat is approximately 1.00943 and minimum bulk ESS approximately 499.83. Nominal 80% intervals cover 18/30 headline records. These checks do not establish calibrated predictive uncertainty.

## Software and provenance

The numerical GP, Monte Carlo and MCMC engines are recovered from v6. The later comparison and sensitivity scripts are recovered from the 2026-11-01 supporting-data package. New release files organise the inputs, verify reconstruction, export tables, build the English interface and provide portable entry points. `docs/source_provenance.json` fingerprints supplied source archives. No new scientific score is substituted for an unrecovered historical result.
