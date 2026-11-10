# Reproducibility scope

This repository preserves the available data, code and outputs supporting the accompanying hydrogel study. The manuscript document and private editorial comments are not distributed here.

| Analysis | Supplied materials | Reproducibility scope |
|---|---|---|
| Dataset and cleaning | Source tables, retained records, exclusions, conflicts, corrections and reconstruction code | All 205 retained records reconstructed and verified |
| Historical benchmarks, Tables 2–3 | Unrounded aggregate report; recovered Ca/no-ion GP code, predictions and folds | Ca/no-ion GP can be refitted; six other subset-A families and all three original subset-B pipelines are unavailable |
| Source-aware GP and alternative models | Training grids, nested validation, per-record predictions and final models | Executable refitting and score reproduction |
| RBF-SVR and histogram gradient boosting | Audit code, fold manifest, predictions and scores | Executable additional comparator analyses |
| Four independent controls | Paired responses, exact doses and ratios | Recalculation available; four pairs do not establish an experimental error floor |
| Candidate search and preparation Monte Carlo | Search engine, candidate lists, stock recipes, seeds and perturbation outputs | Executable conditional simulation; candidates remain experimentally untested |
| Exploratory MCMC | Sampler, all 17 posterior archives, diagnostics and calibration | Archived-chain reproduction and separate fresh-sampling mode |
| Rheological target extraction | Imported G values | Original rheometer sweeps, Maxwell-fitting code/settings and fitting uncertainty unavailable |

## Distinct validation questions

The historical subset-A Ca/no-ion GP uses five-fold composition-grouped validation. Its kernel parameters are fitted by log marginal likelihood; an inner four-fold loop calibrates intervals. Its original pre-crosslinker polymer coordinates are preserved. This differs from the new-series source-aware GP, whose kernel settings are selected by an inner grouped grid search.

The new-series outer validation leaves out polymer compositions and removes matching historical compositions from training. Other new-series compositions remain in training. Headline metrics cover 30 positive records in 15 groups, with four independent repeat/near-repeat controls reported separately. These scores do not describe an untouched prospective test.

The historical-only transfer analysis is retrospective and separate from both the original historical benchmarks and new-series composition-blocked validation. Its RMSE(log10 G) is 0.3612048702, with median error factor 1.7392353 and 19/30 predictions within a factor of two.

The extended comparison includes RBF-SVR and scikit-learn histogram gradient boosting. Their RMSE(log10 G) values are 0.3093206632 and 0.2788343621, respectively. The retained point GP does not have the lowest RMSE. Histogram gradient boosting is distinct from the historical XGBoost implementation. All model-family comparisons remain exploratory after inspection of the same validation results.

## Interpretation of the data and simulations

The historical positives include 12 Fe/Al records. These are excluded from the historical divalent subset B but included in the full 198-positive GP fit. Fe and Al therefore occur in both data sources.

Nominal concentration coordinates use liquid-volume fractions and additive 1:1 polymer-liquid/crosslinker volumes. Experimental spraying or external-bath crosslinking does not establish a homogeneous twofold dilution inside the gel. The coordinates are model conventions, not direct measurements of gel composition or free-ion concentration. Ion valence does not change the volume-dilution factor.

Preparation Monte Carlo propagates specified input perturbations through a fixed GP. It omits model error, measurement error and gelation probability. At the 15 kPa target, 81.38% of simulated preparations leave the measured input support. Exploratory MCMC 80% intervals cover 18/30 headline observations; satisfactory chain diagnostics do not establish predictive calibration.

## Figure and source provenance

[Manuscript-source images](../figures/manuscript_as_supplied/) are preserved separately from [regenerated figures](../revision_figures/). Figure 7 uses a revised presentation of conditional spread. Regenerated figures are not claimed to be pixel-identical to the source images. The [manuscript map](manuscript_map.csv) identifies the corresponding inputs and scripts.

[Source fingerprints](source_provenance.json) retain exact input filenames and SHA-256 hashes for traceability. [Original release notes](original_v6_notes.txt) are an archival record; current commands and interpretation are given in the [README](../README.md) and [methods](methods.md).

## Materials required for broader reproduction

Complete reproduction of the remaining historical model families requires their original implementations, software versions, preprocessing and tuning procedures, exact fold assignments, and record-level predictions. Aggregate scores cannot recover these materials. Upstream reproduction of G additionally requires the original rheometer exports and Maxwell-fitting implementation/settings. Complete preparation/batch metadata and prospective candidate measurements are also absent from this deposit.
