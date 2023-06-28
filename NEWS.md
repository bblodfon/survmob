# survmob 0.1.2

* add a `NEWS.md` file to track changes to the package
* add `reshape_results()` to `MOBenchmark` class
* add `fit_blme_model_cmp()` to fit Bayesian LME models using benchmarking results
* extend `quiet` parameter scope in `MOBenchmark` to silence messages
* add `use_callr` parameter in `SurvLPS` to control encapsulation of `aorsf` learner

# survmob 0.1.1

* update `mlr3` to `>=0.15.0` (oob_error fixed)
* update `mlr3mbo` to `>=v0.2.0` (default result assigner is best average from the archive)
* changes in some tuning hps for RSFs
* Added `parallel_run()` for faster execution of eFS
* add try/catch for `run()` method of eFS (`aorsf` issue)
* add `callr` encapsulation (guard against `aorsf` segmentation faults when #features ~ 30K)

# survmob 0.1.0

* First release version
