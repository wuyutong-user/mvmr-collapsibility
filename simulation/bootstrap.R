# bootstrap.R
# Source all simulation modules in dependency order.

.sim_project_dir <- local({
  # If this file was sourced with chdir=TRUE, getwd() is the project directory.
  # Otherwise, users should set option(mvmr.sim.project_dir="...").
  getOption("mvmr.sim.project_dir", getwd())
})

source(file.path(.sim_project_dir, "00_config.R"))
source(file.path(.sim_project_dir, "01_role_registry.R"))
source(file.path(.sim_project_dir, "02_instrument_geometry.R"))
source(file.path(.sim_project_dir, "03_dgp_continuous.R"))
source(file.path(.sim_project_dir, "04_dgp_binary.R"))
source(file.path(.sim_project_dir, "05_gwas_summary.R"))
source(file.path(.sim_project_dir, "06_ivw_estimators.R"))
source(file.path(.sim_project_dir, "07_conditional_strength.R"))
source(file.path(.sim_project_dir, "08_theorem_calibration.R"))
source(file.path(.sim_project_dir, "09_binary_scale_calibration.R"))
source(file.path(.sim_project_dir, "10_run_six_dag_primary.R"))
source(file.path(.sim_project_dir, "11_run_sensitivities.R"))
source(file.path(.sim_project_dir, "12_summarise_results.R"))
source(file.path(.sim_project_dir, "13_validation_tests.R"))

source(file.path(.sim_project_dir, "14_egger_sensitivity.R"))
