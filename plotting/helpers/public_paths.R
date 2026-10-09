# Package-relative plotting defaults. Environment variables may override them.
ije_simulation_root <- function(package_root) {
  Sys.getenv("MVMR_PROJECT_ROOT", unset = file.path(package_root, "simulation"))
}
ije_simulation_output <- function(package_root) {
  Sys.getenv("MVMR_SIM_OUTPUT_ROOT",
             unset = file.path(ije_simulation_root(package_root), "simulation_output"))
}
