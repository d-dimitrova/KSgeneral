#!/usr/bin/env Rscript

# Compare this checkout against the currently published CRAN build of KSgeneral.
# The script installs CRAN and local builds into isolated temporary libraries,
# evaluates representative R-level entry points in separate R sessions, and
# fails if numeric results differ beyond floating-point tolerance.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script <- sub(file_arg, "", args[startsWith(args, file_arg)][1], fixed = TRUE)
repo <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
cran_lib <- tempfile("ksg-cran-lib-")
local_lib <- tempfile("ksg-local-lib-")
out_dir <- tempfile("ksg-compare-")
dir.create(cran_lib)
dir.create(local_lib)
dir.create(out_dir)

message("Installing CRAN KSgeneral into ", cran_lib)
install.packages("KSgeneral", lib = cran_lib, repos = "https://cloud.r-project.org", dependencies = TRUE, quiet = TRUE)
message("Installing local KSgeneral into ", local_lib)
.libPaths(c(cran_lib, .libPaths()))
install.packages(repo, lib = local_lib, repos = NULL, type = "source", quiet = TRUE)

case_script <- file.path(out_dir, "cases.R")
writeLines(c(
  'suppressPackageStartupMessages(library(KSgeneral))',
  'options(digits = 17)',
  'set.seed(20260722)',
  'discrete_cdf <- stepfun(c(0, 1, 2, 3), c(0, 0.1, 0.4, 0.8, 1.0))',
  'Mixed_cdf_example <- function(x) {',
  '  ifelse(x < 0, 0, ifelse(x < 1, 0.5, ifelse(x == 1, 0.75, 1 - 0.25 * exp(-(x - 1)))))',
  '}',
  'cases <- list(',
  '  cont_ks_cdf = cont_ks_cdf(0.12, 40),',
  '  cont_ks_c_cdf = cont_ks_c_cdf(0.12, 40),',
  '  cont_ks_test = cont_ks_test(c(0.1, 0.4, 0.9, 1.4, 2.0), "pexp")$p.value,',
  '  disc_ks_c_cdf = disc_ks_c_cdf(0.2, 8, discrete_cdf),',
  '  disc_ks_test = disc_ks_test(c(0, 1, 1, 2, 3), discrete_cdf)$p.value,',
  '  mixed_ks_c_cdf = mixed_ks_c_cdf(0.25, 8, c(0, 1), Mixed_cdf_example),',
  '  mixed_ks_test = mixed_ks_test(c(0, 0.2, 1, 1.5, 2.0), c(0, 1), Mixed_cdf_example)$p.value,',
  '  KS2sample_two_sided = KS2sample(c(1, 1, 2, 4), c(1, 3, 3, 4), alternative = "two.sided")$p.value,',
  '  KS2sample_less = KS2sample(c(1, 1, 2, 4), c(1, 3, 3, 4), alternative = "less")$p.value,',
  '  KS2sample_greater = KS2sample(c(1, 1, 2, 4), c(1, 3, 3, 4), alternative = "greater")$p.value,',
  '  Kuiper2sample = Kuiper2sample(c(0.1, 0.1, 0.5, 0.9), c(0.2, 0.5, 0.5, 0.8))$p.value',
  ')',
  'saveRDS(cases, Sys.getenv("KSG_OUT"))'
), case_script)

run_cases <- function(lib, out) {
  r_libs <- paste(c(lib, cran_lib), collapse = .Platform$path.sep)
  cmd <- sprintf('KSG_OUT=%s R_LIBS=%s R_LIBS_USER=%s Rscript --vanilla %s',
                 shQuote(out), shQuote(r_libs), shQuote(r_libs), shQuote(case_script))
  status <- system(cmd)
  if (status != 0) stop("case execution failed for ", lib, call. = FALSE)
}

cran_out <- file.path(out_dir, "cran.rds")
local_out <- file.path(out_dir, "local.rds")
run_cases(cran_lib, cran_out)
run_cases(local_lib, local_out)
cran <- readRDS(cran_out)
local <- readRDS(local_out)

ok <- TRUE
for (name in names(cran)) {
  same <- isTRUE(all.equal(cran[[name]], local[[name]], tolerance = 1e-12, check.attributes = FALSE))
  cat(sprintf("%-22s CRAN=% .17g local=% .17g %s\n", name, cran[[name]], local[[name]], if (same) "OK" else "DIFF"))
  ok <- ok && same
}
if (!ok) stop("local KSgeneral results differ from CRAN", call. = FALSE)
cat("All compared R-level results match CRAN within tolerance 1e-12.\n")
