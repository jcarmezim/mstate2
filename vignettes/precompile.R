## Precompute the vignettes. They use the DIVINE cohort, which is not
## distributed with mstate2, so the .Rmd files are knitted here from the
## .Rmd.orig sources and shipped with their output.
## From the package root:
##   MSTATE2_DIVINE_DATA=/path/to/MSM_Data.RData Rscript vignettes/precompile.R
local({
  old <- setwd("vignettes"); on.exit(setwd(old))
  unlink("figures", recursive = TRUE)
  for (f in c("mstate2", "reference"))
    knitr::knit(paste0(f, ".Rmd.orig"), paste0(f, ".Rmd"), quiet = TRUE)
})
