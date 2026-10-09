## Precompute the vignettes. They use the DIVINE cohort, which is not
## distributed with mstate2, so the .Rmd files are knitted here from the
## .Rmd.orig sources and shipped with their output.
## From the package root:
##   MSTATE2_DIVINE_DATA=/path/to/MSM_Data.RData Rscript vignettes/precompile.R
local({
  old <- setwd("vignettes"); on.exit(setwd(old))
  unlink("figures", recursive = TRUE)
  for (f in c("mstate2", "reference", "paper"))
    knitr::knit(paste0(f, ".Rmd.orig"), paste0(f, ".Rmd"), quiet = TRUE)
})

## The home page of the pkgdown site is the main vignette (not the README):
## pkgdown/index.md = a short header + the body of vignettes/mstate2.Rmd, and
## its figures are copied to pkgdown/assets/figures so the paths still work.
local({
  body <- readLines("vignettes/mstate2.Rmd", encoding = "UTF-8")
  yaml_end <- which(body == "---")[2]
  body <- body[-seq_len(yaml_end)]
  body <- body[!grepl("^<!-- This vignette is precomputed|^     vignettes/precompile.R", body)]
  header <- c(
    "# mstate2 <img src=\"man/figures/logo.png\" align=\"right\" height=\"139\" alt=\"mstate2 logo\"/>",
    "",
    "**Second-order Markov multistate models in R**, illustrated on the DIVINE",
    "cohort. Install the development version from GitHub:",
    "",
    "```r",
    "# install.packages(\"remotes\")",
    "remotes::install_github(\"jcarmezim/mstate2\", build_vignettes = TRUE)",
    "```",
    ""
  )
  dir.create("pkgdown/assets/figures", recursive = TRUE, showWarnings = FALSE)
  unlink(list.files("pkgdown/assets/figures", full.names = TRUE))
  file.copy(list.files("vignettes/figures", pattern = "^guide-", full.names = TRUE),
            "pkgdown/assets/figures", overwrite = TRUE)
  body <- sub("^This vignette shows", "This guide shows", body)
  writeLines(c(header, body), "pkgdown/index.md", useBytes = TRUE)
})
