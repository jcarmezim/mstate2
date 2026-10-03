## Reproduction of Table 2 and the Section 6.3 evolution intervals of
## Najera-Zuloaga, Besalu and Gomez Melis (2025), "Second-Order Markov
## Multistate Models: Nonparametric Estimation and Inference", using the
## real, published DIVINE cohort data:
##   https://github.com/bruigtp/DIVINE/tree/main/msm_data  (MSM_Data.RData)
##
## This script is intentionally NOT part of the automated testthat suite:
## MSM_Data.RData belongs to the DIVINE project, not mstate2, and is not
## bundled here. R CMD check runs every top-level *.R file under tests/, so
## this script checks for the file locally and exits quietly (no error, no
## network access) if it isn't found -- safe to keep in the package for CRAN.
##
## To run it: download MSM_Data.RData from the URL above and either place it
## in the current working directory / a msm_data/ subdirectory, or point
## MSTATE2_DIVINE_DATA at its path, e.g.
##   MSTATE2_DIVINE_DATA=/path/to/MSM_Data.RData Rscript tests/DIVINE_reproduction.R

local({
  candidates <- c(
    Sys.getenv("MSTATE2_DIVINE_DATA", unset = NA),
    "MSM_Data.RData",
    file.path("msm_data", "MSM_Data.RData")
  )
  candidates <- candidates[!is.na(candidates)]
  found <- candidates[file.exists(candidates)]

  if (!length(found)) {
    message("DIVINE MSM_Data.RData not found locally; skipping reproduction.\n",
            "Get it from https://github.com/bruigtp/DIVINE/tree/main/msm_data ",
            "and set MSTATE2_DIVINE_DATA=<path> to run this script.")
    return(invisible())
  }

  library(mstate2)

  ## --- load the real cohort: 2076 patients, one row each ---------------
  e <- new.env()
  load(found[1], envir = e)
  MSM <- e$MSM   # columns: id, inistat, t.nosp, t.sp, t.recov, t.nimv, t.mv, disch.s, death.s

  ## --- sojourn durations -> daily panel, exactly ?sojourn_to_panel's example
  segs <- c(NSP = "t.nosp", SP = "t.sp", NIMV = "t.nimv", IMV = "t.mv", Recov = "t.recov")
  absb <- c(Disch = "disch.s", Death = "death.s")
  panel <- sojourn_to_panel(MSM, "id", segs, absb)

  st  <- c("NSP", "SP", "Recov", "NIMV", "IMV", "Disch", "Death")
  d   <- prep2(panel, states = st, check.consecutive = FALSE)
  fit <- P2est(d)

  get_p <- function(h, l) fit$estimate$p[fit$estimate$h == h & fit$estimate$j == "SP" &
                                          fit$estimate$l == l]

  published <- c(p124 = 0.224, p224 = 0.025, p125 = 0.151, p225 = 0.018)
  mine <- c(p124 = get_p("NSP", "NIMV"), p224 = get_p("SP", "NIMV"),
           p125 = get_p("NSP", "IMV"),  p225 = get_p("SP", "IMV"))

  cat("Table 2 reproduction (mstate2 on real DIVINE data vs. the published paper):\n")
  print(data.frame(published = published, mstate2 = round(mine, 4)))

  ## RPE is exact given the reconstructed panel; published values are
  ## rounded to 3dp, so allow a small tolerance for the panel-reconstruction
  ## step (sojourn_to_panel assumes a fixed NSP -> SP -> NIMV/IMV -> Recov
  ## visiting order, since only per-state total durations are published, not
  ## timestamped transitions).
  stopifnot(all(abs(mine - published) < 0.01))
  cat("\nAll four RPE estimates match the published Table 2 within tolerance.\n")

  ## First-order Markov assumption: the paper concludes it is violated for
  ## both SP -> NIMV and SP -> IMV (non-overlapping 95% CIs across NSP vs. SP).
  e  <- fit$estimate
  ci <- function(hh, ll) unlist(e[e$h == hh & e$j == "SP" & e$l == ll, c("lower", "upper")])
  for (ll in c("NIMV", "IMV"))
    stopifnot(ci("NSP", ll)[1] > ci("SP", ll)[2])
  cat("\n95% CIs of NSP -> SP and SP -> SP do not overlap for NIMV and IMV.\n")

  ## Evolution intervals / "memory horizon" (paper's Section 6.3: at least
  ## ~5 days for SP -> NIMV, ~6 days for SP -> IMV).
  cmp_nimv <- compare2(fit, h = c("NSP", "SP"), j = "SP", l = "NIMV", nsteps = 9)
  cmp_imv  <- compare2(fit, h = c("NSP", "SP"), j = "SP", l = "IMV",  nsteps = 9)
  cat("\nEvolution intervals (SP -> NIMV):\n"); summary(cmp_nimv)
  cat("\nEvolution intervals (SP -> IMV):\n");  summary(cmp_imv)

  cat("\nDONE: reproduction of the real-data illustration in Najera-Zuloaga,",
      "Besalu and Gomez Melis (2025) succeeded.\n")
})
