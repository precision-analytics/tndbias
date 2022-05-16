#' One iteration of the simulation.
#'
#' @inheritParams estimate_bias_parallel
#'
#' @export
one_iter <- function(sampleN,
                     rr_voi.given.cv.status = 8.0,
                     coverage_voi = 0.7,
                     coverage_cv = 0.55,
                     ve_cv = 0.4,
                     ve_voi = 0.9,
                     risk_cv.dx.given.no.vx = 0.05,
                     risk_voi.dx.given.no.vx = 0.05,
                     weight = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
                     fixed = FALSE) {

  # If the 2x2 coverage table is fixed, use it to determine the probability of
  # VOI vaccination.

  a <-
    (coverage_voi * coverage_cv * rr_voi.given.cv.status) /
    (1 - coverage_cv + (rr_voi.given.cv.status * coverage_cv))

  prob = c(a,                  coverage_cv - a,
           coverage_voi - a, 1+a-(coverage_cv+coverage_voi))

  if(!fixed){
    two.by.two <-
      table(factor(
        sample(1:4, size = sampleN,
               prob = prob,
               replace = TRUE), levels = 1:4))
  } else {
    two.by.two <- prob * sampleN
  }
  # Additionally, the probability of VOI vaccination given CV vaccination are
  # fixed given the three parameters (coverage_cv, coverage_voi, rr_voi.given.cv.status).

    risk_voi.vx.given.cv.vx <- prob[1] / (prob[1] + prob[2])
  risk_voi.vx.given.no.cv.vx <- prob[3] / (prob[3] + prob[4])

  two.by.two <- matrix(two.by.two, byrow = TRUE, nrow = 2)
  dimnames(two.by.two) <- list(CV = c("Vx", "Unvx"), VOI = c("Vx", "Unvx"))

  num_cv.vx <- sum(two.by.two["Vx",])
  num_no.cv.vx <- sum(two.by.two["Unvx",])

  num_voi.vx <- sum(two.by.two[,"Vx"])
  num_no.voi.vx <- sum(two.by.two[,"Unvx"])

  ### Risk of CV disease
  # Given the estimated probability of CV disease given no vaccine specified in the
  # argument in the function, use the CV VE to back calculate the probability
  # of CV disease given the CV vaccine.

   risk_cv.dx.given.cv.vx <-
    risk_cv.dx.given.no.vx * (1-ve_cv)

  if(!fixed){
    no_cv.dx.given.cv.vx <-
      rbinom(1, num_cv.vx, risk_cv.dx.given.cv.vx)
    no_cv.dx.given.no.cv.vx  <-
      rbinom(1, num_no.cv.vx, risk_cv.dx.given.no.vx)
    ctvx_cv <-
      rbinom(1, no_cv.dx.given.cv.vx, risk_voi.vx.given.cv.vx) +
      rbinom(1, no_cv.dx.given.no.cv.vx, risk_voi.vx.given.no.cv.vx)
    ctvx_cv_unbiased <-
      rbinom(1,
             no_cv.dx.given.cv.vx + no_cv.dx.given.no.cv.vx,
             coverage_voi)
  } else {
    no_cv.dx.given.cv.vx <-
      num_cv.vx * risk_cv.dx.given.cv.vx
    no_cv.dx.given.no.cv.vx  <-
      num_no.cv.vx * risk_cv.dx.given.no.vx
    ctvx_cv <-
      no_cv.dx.given.cv.vx * risk_voi.vx.given.cv.vx +
      no_cv.dx.given.no.cv.vx * risk_voi.vx.given.no.cv.vx
  }
  no_cv.dx <- no_cv.dx.given.cv.vx + no_cv.dx.given.no.cv.vx
  ctnovx_cv <- no_cv.dx - ctvx_cv

  ### Risk of VOI disease
  # Given the estimated probability of VOI disease given no VOI vaccine specified in
  # the argument in the function, use the VOI VE to back calculate the
  # probability of VOI disease given the VOI vaccine.

   risk_voi.dx.given.voi.vx <-
    risk_voi.dx.given.no.vx * (1-ve_voi)

  if(!fixed) {
    no_voi.dx.given.voi.vx <- rbinom(1, num_voi.vx, risk_voi.dx.given.voi.vx)
    no_voi.dx.given.no.voi.vx <- rbinom(1, num_no.voi.vx, risk_voi.dx.given.no.vx)
  } else {
    no_voi.dx.given.voi.vx <- num_voi.vx * risk_voi.dx.given.voi.vx
    no_voi.dx.given.no.voi.vx <- num_no.voi.vx * risk_voi.dx.given.no.vx
  }

  # Unbiased counts


  ### Calculate bias
  odds_voi <- no_voi.dx.given.voi.vx / no_voi.dx.given.no.voi.vx
  odds_cv <- ctvx_cv / ctnovx_cv
  if(!fixed) {
    unbiased_odds_cv <- ctvx_cv_unbiased / (no_cv.dx - ctvx_cv_unbiased)
  } else {
    unbiased_odds_cv <- coverage_voi / (1-coverage_voi)
  }

  VE <- 1 - (odds_voi / (odds_cv*weight + unbiased_odds_cv*(1-weight)))

  data.frame(calculated_VE = VE - ve_voi, bias_weight = weight)
}



#' estimate_bias_parallel
#'
#' @param iters Number of iterations simulation will run
#' @param sampleN The total number in the population
#' @param rr_voi.given.cv.status Risk ratio of vaccine of interest given confounding vaccination status.
#' \deqn{\frac{Pr(vaccine of interest | vaccinated with confounding vaccine)}{Pr(vaccine of interest | not vaccinated with confounding vaccine)}}{Pr(vaccine of interest | vaccinated with confounding vaccine) / Pr(vaccine of interest | not vaccinated with confounding vaccine}
#' @param coverage_voi Vaccine of interest coverage. Pr(vaccine of interest).
#' @param coverage_cv Confounding vaccine coverage. Pr(confounding vaccine).
#' @param ve_cv Effectiveness of confounding vaccine.
#' 1-(Pr(confounding vaccine disease | vaccinated with confounding vaccine) / Pr(confounding vaccine disease | not vaccinated with confounding vaccine))
#' @param ve_voi Effectiveness of vaccine of interest (VOI).
#' 1-(Pr(vaccine of interest disease | vaccinated with vaccine of interest) / Pr(vaccine of interest disease | not vaccinated with vaccine of interest))
#' @param risk_cv.dx.given.no.vx Risk of confounding vaccine disease among the unvaccinated. Pr(confounding vaccine disease | not vaccinated with confounding vaccine).
#' @param risk_voi.dx.given.no.vx Risk of vaccine of interest disease among the unvaccinated. Pr(vaccine of interest disease | not vaccineated with vaccine of interest).
#' @param weight Pr(confounding vaccine disease | control). A vector of numeric may be specified.
#'
#' @export
#'
#' @importFrom parallel mcmapply
estimate_bias_parallel <- function(iters,
                                   sampleN,
                                   rr_voi.given.cv.status = 1.5,
                                   coverage_voi = 0.7,
                                   coverage_cv = 0.55,
                                   ve_cv = 0.4,
                                   ve_voi = 0.9,
                                   risk_cv.dx.given.no.vx = 0.0215,
                                   risk_voi.dx.given.no.vx = 0.029,
                                   weight = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
) {
  sims <- mcmapply(
    function(x) one_iter(
      sampleN,
      rr_voi.given.cv.status = rr_voi.given.cv.status,
      coverage_voi = coverage_voi,
      coverage_cv = coverage_cv,
      ve_cv = ve_cv,
      ve_voi = ve_voi,
      risk_cv.dx.given.no.vx = risk_cv.dx.given.no.vx,
      risk_voi.dx.given.no.vx = risk_voi.dx.given.no.vx,
      weight = weight,
      fixed = FALSE
    )$calculated_VE,
    1:iters,
    mc.set.seed = TRUE
  )
  rownames(sims) <- weight
  sims
}


#' ci_from_sims
#'
#' @param sims sims
#'
#' @export
#'
ci_from_sims <- function(sims) {
  sim_quantiles <- data.frame(t(apply(sims, 1, quantile, c(0.5, 0.025, 0.975))),
                              check.names = FALSE)

  rnd <- function(x) sprintf("% .1f", round(x, 2))

  sim_quantiles$bias_ci <-
    with(sim_quantiles,
         paste0(rnd(`50%`*100), " (", rnd(`2.5%`*100), ", ", rnd(`97.5%`*100), ")"))

  sim_quantiles["bias_ci"]
}



######## Example simulation ###################

# inputs for example

# population size = 200,000
# iterations = 10,000

# RR VOI vaccination|CV vaccination status = 8.0

# VOI coverage = 55%
# VOI true VE = 40%
# VOI disease incidence (risk) among unvaccinated = 5%

# CV coverage = 70%
# CV true VE = 90%
# CV disease incidence (risk) among unvaccinated = 5%

#weight = proportion of controls with CV disease included c(0%, 10%, 25%, 50%, 75%, 100%)

# output returns mean bias in VE estimates (estimated VOI VE - true VOI VE) and 95% CI with different proportions of CV disease controls included



