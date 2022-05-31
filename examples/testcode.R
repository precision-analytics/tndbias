library(ggplot)
options(mc.cores = 4) # Tune to your machine. Probably 4.
# Set the RNG so you can replicate parallel runs.
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

ve_voi <- seq(0,0.6,0.1)
debugonce(get_plot_data)
ci_res <-
  get_plot_data(ve_voi = 0.1,
            iters = 1e3,
            sampleN = 200000,
            rr_voi.given.cv.status = 8.0,
            coverage_voi = 0.55,
            coverage_cv = 0.7,
            ve_cv = 0.9,
            risk_cv.dx.given.no.vx = 0.05,
            risk_voi.dx.given.no.vx = 0.05,
            weight = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  )

ggplot(ci_res %>% filter(weight %in% c(0,0.5,0.25,0.1)),
       aes(x = ve_voi, y = ve_voi + bias, colour = weight)) +
  geom_line() + geom_point()





estimate_bias_parallel(
  iters = 1e3,
  sampleN = 200000,
  rr_voi.given.cv.status = 8.0,
  coverage_voi = 0.55,
  coverage_cv = 0.7,
  ve_cv = 0.9,
  ve_voi = 0.4,
  risk_cv.dx.given.no.vx = 0.05,
  risk_voi.dx.given.no.vx = 0.05,
  weight = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
)



sims <-
  estimate_bias_parallel(
    iters = 1e3,
    sampleN = 200000,
    rr_voi.given.cv.status = 8.0,
    coverage_voi = 0.55,
    coverage_cv = 0.7,
    ve_cv = 0.9,
    ve_voi = 0.4,
    risk_cv.dx.given.no.vx = 0.05,
    risk_voi.dx.given.no.vx = 0.05,
    weight = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  )

# True VE is an input on the bottom.
# So I need to replicated for five?









ci_res





ci_from_sims(sims)




# First, get my graph to be produced.


\deqn{p(x) =
  \frac{\lambda^x e^{-\lambda}}{x!}}


{p(x) = \lambda^x exp(-\lambda)/x!}

# Weight = % Controls with Confounding Dx
# sampleN = E6 = Total in population
# rr_voi.given.cv.status = F10 = D10/D11 =
# (Risk ratio of the vaccine of interest comparing confounding vaccine vs not confounding vaccine)





