library(tidyverse)
library(scales)
# library(tndbias)
set.seed(1)

ve_voi <- seq(0,1,0.1)
ci_res <-
  get_plot_data(ve_voi = ve_voi,
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

ci_res2 <-
  get_plot_data(ve_voi = ve_voi,
                iters = 1e3,
                sampleN = 200000,
                rr_voi.given.cv.status.min = 5.0,
                rr_voi.given.cv.status.max = 11.0,
                coverage_voi = 0.55,
                coverage_cv = 0.7,
                ve_cv = 0.9,
                risk_cv.dx.given.no.vx = 0.05,
                risk_voi.dx.given.no.vx = 0.05,
                weight = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  )



blacken_first <- function(n) {
  c("#000000", scales::hue_pal()(n-1))
}

ggplot(ci_res %>%
         filter(weight %in% c(0,0.5,0.25,0.1)) %>%
         mutate(baseline = weight == 0) %>%
         mutate(weight = as.character(weight)) %>%
         mutate(weight2 = ifelse(baseline, NA, weight)) %>%
         mutate(`97.5%` = ifelse(baseline, NA, `97.5%`)) %>%
         mutate(`2.5%` = ifelse(baseline, NA, `2.5%`)),
       aes(x = ve_voi, y = ve_voi + bias)) +
  geom_ribbon(aes(ymin = ve_voi + `2.5%`,
                  ymax = ve_voi + `97.5%`,
                  fill = weight),
              alpha = 0.25, show.legend = FALSE) +
  geom_line(aes(colour = weight),
            size = 1.5, na.rm = TRUE) +
  geom_point(aes(colour = weight2),
             size = 4, show.legend = FALSE, na.rm = TRUE) +
  scale_x_continuous(name = "True VE",
                     labels = label_percent(),
                     breaks = seq(0,1,0.2)) +
  scale_y_continuous(name = "Observed VE (estimated)\nand 95% CI",
                     labels = label_percent(),
                     breaks = seq(-1,1,0.2)) +
  scale_fill_manual(palette = blacken_first,
                    na.value = NA) +
  scale_colour_manual(name = "Percentage of controls\nincluded with the\nconfounding vaccine\ndisease",
                      labels = function(x) label_percent()(as.numeric(x)),
                      palette = blacken_first,
                      na.value = NA) +
  theme(text = element_text(size = 30),
        legend.title = element_text(size = 15)) +
  guides(size = 'none')



x=ci_res %>% filter(weight %in% c(0,0.5,0.25,0.1)) %>%
  mutate(baseline = weight == 0) %>%
  mutate(weight = as.character(weight)) %>%
  mutate(weight2 = ifelse(baseline, NA, weight))

colour = weight

ggplot(x,
       aes(x = ve_voi, y = ve_voi + bias, colour = weight)) +
  geom_line(aes(size = baseline)) +
  geom_point(aes(colour = weight2), show.legend = FALSE) +
  scale_x_continuous(name = "True VE", labels = label_percent()) +
  scale_y_continuous(name = "Observed VE (estimated)", labels = label_percent()) +
  scale_size_manual(values = c(`TRUE` = 1.5, `FALSE` = 1.25)) +
  scale_colour_manual(name = "Percentage of controls\nincluded with the\nconfounding vaccine\ndisease",
                      labels = function(x) label_percent()(as.numeric(x)),
                      palette = blacken_first,
                      na.value = NA) +
  theme(text = element_text(size = 30),
        legend.title = element_text(size = 15),
        legend.box = "horizontal")





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





