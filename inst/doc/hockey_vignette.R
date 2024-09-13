## ----include = FALSE----------------------------------------------------------
# Save original options and par to reset at the end
original_options <- options()
original_par <- par(no.readonly=TRUE)

# save the built-in output hook
hook_output <- knitr::knit_hooks$get("output")

# set a new output hook to truncate text output
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(ns <- c(options$out.lines.head,
                       options$out.lines.tail))) {
    x <- xfun::split_lines(x)
    if (length(x) > sum(ns)) {
      # truncate the output
      x <- c(head(x, ns[1]), "...\n", tail(x, ns[2]))
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ">"
)

library(ggplot2)
library(devtools)
library(gridExtra)

## ----eval=FALSE---------------------------------------------------------------
#  # Install via github
#  devtools::install_github("apguthrie/BRcal")

## ----eval=FALSE---------------------------------------------------------------
#  # Install via tarball
#  install.packages("path_to_file\BRcal_1.0.0.tar.gz", repos = NULL, type="source")

## ----eval=FALSE---------------------------------------------------------------
#  # Install via CRAN
#  install.packages("BRcal")

## ----setup--------------------------------------------------------------------
library(BRcal)

## -----------------------------------------------------------------------------
data("hockey")

## -----------------------------------------------------------------------------
bt538 <- bayes_ms(hockey$x, hockey$y)
bt538

## -----------------------------------------------------------------------------
bayes_ms(hockey$x, hockey$y, Pmc=0.2)

## -----------------------------------------------------------------------------
llo_lrt(hockey$x, hockey$y)

## ----out.lines.head=4, out.lines.tail=23--------------------------------------
mle_recal(hockey$x, hockey$y)

## ----out.lines.head=3, out.lines.tail=3---------------------------------------
mle_recal(hockey$x, hockey$y, probs_only=TRUE)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(hockey$x, hockey$y)

## ----br90, cache=TRUE, out.lines.head=12, out.lines.tail=9--------------------
br90 <- brcal(hockey$x, hockey$y, t=0.9)

## ----out.lines.head=43, out.lines.tail=3--------------------------------------
br90

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(x=hockey$x, y=hockey$y, t_levels = c(0.95, 0.9))

## ----br90_tau, cache=TRUE, out.lines.head=12, out.lines.tail=9----------------
br90_tau <- brcal(hockey$x, hockey$y, t=0.9, tau=TRUE)

## ----out.lines.head=43, out.lines.tail=3--------------------------------------
br90_tau

## ----out.lines.head=12, out.lines.tail=14-------------------------------------
try(brcal(hockey$x, hockey$y, x0 = c(10, -5), start_at_MLEs = FALSE))

## ----out.lines.head=8, out.lines.tail=9---------------------------------------
br90_inMMA <- brcal(hockey$x, hockey$y, t=0.9, 
                    opts=list(local_opts=list(algorithm="NLOPT_LD_MMA")))

## -----------------------------------------------------------------------------
br90_inMMA$nloptr

## ----out.lines.head=8, out.lines.tail=9---------------------------------------
br90_outMMA <- brcal(hockey$x, hockey$y, t=0.9, x0=c(1, 2), start_at_MLEs = FALSE,
                     xtol_rel_outer = 1e-05, 
                     opts=list(algorithm="NLOPT_LD_MMA"))

## -----------------------------------------------------------------------------
br90_outMMA$nloptr

## ----out.lines.head=8, out.lines.tail=9---------------------------------------
br90_bound <- brcal(hockey$x, hockey$y, t=0.9, lb=c(0.25, -Inf), ub=c(5, Inf))

## -----------------------------------------------------------------------------
br90_bound$nloptr

## ----out.lines.head=8, out.lines.tail=9---------------------------------------
br90_stop <- brcal(hockey$x, hockey$y, t=0.9, maxeval=100, maxtime = 30, 
      xtol_rel_outer = 0.001, xtol_rel_inner = 0.001)

## -----------------------------------------------------------------------------
br90_stop$nloptr

## ----br95_2, cache=TRUE, out.lines.head=6, out.lines.tail=7-------------------
br90 <- brcal(hockey$x, hockey$y, t=0.9, print_level=1)

## ----br95_3, cache=TRUE, out.lines.head=9, out.lines.tail=10------------------
br90 <- brcal(hockey$x, hockey$y, t=0.9, print_level=2)

## ----br95, cache=TRUE, out.lines.head=12, out.lines.tail=9--------------------
br90 <- brcal(hockey$x, hockey$y, t=0.9, print_level=0)

## -----------------------------------------------------------------------------
bayes_ms(hockey$x, hockey$y, par = c(10, -5))

## -----------------------------------------------------------------------------
bayes_ms(hockey$x, hockey$y, control=list(reltol=1e-10))

## -----------------------------------------------------------------------------
bayes_ms(hockey$x, hockey$y, method = "L-BFGS-B", lower = c(0, -1), upper = c(10, 25))

## -----------------------------------------------------------------------------
bayes_ms(hockey$x, hockey$y, method = "L-BFGS-B", 
         lower = c(0, -Inf), upper = c(Inf, Inf))

## -----------------------------------------------------------------------------
bayes_ms(hockey$x, hockey$y, optim_details=FALSE)

## -----------------------------------------------------------------------------
br <- brcal(hockey$x, hockey$y, t=0.9, print_level = 0,
            optim_options=list(method="L-BFGS-B", 
                               lower=c(0.0001, 10), 
                               upper=c(0.0001, 10), 
                               control=list(maxit=200)))
br$nloptr

## ----out.lines.head=3, out.lines.tail=3---------------------------------------
hockey$x_mle <- LLO(x=hockey$x, delta=bt538$MLEs[1], gamma=bt538$MLEs[2])
hockey$x_mle

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
# plot original vs LLO-adjusted via MLEs
ggplot2::ggplot(data=hockey, mapping=aes(x=x, y=x_mle)) +
  stat_function(fun=LLO,
                args=list(delta=bt538$MLEs[1], 
                          gamma=bt538$MLEs[2]),
                geom="line",
                linewidth=1,
                color="black",
                xlim=c(0, 1)) +
  geom_point(aes(color=winner), alpha=0.75, size=2) +
  lims(x=c(0,1), y=c(0,1)) +
  labs(x="Original", y="LLO(x, 2, 2)", color="Winner") +
  theme_bw()

## ----out.lines.head=3, out.lines.tail=3---------------------------------------
LLO(x=hockey$x, br90$BR_params[1], br90$BR_params[2])

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
# LLO-adjust using delta=2, gamma=1000
hockey$x2 <- LLO(hockey$x, delta=2, gamma=1000)

# plot original vs LLO-adjusted via 2,2
ggplot2::ggplot(data=hockey, mapping=aes(x=x, y=x2)) +
  stat_function(fun=LLO,
                args=list(delta=2, gamma=1000),
                geom="line",
                linewidth=1,
                color="black",
                xlim=c(0, 1)) +
  geom_point(aes(color=winner), alpha=0.75, size=2) +
  lims(x=c(0,1), y=c(0,1)) +
  labs(x="Original", y="LLO(x, 2, 1000)", color="Winner") +
  theme_bw()

## -----------------------------------------------------------------------------
bayes_ms(hockey$x, hockey$winner, event="home", optim_details=FALSE)

## -----------------------------------------------------------------------------
bayes_ms(1-hockey$x, hockey$y, event=0, optim_details=FALSE)

## -----------------------------------------------------------------------------
bayes_ms(1-hockey$x, hockey$winner, event="away", optim_details=FALSE)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
plot_params(hockey$x, hockey$y)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
plot_params(hockey$x, hockey$y, dlim=c(0.5, 1.5), glim=c(0.25, 2.75))

## ----eval=FALSE, fig.align='center', fig.height=4, fig.width=6, include=TRUE----
#  plot_params(hockey$x, hockey$y, k=200, dlim=c(0.5, 1.5), glim=c(0.25, 2.75))

## ----echo=FALSE, fig.align='center', fig.height=4, fig.width=6----------------
zmat_list <- plot_params(hockey$x, hockey$y, k=200, dlim=c(0.5, 1.5), glim=c(0.25, 2.75), 
                         return_z = TRUE)

## ----save_z, eval=FALSE, fig.align='center', fig.height=4, fig.width=6, include=TRUE----
#  zmat_list <- plot_params(hockey$x, hockey$y, k=200, dlim=c(0.5, 1.5), glim=c(0.25, 2.75),
#                           return_z = TRUE)

## ----echo=FALSE, fig.align='center', fig.height=4, fig.width=6----------------
plot_params(z=zmat_list$z, k=200, dlim=c(0.5, 1.5), glim=c(0.25, 2.75))

## ----out.lines.head=4, out.lines.tail=9---------------------------------------
zmat_list

## -----------------------------------------------------------------------------
class(zmat_list$z)
dim(zmat_list$z)

## ----fig.width=8, fig.height=4, fig.align='center'----------------------------
par(mfrow=c(1,3))
plot_params(z=zmat_list$z, k=200, dlim=c(0.5, 1.5), glim=c(0.25, 2.75), main="Same range as z")
plot_params(z=zmat_list$z, k=200, dlim=c(0.7, 1.3), glim=c(0.5, 2.5), main="Smaller range than z")
plot_params(z=zmat_list$z, k=200, dlim=c(1e-04, 5), glim=c(1e-04, 5), main="Larger range than z")

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
plot_params(z=zmat_list$z, t_levels=c(0.95, 0.9, 0.8),
            k=200, dlim=c(0.5, 1.5), glim=c(0.25, 2.75))

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
plot_params(z=zmat_list$z, t_levels=c(0.99, 0.95, 0.9, 0.8, 0.7),
            k=200, dlim=c(0.5, 1.5), glim=c(0.25, 2.75), contours_only=TRUE)

## ----fig.width=6, fig.height=6, fig.align='center'----------------------------
plot_params(z=zmat_list$z, k=200, dlim=c(0.5, 1.5), glim=c(0.25, 2.75),
            imgplt_options=list(horizontal = TRUE, nlevel=10, legend.lab="Posterior Model Prob"))

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
plot_params(z=zmat_list$z, k=200, dlim=c(0.5, 1.5), glim=c(0.25, 2.75), t_levels = c(0.99, 0.1), 
            contour_options=list(lty = "dotted", col="hotpink", lwd=2))

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(x=hockey$x, y=hockey$y)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(x=hockey$x, y=hockey$y, t_levels = c(0.95, 0.9, 0.8))

## -----------------------------------------------------------------------------
lp <- lineplot(x=hockey$x, y=hockey$y, return_df = TRUE, t_levels = c(0.95, 0.9, 0.8))

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lp$plot

## ----out.lines.head=7, out.lines.tail=6---------------------------------------
lp$df

## -----------------------------------------------------------------------------
summary(lp$df)

## -----------------------------------------------------------------------------
nrow(lp$df)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(df=lp$df, t_levels=c(0.95, 0.9, 0.8))

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(df=lp$df, t_levels=c(0.95))

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(df=lp$df, t_levels=c(0.95, 0.9, 0.8), plot_original=FALSE)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(df=lp$df, t_levels=c(0.95, 0.9, 0.8), ggpoint_options=list(size=3, shape=1))

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(df=lp$df, t_levels=c(0.95, 0.9, 0.8), ggline_options=list(linewidth=2, alpha=0.1))

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(df=lp$df, t_levels=c(0.95, 0.9, 0.8)) + 
  scale_color_manual(values = c("orchid", "darkgreen"))

## ----fig.width=10, fig.height=4, fig.align='center'---------------------------
lp1 <- lineplot(df=lp$df, t_levels=c(0.95, 0.9, 0.8), thin_to=500)
lp2 <- lineplot(df=lp$df, t_levels=c(0.95, 0.9, 0.8), thin_prop=0.5)
lp3 <- lineplot(df=lp$df, t_levels=c(0.95, 0.9, 0.8), thin_by=2)
gridExtra::grid.arrange(lp1, lp2, lp3, ncol=3)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(df=lp$df, t_levels=c(0.95, 0.9, 0.8), thin_to=500, seed = 47)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
lineplot(df=lp$df, t_levels=c(0.95, 0.9, 0.8), thin_to=500, seed=NULL)

## ----include=FALSE------------------------------------------------------------
# Reset user options and par
options(original_options)
par(original_par)

