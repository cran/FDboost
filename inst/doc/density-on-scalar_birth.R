## ----plot-data, echo=TRUE, warning=FALSE, message=FALSE, fig.cap='Densities of births in Germany per year and sex in $B^2(\\delta)$. Years are coded by different colors and line types, see Figure \\ref{fig:legend}.', fig.height=3.5, fig.width=11, fig.pos='!h'----
# load FDboost package
library(FDboost)
# load birth_densities
data("birthDistribution", package = "FDboost")

# function to plot a matrix or vector containing functions in B^2(delta) or L^2_0(delta);
# Is used for densities, effects, predictions (also clr transformed)
plot_function <- function(plot_matrix, ...) {
  funplot(1:12, plot_matrix, xlab = "month", xaxp = c(1, 12, 11), pch = 20, ...)
  abline( h = 0, col = "grey", lwd = 0.5)
}

# function to create two plots (by sex) from a matrix containing densities or predictions
# (also clr transformed) for males in first half of rows and females in second half
plot_birth_densities <- function(birth_matrix, ylim = range(birth_matrix), ...) {
  par(mfrow = c(1, 2))
  for (k in 1:2) {
    n_obs <- nrow(birth_matrix) / 2
    obs <- 1:n_obs + (k - 1) * n_obs
    plot_function(birth_matrix[obs, ], main = c("Male", "Female")[k], 
                   ylim = ylim, col = rainbow(n_obs, start = 0.5, end = 1), 
                   lty = c(1, 2, 4, 5), ...)
  }
}

# Plot densities
plot_birth_densities(birthDistribution$birth_densities, ylab = "densities")

## ----legend, echo=TRUE, fig.cap='Coding of the years.', fig.height=1.2, fig.width=8, fig.pos='!h'----
# legend (also for later plots)
year_col <- rainbow(70, start = 0.5, end = 1)
year_lty <- c(1, 2, 4, 5)
par(mar = c(0, 0, 0, 0) + 0.1)
plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1, ylim = 0:1)
legend("top", xpd = TRUE, legend = 1950:2019, lty = year_lty, ncol = 10, bty = "n",
       text.col = year_col, col = year_col, cex = 0.7)

## ----clr-data, echo=TRUE, fig.cap='Clr transformed densities in $L^2_0(\\delta)$. Years are coded by different colors and line types, see Figure \\ref{fig:legend}.', fig.height=3.5, fig.width=11, fig.pos='!h'----
# The function clr() can be used to compute the clr-densities; Our reference measure delta
# corresponds to equal integration weights w = 1 for all density values
birth_densities_clr_test <- t(apply(birthDistribution$birth_densities, 1, clr, w = 1))
# Compare with clr-densities contained in data set
sum(birth_densities_clr_test != birthDistribution$birth_densities_clr)
# Plot clr-densities
plot_birth_densities(birthDistribution$birth_densities_clr, ylab = "clr-densities")

## ----fit-model, echo = TRUE, message=FALSE------------------------------------
model <- FDboost(birth_densities_clr ~ 1 + bolsc(sex, df = 1) + 
                   bbsc(year, df = 1, differences = 1),
                 # use bbsc() in timeformula to ensure integrate-to-zero constraint
                 timeformula = ~bbsc(month, df = 4, 
                                     # December is followed by January of subsequent year
                                     cyclic = TRUE, 
                                     # knots = {1, ..., 12} with additional boundary knot
                                     # 0 (coinciding with 12) due to cyclic = TRUE
                                     knots = 1:11, boundary.knots = c(0, 12), 
                                     # degree = 1 with these knots yields identity matrix 
                                     # as design matrix
                                     degree = 1),
                 data = birthDistribution, offset = 0, 
                 control = boost_control(mstop = 1000))

## ----get-mstop, echo = TRUE, message=FALSE------------------------------------
# set.seed(1708)
# folds <- applyFolds(model, folds = cv(rep(1, model$ydim[1]), type = "bootstrap", B = 10))
# ms <- mstop(folds) # = 999
ms <- 999 
model <- model[ms] 

## ----plot_effects_clr, echo = TRUE, fig.cap='Estimated effects in $L_0^2(\\delta)$. Years are coded by different colors and line types, see Figure \\ref{fig:legend}.', fig.height=3, fig.width=11, fig.pos='!h'----
# Plotting 'model' yields the clr-transformed effects
par(mfrow = c(1, 3))
plot(model, n1 = 12, n2 = 12)

## ----get_and_plot_effects, echo = TRUE, fig.cap='Estimated effects in $B^2(\\delta)$. Years are coded by different colors and line types, see Figure \\ref{fig:legend}.', fig.height=3, fig.width=11, fig.pos='!h'----
# Get estimated clr transformed effects; we use predict(), which returns a matrix of the
# same dimension as the response (140 x 12), i.e., we have to extract the respective rows;
# Alternatively, one could use coef(), but has to specify n1 = 12, n2 = 12 to get the den-
# sities at 1, ..., 12, which also only yields the year effect on a grid of 12 years

# all rows contain intercept
intercept_clr <- predict(model, which = 1)[1, ] 
# first 70 rows contain effect for sex = male, second 70 rows for sex = female
sex_clr <- predict(model, which = 2)[c(1, 71), ] 
# first 70 rows contain effect for years from 1950 to 2019, second 70 rows are repetition
year_clr <- predict(model, which = 3)[1:70, ]

sex_col <- c("blue", "red")
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2) + 0.1)

# Retransform to Bayes Hilbert space using clr(..., inverse = TRUE); Our reference measure
# delta corresponds to equal integration weights w = 1 for all function values
intercept <- clr(intercept_clr, w = 1, inverse = TRUE)
sex <- t(apply(sex_clr, 1, clr, w = 1, inverse = TRUE))
year <- t(apply(year_clr, 1, clr, w = 1, inverse = TRUE))

# Plot retransformed effects
plot_function(intercept, main = "Intercept", ylab = expression(hat(beta)[0]), 
              id = rep(1, 12)) # id is passed to funplot since intercept is a vector
plot_function(sex, main = "Effect of sex", col = sex_col, 
              ylab = expression(hat(beta)["sex"]))
legend("topleft", legend = c("sex = male", "sex = female"), text.col = sex_col, bty = "n")
plot_function(year, main = "Effect of year", col = year_col, 
              ylab = expression(hat(g)("year")), lty = year_lty)

## ----get_and_plot_predictions, echo = TRUE, fig.cap='Predicted densities in $B^2(\\delta)$. Years are coded by different colors and line types, see Figure \\ref{fig:legend}.', fig.height=3.5, fig.width=11, fig.pos='!h'----
predictions_clr <- predict(model)
predictions <- t(apply(predictions_clr, 1, clr, inverse = TRUE))
plot_birth_densities(predictions, ylim = range(birthDistribution$birth_densities), 
                     ylab = "predictions")

