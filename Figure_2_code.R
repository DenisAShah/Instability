### Figure 2 code (it's not pretty and it's not quick)
#-- No, it's not, on both counts.  I made changes to (i) reduce the repetitiveness (follow DRY - do not repeat yourself), (ii) update broken rms code to work with rms ver. 8.2-1, (iii) remove the dependency on the reshape2 package, which is retired, (iv) use the native pipe operator (|>) instead of the magrittr pipe operator (%>%). 
# The massive nested for loop is a bottleneck (takes hours to run), and that is where future improvements could be made.

### Looking at model stability in terms of variability in predictions 
### as a function of sample size.
### Backwards Elimination models with and without shrinkage
### compared against elastic net

rm(list = ls())
library(MASS)
library(glmnet)
library(rms)
library(progress)
library(tidyverse)
library(dcurves)

# No longer needed:
# library(data.table)

# Not used in the script:
#library(doMC)
#library(ggh4x)
#library(ranger)
#library(doParallel)
# library(ragg)
#library(facetscales)
#library(apricom)
#library(patchwork)

# I'm working in Positron now:
# if(!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])

### define some functions
# function to generate the data
generate_data_1_predictor_for_paper <- function(
  NN,
  n.true.predictors = n.true.predictors,
  n.noise.predictors = n.noise.predictors,
  cor0 = 0,
  cor1 = 0,
  beta.0 = 0,
  var0 = 4
) {
  n.predictors <- n.true.predictors + n.noise.predictors
  mu0 <- rep(0, n.predictors)

  # Specify correlation matrix
  Sigma0 <- matrix(0, nrow = n.predictors, ncol = n.predictors)
  Sigma0[1:n.true.predictors, 1:n.true.predictors] <- cor0
  Sigma0[
    (n.true.predictors + 1):n.predictors,
    (n.true.predictors + 1):n.predictors
  ] <- cor1
  diag(Sigma0) <- c(var0, rep(1.0, n.noise.predictors))

  x <- mvrnorm(NN, mu0, Sigma0)

  #beta <- c(0.5, 0.3, 0.3, 0.25, 0.25, rep(0, n.noise.predictors))
  beta <- c(1, rep(0, n.noise.predictors))

  y <- runif(NN) < 1 / (1 + exp(-beta.0 - x %*% beta))

  DATA <- data.frame(x)
  DATA$y <- y * 1
  DATA
} # end generate_data_1_predictor_for_paper

generate_data <- function(
  NN,
  n.true.predictors = n.true.predictors,
  n.noise.predictors = n.noise.predictors,
  cor0 = 0.1,
  cor1 = 0.05,
  beta.0 = 0
) {
  n.predictors <- n.true.predictors + n.noise.predictors
  mu0 <- rep(0, n.predictors)

  # Specify correlation matrix
  Sigma0 <- matrix(0, nrow = n.predictors, ncol = n.predictors)
  Sigma0[1:n.true.predictors, 1:n.true.predictors] <- cor0
  Sigma0[
    (n.true.predictors + 1):n.predictors,
    (n.true.predictors + 1):n.predictors
  ] <- cor1
  diag(Sigma0) <- 1.0

  x <- mvrnorm(NN, mu0, Sigma0)

  beta <- c(0.5, 0.3, 0.3, 0.25, 0.25, rep(0, n.noise.predictors))

  y <- runif(NN) < 1 / (1 + exp(-beta.0 - x %*% beta))

  DATA <- data.frame(x)
  DATA$y <- y * 1
  DATA
} # end generate_data

val.perf <- function(data = DATA.val, fit = fit.full) {
  # Calibration slope model:
  fit.val       <- lrm(data$y~predict(fit, data))
  cal.slope     <- as.numeric(coef(fit.val)[2])
  # Calibration-in-the-large model:
  fit.cal       <- glm(data$y~offset(predict(fit, data)), family = 'binomial')
  cal.intercept <- as.numeric(coef(fit.cal)[1])
  c.val         <- as.numeric(fit.val$stats['C'])
  R2.val        <- as.numeric(fit.val$stats['R2'])
  return(list(c.val = c.val, 
              R2.val = R2.val, 
              cal.intercept = cal.intercept,
              cal.slope = cal.slope))
} # end val.perf

# apply the shrinkage factor to the model coefficients and re-estimate the intercept
# NOTE: the original function will not work because of code breaking changes in the rms package. I will rewrite the function to work with rms 8.2-1. 
# recalibrate  <- function(shrinkage, fit, data = DATA){
#   shrunk.coef   <- shrinkage * coef(fit)
#   fit.formula   <- formula(fit)
#   fit.shrunk    <- lrm(fit.formula, data = data, init = shrunk.coef, maxit = 1)
#   lrm.offset    <- lrm.fit(y = data$y, offset = predict(fit.shrunk, data, type = 'lp'))
#   new.intercept <- fit.shrunk$coefficients[1] + lrm.offset$coefficients[1]
#   fit.shrunk    <- lrm(fit.formula, data = data, init = c(new.intercept, shrunk.coef[-1]), maxit = 1)
#   return(fit.shrunk)
# }

# Rewritten to work with rms 8.2-1
recalibrate <- function(shrinkage, fit, data = DATA) {
  b <- coef(fit)
  shrunk.coef <- shrinkage * b

  fit.shrunk <- fit
  fit.shrunk$coefficients <- shrunk.coef
  lp.shrunk <- as.vector(predict(fit.shrunk, newdata = data, type = "lp"))

  fit.off <- rms::lrm.fit(y = data$y, offset = lp.shrunk)
  new.intercept <- unname(shrunk.coef[1] + coef(fit.off)[1])

  fit.final <- fit
  fit.final$coefficients <- c(new.intercept, shrunk.coef[-1])
  names(fit.final$coefficients) <- names(b)
  fit.final$linear.predictors <-
    as.vector(predict(fit.final, newdata = data, type = "lp"))
  fit.final
}

# This base-R version doesn't depend on any rms internals:
recalibrate <- function(shrinkage, fit, data = DATA) {
  b  <- coef(fit)
  lp <- shrinkage * as.vector(predict(fit, newdata = data, type = "lp"))
  c0 <- unname(coef(glm(data$y ~ offset(lp), family = binomial))[1])
  fit.final <- fit
  fit.final$coefficients <- c(c0 + 0, shrinkage * b[-1])
  fit.final$coefficients[1] <- shrinkage * b[1] + c0
  names(fit.final$coefficients) <- names(b)
  fit.final
}

full_model   <- function(data.d = DATA, B = B, data.v = DATA.val, n.cal.grid = 100) {
  ### fit maximum likelihood full model (no variable selection)
  fit.full <- lrm(y ~ ., data = data.d, x = T, y = T)
  
  ### get apparent performance of full model (full.apparent)
  apparent <- as.numeric(val.perf(data = data.d, fit = fit.full))
  val      <- as.numeric(val.perf(data = data.v, fit = fit.full))
  
  # bootstrap internal validation of the full model to get the optimism-corrected bootstrap
  # calibration slope
  val.full  <- rms::validate(fit.full, bw = F, B = B, pr = F, estimates = F)
  xx        <- c(as.vector(val.full[c("Dxy", "R2", "Intercept","Slope"), 'index.corrected']))
  boot.opt  <- c((1 +xx[1])/2, xx[-1])
  
  ## apply the bootstrap corrected slope shrinkage factor (and re-estimate the intercept)
  fit.shrunk <- recalibrate(shrinkage = boot.opt[4], fit = fit.full, data = data.d)
  shrunk.apparent <- as.numeric(val.perf(data = data.d, fit = fit.shrunk))
  shrunk.val      <- as.numeric(val.perf(data = data.v, fit = fit.shrunk))
  
  pred.val <- predict(fit.full, newdata = data.v, type = 'fitted')
  Sm       <- lowess(pred.val, data.v$y, iter = 0)
  pp.full  <- seq(min(pred.val), max(pred.val), length = n.cal.grid)
  Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y
  
  pred.shrunk <- predict(fit.shrunk, newdata = data.v, type = 'fitted')
  Sm          <- lowess(pred.shrunk, data.v$y, iter = 0)
  pp.shrunk   <- seq(min(pred.shrunk), max(pred.shrunk), length = n.cal.grid)
  Sm.shrunk   <- approx(Sm, xout = pp.shrunk, ties = mean)$y
  
  return(list(apparent        = apparent, 
              val             = val, 
              opt.cor         = boot.opt,
              shrunk.apparent = shrunk.apparent,
              shrunk.val      = shrunk.val,
              #fit.cal.x       = Sm.full$x,
              #fit.cal.y       = Sm.full$y,
              #shrunk.cal.x    = Sm.shrunk$x,
              #shrunk.cal.y    = Sm.shrunk$y,
              fit.cal.x       = pp.full,
              fit.cal.y       = Sm.full,
              shrunk.cal.x    = pp.shrunk,
              shrunk.cal.y    = Sm.shrunk,
              fit             = fit.full,
              fit.shrunk      = fit.shrunk, 
              shrinkage.f     = boot.opt[4]))
} # end full_model

EN_model <- function(
  data.d = data,
  data.v = data.v,
  mixing = 1,
  n.cal.grid = 100
) {
  P <- ncol(data.d) - 1
  fit.glmnet.elasticnet <- cv.glmnet(
    x = as.matrix(data.d[, 1:P]),
    y = data.d$y,
    family = 'binomial',
    alpha = mixing,
    nfolds = 5,
    type.measure = 'deviance',
    parallel = F
  )

  fit.full <- lrm(y ~ ., data = data.d, x = T, y = T)
  ### Pull out lambda.min
  cv.elasticnet <- fit.glmnet.elasticnet$lambda.min

  ### Get betas for lambda.min
  optimal.fullta.elasticnet <- as.numeric(predict(
    fit.glmnet.elasticnet,
    type = 'coefficients',
    s = "lambda.min"
  ))

  ### create lrm object with coefficients from the glmnet
  # fit.elasticnet <- lrm(
  #   y ~ .,
  #   data = data.d,
  #   init = optimal.fullta.elasticnet,
  #   maxit = 1
  # )
  # UPDATE: for rms 8.2-1. replace with a lrm object carrying the glmnet coefficients (replaces lrm(init=, maxit=1))
  fit.elasticnet <- fit.full
  b <- setNames(as.numeric(predict(fit.glmnet.elasticnet, type = "coefficients", s = "lambda.min")),
              c("Intercept", colnames(data.d)[1:P]))
  stopifnot(setequal(names(b), names(coef(fit.full))))
  fit.elasticnet$coefficients <- b[names(coef(fit.full))]
  fit.elasticnet$linear.predictors <-
    as.vector(predict(fit.elasticnet, newdata = data.d, type = "lp"))

  if (all(coef(fit.elasticnet)[-1] == 0)) {
    apparent <- rep(NA, 4)
    val <- rep(NA, 4)
    return(list(
      apparent = apparent,
      val = val,
      EN.cal.x = rep(NA, n.cal.grid),
      EN.cal.y = rep(NA, n.cal.grid)
    ))
  } else {
    ### Get various performance measures (apparent)
    apparent <- as.numeric(val.perf(data = data.d, fit = fit.elasticnet))
    val <- as.numeric(val.perf(data = data.v, fit = fit.elasticnet))

    pred <- predict(fit.elasticnet, newdata = data.v, type = 'fitted')
    Sm <- lowess(pred, data.v$y, iter = 0)
    pp.EN <- seq(min(pred), max(pred), length = 100)
    Sm.EN <- approx(Sm, xout = pp.EN, ties = mean)$y

    return(list(
      apparent = apparent,
      val = val,
      EN.cal.x = pp.EN,
      EN.cal.y = Sm.EN,
      fit = fit.full,
      fit.shrunk = fit.elasticnet
    ))
  }
} # end EN_model

NBOOT              <- 200
n.true.predictors  <- 1
n.noise.predictors <- 10
N.VAL              <- 5000
N.POP              <- 1000000

# fit model to a 'super' population (large sample model)
DATA.POP <- generate_data_1_predictor_for_paper(NN = N.POP, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors)
fit.pop  <- lrm(y ~ ., data = DATA.POP, x = T, y = T)
phi      <- as.numeric(prop.table(table(DATA.POP$y))[2])
R2_CS    <- as.numeric(fit.pop$stats['R2']) * (1 - (phi^phi * (1 - phi)^(1 - phi))^2)

P <- ncol(DATA.POP) - 1
mixing <- 1

fit.glmnet.elasticnet <- cv.glmnet(x            = as.matrix(DATA.POP[,1:P]), 
                                   y            = DATA.POP$y, 
                                   family       = 'binomial', 
                                   alpha        = mixing, 
                                   nfolds       = 5,
                                   type.measure = 'deviance', 
                                   parallel     = F)
cv.elasticnet <- fit.glmnet.elasticnet$lambda.min
optimal.fullta.elasticnet <- as.numeric(predict(fit.glmnet.elasticnet, type = 'coefficients', s = "lambda.min"))
# Due to changes in the rms package, this next line does not work:
# fit.EN.pop <- lrm(y ~ ., data = DATA.POP, init = optimal.fullta.elasticnet, maxit = 1)
# Workaround for rms 8.2-1:
b.en <- as.numeric(predict(fit.glmnet.elasticnet, type = "coefficients", s = "lambda.min"))

fit.pop.ml <- rms::lrm(y ~ ., data = DATA.POP, x = TRUE, y = TRUE)

stopifnot(length(b.en) == length(coef(fit.pop.ml)))
fit.EN.pop <- fit.pop.ml
fit.EN.pop$coefficients <- setNames(b.en, names(coef(fit.pop.ml)))
fit.EN.pop$linear.predictors <-
  as.vector(predict(fit.EN.pop, newdata = DATA.POP, type = "lp"))

# calculate required sample for model development
xx       <- pmsampsize::pmsampsize(type       = 'b', 
                                   parameters = n.true.predictors+n.noise.predictors,
                                   prevalence = as.numeric(prop.table(table(DATA.POP$y))[2]),
                                   csrsquared   = R2_CS,
                                   shrinkage  = 0.9)

N.Riley <- xx$sample_size

DATA.val              <- generate_data_1_predictor_for_paper(NN = N.VAL, n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors)
N.SIM                 <- 500
n.cal.grid            <- 100
N                     <- sort(c(100, 200, 500, 1000, 5000, N.Riley))
p.full                <- array(0, dim = c(nrow(DATA.val), N.SIM, length(N)))
p.shrunk              <- array(0, dim = c(nrow(DATA.val), N.SIM, length(N)))
p.EN                  <- array(0, dim = c(nrow(DATA.val), N.SIM, length(N)))
cal.shrunk            <- matrix(ncol = length(N), nrow = N.SIM)
cal.EN                <- matrix(ncol = length(N), nrow = N.SIM)
full.val.cal.x        <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
full.val.cal.y        <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
full.val.shrunk.cal.x <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
full.val.shrunk.cal.y <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
EN.val.cal.x          <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
EN.val.cal.y          <- array(NA, dim = c(n.cal.grid, N.SIM, length(N)))
dca.full              <- array(NA, dim = c(100, N.SIM, length(N)))
dca.shrunk            <- array(NA, dim = c(100, N.SIM, length(N)))
dca.EN                <- array(NA, dim = c(100, N.SIM, length(N)))

pb <- progress::progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = length(N) * N.SIM, clear = FALSE, width = 60)

# NOTE: this takes several hours to run, so below I added a line to save the results in a .RData file. 
for(j in 1:length(N)){
  cat("\n loop = ", j)
  for(i in 1:N.SIM){
    pb$tick()
    DATA <- generate_data_1_predictor_for_paper(NN = N[j], n.true.predictors = n.true.predictors, n.noise.predictors = n.noise.predictors)
    
    # full model + shrinkage
    fit.full <- full_model(data.d = DATA, B = NBOOT, data.v = DATA.val)
    
    ### get the loess calibration on the validation data (full model)
    full.val.cal.x[,i,j] <- fit.full$fit.cal.x
    full.val.cal.y[,i,j] <- fit.full$fit.cal.y
    
    ### get the loess calibration on the validation data (shrunken model)
    full.val.shrunk.cal.x[,i,j] <- fit.full$shrunk.cal.x
    full.val.shrunk.cal.y[,i,j] <- fit.full$shrunk.cal.y
    
    if(length(class(fit.full$fit))>1){
      p.full[,i,j]    <- predict(fit.full$fit, newdata = DATA.val, type = 'fitted')
      p.shrunk[,i,j]  <- predict(fit.full$fit.shrunk, newdata = DATA.val, type = 'fitted')
      cal.shrunk[i,j] <- as.numeric(val.perf(data = DATA.val, fit = fit.full$fit.shrunk))[4]
      DAT.DCA         <- data.frame(pred.full   = predict(fit.full$fit,        newdata = DATA.val, type = 'fitted'),
                                    pred.shrunk = predict(fit.full$fit.shrunk, newdata = DATA.val, type = 'fitted'),
                                    y           = DATA.val$y)
      fit.dca.full     <- dca(y~pred.full, data = DAT.DCA)
      fit.dca.shrunk   <- dca(y~pred.shrunk, data = DAT.DCA)
      dca.full[,i,j]   <- fit.dca.full$dca$net_benefit[fit.dca.full$dca$label == 'pred.full']
      dca.shrunk[,i,j] <- fit.dca.shrunk$dca$net_benefit[fit.dca.shrunk$dca$label == 'pred.shrunk']
    } else {
      p.full[,i, j]        <- rep(NA, n.cal.grid)
      p.shrunk[,i, j] <- rep(NA, n.cal.grid)
      cal.shrunk[i, j]     <- NA
      dca.full[,i, j]      <- rep(NA, 100)
      dca.shrunk[,i, j]    <- rep(NA, 100)
    }
  
    # Elastic net
    fit.EN   <- EN_model(data.d = DATA, data.v = DATA.val)
    
    ### get the loess calibration on the validation data (elastic net model)
    EN.val.cal.x[,i,j] <- fit.EN$EN.cal.x
    EN.val.cal.y[,i,j] <- fit.EN$EN.cal.y
    
    if(length(class(fit.EN$fit))>1){
      p.EN[,i,j]   <- predict(fit.EN$fit, newdata = DATA.val, type = 'fitted')
      cal.EN[i,j]  <- as.numeric(val.perf(data = DATA.val, fit = fit.EN$fit.shrunk))[4]
      DAT.DCA      <- data.frame(pred.EN = predict(fit.EN$fit.shrunk, newdata = DATA.val, type = 'fitted'),y = DATA.val$y)
      fit.dca.EN   <- dca(y~pred.EN,     data = DAT.DCA)
      dca.EN[,i,j] <- fit.dca.EN$dca$net_benefit[fit.dca.EN$dca$label == 'pred.EN']
    } else {
      p.EN[,i, j]  <- rep(NA, n.cal.grid)
      cal.EN[i, j] <- NA
      dca.EN[,i,j] <- rep(NA, 100)
    }
  }
}

### get predictions of the 'population' model on the validation data
p.POP    <- predict(fit.pop, newdata = DATA.val, type = 'fitted')
p.POP.EN <- predict(fit.EN.pop, newdata = DATA.val, type = 'fitted')

# Save the results so that you don't have to run the simulation again (it takes several hours to run!):
save.image(file = "Figure_2_code.RData")

# find the randomly chosen individual with large sample predictions of 0.1:0.9[0.1]
which.individual <- c(which(as.numeric(round(p.POP, 2)) == 0.1)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.1)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.2)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.2)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.3)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.3)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.4)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.4)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.5)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.5)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.6)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.6)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.7)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.7)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.8)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.8)))[1]],
                      which(as.numeric(round(p.POP, 2)) == 0.9)[sample(1:length(which(as.numeric(round(p.POP, 2)) == 0.9)))[1]])

which.individual.EN <- c(which(as.numeric(round(p.POP.EN, 2)) == 0.1)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.1)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.2)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.2)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.3)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.3)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.4)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.4)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.5)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.5)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.6)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.6)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.7)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.7)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.8)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.8)))[1]],
                         which(as.numeric(round(p.POP.EN, 2)) == 0.9)[sample(1:length(which(as.numeric(round(p.POP.EN, 2)) == 0.9)))[1]])

# This code using gather() is outdated, see the tidyr package documentation for the new pivot_longer() function. I will update this code to use pivot_longer() instead of gather().
# OUT.orig.1 <- gather(as_tibble(cbind(p.full[which.individual[1],,1], 
#                                      p.full[which.individual[1],,2],
#                                      p.full[which.individual[1],,3], 
#                                      p.full[which.individual[1],,4], 
#                                      p.full[which.individual[1],,5], 
#                                      p.full[which.individual[1],,6])))

# Additionally, the objects OUT.orig.x and others follow a repeating structure, so we shall use help functions:
long_slice <- function(arr, i) {
  m <- arr[i, , ]
  colnames(m) <- paste0("V", seq_len(ncol(m)))
  tibble::as_tibble(m) |>
    tidyr::pivot_longer(everything(), names_to = "key", values_to = "value", cols_vary = "slowest")
}

long_slices <- function(arr, idx, prefix) {
  out <- lapply(idx, long_slice, arr = arr)
  setNames(out, paste0(prefix, ".", seq_along(idx)))
}

# Revised:
OUT.orig.list <- long_slices(p.full, which.individual, "OUT.orig")

# # Verify that the old code and new give the same result (on OUT.orig.1):
# old <- suppressWarnings(
#   tidyr::gather(as_tibble(cbind(p.full[which.individual[1], , 1], p.full[which.individual[1], , 2],
#                                 p.full[which.individual[1], , 3], p.full[which.individual[1], , 4],
#                                 p.full[which.individual[1], , 5], p.full[which.individual[1], , 6])))
# )
# stopifnot(isTRUE(all.equal(old, OUT.orig.list$OUT.orig.1)))  # they are identical

OUT.shrunk.list <- long_slices(p.shrunk, which.individual, "OUT.shrunk")

# old <- suppressWarnings(
#   tidyr::gather(as_tibble(cbind(p.shrunk[which.individual[1],,1], 
#                                        p.shrunk[which.individual[1],,2],
#                                        p.shrunk[which.individual[1],,3], 
#                                        p.shrunk[which.individual[1],,4], 
#                                        p.shrunk[which.individual[1],,5], 
#                                        p.shrunk[which.individual[1],,6])))
#   )
# stopifnot(isTRUE(all.equal(old, OUT.shrunk.list$OUT.shrunk.1)))  # they are identical


OUT.EN.list <- long_slices(p.EN, which.individual.EN, "OUT.EN")

# old <- suppressWarnings(
#   tidyr::gather(as_tibble(cbind(p.EN[which.individual.EN[1],,1], 
#                                    p.EN[which.individual.EN[1],,2],
#                                    p.EN[which.individual.EN[1],,3], 
#                                    p.EN[which.individual.EN[1],,4], 
#                                    p.EN[which.individual.EN[1],,5], 
#                                    p.EN[which.individual.EN[1],,6])))
#   )
# stopifnot(isTRUE(all.equal(old, OUT.EN.list$OUT.EN.1)))  # they are identical

OUT.orig <- dplyr::bind_rows(OUT.orig.list)

OUT.shrunk <- dplyr::bind_rows(OUT.shrunk.list)

OUT.EN <- dplyr::bind_rows(OUT.EN.list)

OUT.orig   <- tibble::add_column(OUT.orig,   model = rep("full model", nrow(OUT.orig)))
OUT.shrunk <- tibble::add_column(OUT.shrunk, model = rep("shrunk LR",  nrow(OUT.shrunk)))
OUT.EN     <- tibble::add_column(OUT.EN,     model = rep("lasso",      nrow(OUT.EN)))

OUT     <- dplyr::bind_rows(OUT.orig, OUT.shrunk, OUT.EN)
OUT$key <- factor(OUT$key, levels = paste("V", seq(1:length(N)), sep = ''), labels = as.character(N))
OUT     <- dplyr::mutate(OUT, model = factor(model, levels = c("full model", "shrunk LR", "lasso")))
OUT     <- dplyr::mutate(OUT, key = factor(key, levels = as.character(N)))
OUT     <- tibble::add_column(OUT, prob = c(rep("p=0.1", N.SIM * length(N)),
                                    rep("p=0.2", N.SIM * length(N)),
                                    rep("p=0.3", N.SIM * length(N)),
                                    rep("p=0.4", N.SIM * length(N)),
                                    rep("p=0.5", N.SIM * length(N)),
                                    rep("p=0.6", N.SIM * length(N)),
                                    rep("p=0.7", N.SIM * length(N)),
                                    rep("p=0.8", N.SIM * length(N)),
                                    rep("p=0.9", N.SIM * length(N)),
                                    rep("p=0.1", N.SIM * length(N)),
                                    rep("p=0.2", N.SIM * length(N)),
                                    rep("p=0.3", N.SIM * length(N)),
                                    rep("p=0.4", N.SIM * length(N)),
                                    rep("p=0.5", N.SIM * length(N)),
                                    rep("p=0.6", N.SIM * length(N)),
                                    rep("p=0.7", N.SIM * length(N)),
                                    rep("p=0.8", N.SIM * length(N)),
                                    rep("p=0.9", N.SIM * length(N)),
                                    rep("p=0.1", N.SIM * length(N)),
                                    rep("p=0.2", N.SIM * length(N)),
                                    rep("p=0.3", N.SIM * length(N)),
                                    rep("p=0.4", N.SIM * length(N)),
                                    rep("p=0.5", N.SIM * length(N)),
                                    rep("p=0.6", N.SIM * length(N)),
                                    rep("p=0.7", N.SIM * length(N)),
                                    rep("p=0.8", N.SIM * length(N)),
                                    rep("p=0.9", N.SIM * length(N))))
OUT <- dplyr::mutate(OUT, prob = factor(prob))

OUT_mean <- dplyr::filter(OUT, key == '5000') |> dplyr::group_by(prob) |> dplyr::summarize(mean_val = mean(value, na.rm = T))
OUT_mean$mean_val <- seq(0.1, 0.9, 0.1)

# The original code's dplyr steps left the output still grouped by prob, which I don't think is the intent. 
# In revising the code, bear in mind, group_by() sorts the output by the grouping keys. .by keeps groups in order of first appearance in the data. If a plot or comparison assumes sorted rows (which it may in this case -- would have to determine that via further testing), add arrange(prob, key)
# x1 <- dplyr::filter(OUT, model == 'full model') |> dplyr::group_by(prob, key) |> dplyr::summarize(sd = sd(value, na.rm = T))
# x2 <- dplyr::filter(OUT, model == 'shrunk LR') |> dplyr::group_by(prob, key) |> dplyr::summarize(sd = sd(value, na.rm = T))
# x3 <- dplyr::filter(OUT, model == 'lasso') |> dplyr::group_by(prob, key) |> dplyr::summarize(sd = sd(value, na.rm = T))

x1 <- OUT |>
  dplyr::filter(model == "full model") |>
  dplyr::summarize(sd = sd(value, na.rm = TRUE), .by = c(prob, key)) |>
  dplyr::arrange(prob, key) |>
  as.data.frame()

x2 <- OUT |>
  dplyr::filter(model == "shrunk LR") |>
  dplyr::summarize(sd = sd(value, na.rm = TRUE), .by = c(prob, key)) |>
  dplyr::arrange(prob, key) |>
  as.data.frame()

x3 <- OUT |>
  dplyr::filter(model == "lasso") |>
  dplyr::summarize(sd = sd(value, na.rm = TRUE), .by = c(prob, key)) |>
  dplyr::arrange(prob, key) |>
  as.data.frame()

# Don't need these lines anymore:
# x1 <- as.data.frame(x1)
# x2 <- as.data.frame(x2)
# x3 <- as.data.frame(x3)
# colnames(x1) <- c("prob", "key", "sd")
# colnames(x2) <- c("prob", "key", "sd")
# colnames(x3) <- c("prob", "key", "sd")

xx.sd <- cbind(x1, x2$sd, x3$sd)
colnames(xx.sd)<- c("prob", "N", "sd(full model)", "sd(shrunk)", "sd(lasso)")

### plot calibration
########################################### REPLACE THIS CODE BLOCK WITH A FUNCTION TO REDUCE REPETITION AND IMPROVE MAINTAINABILITY
# NOTES: reshape2::melt() on an array turns every cell into one row. The columns are an index for each dimension, plus the cell's value. full.val.cal.x for example is of dim [100, 500, 6]
# Var1	Index along dimension 1 (1 to 100), the position on the calibration grid
# Var2	Index along dimension 2 (1 to 500), the replicate
# Var3	Index along dimension 3 (1 to 6), the sixth dimension, probably the sample-size or scenario level
# value	The array entry at [Var1, Var2, Var3]
# These names are defaults.
# reshape2 is retired.

# cal.full.xx <- reshape2::melt(full.val.cal.x)
# cal.full.yy <- reshape2::melt(full.val.cal.y)

# The base-R equivalent is:
# cal.full.xx <- data.frame(
#   as.data.frame(arrayInd(seq_along(full.val.cal.x), dim(full.val.cal.x)))  |> setNames(c("Var1", "Var2", "Var3")),
#   value = as.vector(full.val.cal.x)
# )

# cal.full.yy <- data.frame(
#   as.data.frame(arrayInd(seq_along(full.val.cal.y), dim(full.val.cal.y)))  |> setNames(c("Var1", "Var2", "Var3")),
#   value = as.vector(full.val.cal.y)
# )

# a tidyverse version is not as elegant:
# dd <- dim(full.val.cal.x)
# cal.full.xx.3 <- tidyr::expand_grid(
#   Var1 = seq_len(dd[1]),
#   Var2 = seq_len(dd[2]),
#   Var3 = seq_len(dd[3]),
#   .vary = "fastest"
# ) |>
#   dplyr::mutate(value = as.vector(full.val.cal.x)) |>
#   as.data.frame()

# cal.shrunk.xx <- reshape2::melt(full.val.shrunk.cal.x)
# cal.shrunk.xx <- data.frame(
#   as.data.frame(arrayInd(seq_along(full.val.shrunk.cal.x), dim(full.val.shrunk.cal.x)))  |> setNames(c("Var1", "Var2", "Var3")),
#   value = as.vector(full.val.shrunk.cal.x)
# )

# cal.shrunk.yy <- reshape2::melt(full.val.shrunk.cal.y)
# cal.shrunk.yy <- data.frame(
#   as.data.frame(arrayInd(seq_along(full.val.shrunk.cal.y), dim(full.val.shrunk.cal.y)))  |> setNames(c("Var1", "Var2", "Var3")),
#   value = as.vector(full.val.shrunk.cal.y)
# )

# cal.EN.xx <- reshape2::melt(EN.val.cal.x)
# cal.EN.xx <- data.frame(
#   as.data.frame(arrayInd(seq_along(EN.val.cal.x), dim(EN.val.cal.x)))  |> setNames(c("Var1", "Var2", "Var3")),
#   value = as.vector(EN.val.cal.x)
# )

# cal.EN.yy <- reshape2::melt(EN.val.cal.y)
# cal.EN.yy <- data.frame(
#   as.data.frame(arrayInd(seq_along(EN.val.cal.y), dim(EN.val.cal.y)))  |> setNames(c("Var1", "Var2", "Var3")),
#   value = as.vector(EN.val.cal.y)
# )

# cal.full.xx$method   <- rep("full model", nrow(cal.full.xx))
# cal.shrunk.xx$method <- rep("shrunk LR",   nrow(cal.shrunk.xx))
# cal.EN.xx$method     <- rep("lasso", nrow(cal.EN.xx))

# cal <- rbind(cal.full.xx,   
#              cal.shrunk.xx,
#              cal.EN.xx)

# colnames(cal) <- c("Var1", "Sim", "N", "value", "method")
# cal$method    <- factor(cal$method, levels = c('full model', 'shrunk LR', 'lasso'))
# cal$N         <- factor(cal$N, levels = seq(1:length(N)), labels = paste("N=", N, sep=''))

# These next two lines are not needed. Perhaps used to avoid copying a moderately large object in memory when adding a single column. No actual benefit here, because data.table(cal) makes a copy during conversion anyway. Modern R uses copy-on-modify semantics at the column level, so adding a single column to a data frame does not duplicate the existing columns.
# OUT.cal <- data.table(cal)
# OUT.cal[, y:=c(cal.full.yy$value,   
#                cal.shrunk.yy$value,
#                cal.EN.yy$value)]

# You can instead use:
# cal$y <- c(cal.full.yy$value,
#            cal.shrunk.yy$value,
#            cal.EN.yy$value)
###########################################

# The entire code block above for preparing the calibration data can be simplified (following DRY - don't repeat yourself) by creating a helper function to handle the repeated logic for each model type. This would reduce redundancy and improve maintainability.

calibration_df <- function(
  # Prepare calibration curves for plotting
  x,
  y,
  method_name,
  sample_sizes,
  method_levels = c("full model", "shrunk LR", "lasso")
) {
  # The dimension check catches mismatched array shapes before they are combined.
  stopifnot(
    length(dim(x)) == 3L,
    identical(dim(x), dim(y)),
    all(dim(x) > 0L),
    dim(x)[3] == length(sample_sizes),
    !anyNA(sample_sizes),
    anyDuplicated(sample_sizes) == 0L,
    length(method_name) == 1L,
    method_name %in% method_levels
  )

  # Dimensions: calibration-grid position × replicate × sample-size scenario
  out <- as.data.frame(arrayInd(seq_along(x), dim(x)))
  names(out) <- c("Var1", "Sim", "N")

  out$N <- factor(
    out$N,
    levels = seq_along(sample_sizes),
    labels = paste0("N=", sample_sizes)
  )

  out$value <- as.vector(x)
  out$method <- factor(
    rep(method_name, nrow(out)),
    levels = method_levels
  )
  out$y <- as.vector(y)

  out
}

# Each method’s x- and y-coordinates stay together. 
cal <- rbind(
  calibration_df(
    x = full.val.cal.x,
    y = full.val.cal.y,
    method_name = "full model",
    sample_sizes = N
  ),
  calibration_df(
    x = full.val.shrunk.cal.x,
    y = full.val.shrunk.cal.y,
    method_name = "shrunk LR",
    sample_sizes = N
  ),
  calibration_df(
    x = EN.val.cal.x,
    y = EN.val.cal.y,
    method_name = "lasso",
    sample_sizes = N
  )
)
  
p3 <- 
  cal |>
  dplyr::filter(method != 'shrunk LR') %>%
  ggplot(aes(x = value, y = y, group = Sim)) +
  geom_line(alpha = 0.05, colour = 'darkgrey') +
  facet_grid(method ~ N) +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab('estimated risk') +
  ylab('observed') +
  geom_abline(slope = 1, intercept = 0, colour = 'black') +
  theme_bw() +
  theme(axis.text = element_text(size = 6))
p3

