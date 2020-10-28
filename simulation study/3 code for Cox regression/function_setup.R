#####################################################################
###### functions to generate data set from simulation setup ########
#####################################################################

# correlation matrix

sigma <- matrix(0, 15, 15)
diag(sigma) <- 1

sigma[1, 2] <- 0.8
sigma[1, 9] <- 0.3
sigma[3, 9] <- -0.5
sigma[8, 9] <- -0.3
sigma[3, 5] <- 0.3
sigma[8, 11] <- 0.3
sigma[4, 7] <- -0.3
sigma[6, 9] <- 0.3
sigma[4, 6] <- -0.5
sigma[5, 6] <- -0.3
sigma[5, 12] <- 0.5
sigma[6, 7] <- 0.5
sigma[7, 11] <- 0.3
sigma[6, 11] <- 0.5
sigma[6, 14] <- 0.3
sigma[7, 14] <- 0.3
sigma[11, 14] <- 0.5

for (i in 1:14) {
  for (j in (i + 1):15)
    sigma[j, i] <- sigma[i, j]
}


# generate z from a multivariate normal distribution

gen.z <- function(n = 50000, s = sigma) {
  z <- rmvnorm(n, mean = rep(0, 15), s)
  return(z)
}

# generate x from z

gen.x.from.z <- function(z) {
  x <- z * 0
  x[, 1] <- floor(10 * z[, 1] + 55)
  x[, 2] <- (z[, 2] < 0.6) * 1
  x[, 3] <- exp(0.4 * z[, 3] + 3)
  x[, 4] <- (z[, 4] > -1.2) + (z[, 4] >= 0.75)
  x[, 5] <- exp(0.5 * z[, 5] + 1.5)
  x[, 6] <- floor(apply(cbind(100 * exp(z[, 6]) - 20, 0), 1, max))
  x[, 7] <- floor(apply(cbind(80 * exp(z[, 7]) - 20, 0), 1, max))
  x[, 8] <- (z[, 8] < -0.35)
  x[, 9] <- ((z[, 9] >= 0.5) & (z[, 9] < 1.5)) + 2 * (z[, 9] >= 1.5)
  x[, 10] <- 0.01 * (100 * (z[, 10] + 4) ** 2)
  x[, 11] <- floor(10 * z[, 11] + 55)
  x[, 12] <- floor(10 * z[, 12] + 55)
  x[, 13] <- floor(10 * z[, 13] + 55)
  x[, 14] <- (z[, 14] < 0)
  x[, 15] <- (z[, 15] < 0)
  return(x)
}


# partial predictors

partial.predictor <- function(x) {
  pred <- x * 0
  pred[, 1] <- 3.5 * x[, 1] ^ 0.5 - 0.25 * x[, 1]
  pred[, 3] <- 2 * (log((x[, 3] + 10) / 25)) ^ 2
  pred[, 4] <- -0.4 * (x[, 4] == 1)
  pred[, 5] <- -(0.15 * x[, 5] + 0.75 * exp(-((log(
    x[, 5]
  ) - 1.5) ^ 2) / 4))
  pred[, 6] <- 0.25 * log(x[, 6] + 1)
  pred[, 8] <- 0.4 * x[, 8]
  pred[, 10] <- 0.021 * x[, 10]
  pred[, 11] <- 0.04 * x[, 11]
  return(pred)
}


# transform some x

pretrans.x <- function(x = x) {
  x.t <- x
  x.t[, 3] <- log(x[, 3])
  x.t[, 5] <- log(x[, 5] + 1.5)
  x.t[, 6] <- log(x[, 6] + 1)
  return(x.t)
}


# create design variables

x.to.data <- function(x, y) {
  outdat <-
    data.frame(
      y = y,
      x1 = x[, 1],
      x2 = x[, 2],
      x3 = x[, 3],
      x4a = 1 * (x[, 4] == 1),
      x4b = 1 * (x[, 4] == 2),
      x5 = x[, 5],
      x6 = x[, 6],
      x7 = x[, 7],
      x8 = x[, 8],
      x9a = 1 * (x[, 9] == 1),
      x9b = 1 * (x[, 9] == 2),
      x10 = x[, 10],
      x11 = x[, 11],
      x12 = x[, 12],
      x13 = x[, 13],
      x14 = x[, 14],
      x15 = x[, 15]
    )
  return(outdat)
}