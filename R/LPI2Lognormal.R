##############################################################
# LPI2LogNormal
# Formula: (Q_median - L) / (Q3 - Q1)
##############################################################

###############################################################
# Internal: Lognormal MLE Helper
###############################################################

LPI2LogNormal_mle_helper <- function(x){

  lx <- log(x)
  mu_hat <- mean(lx)
  sigma_hat <- sqrt(mean((lx - mu_hat)^2))

  n <- length(x)

  V <- matrix(c(
    sigma_hat^2/n, 0,
    0, sigma_hat^2/(2*n)
  ), 2, 2, byrow=TRUE)

  list(mu=mu_hat, sigma=sigma_hat, vcov=V)
}

###############################################################
# Point Estimation
###############################################################

LPI2LogNormal_point_est <- function(x, p1=0.25, p3=0.75, L){

  out <- LPI2LogNormal_mle_helper(x)
  mu <- out$mu
  s  <- out$sigma

  Q1   <- exp(mu + s*qnorm(p1))
  Qmed <- exp(mu)
  Q3   <- exp(mu + s*qnorm(p3))

  Val <- (Qmed - L)/(Q3 - Q1)

  list(Value=Val, mu=mu, sigma=s, vcov=out$vcov)
}

###############################################################
# Asymptotic CI
###############################################################

LPI2LogNormal_asymptotic <- function(x, p1=0.25, p3=0.75, L, alpha=0.05){

  out <- LPI2LogNormal_point_est(x, p1, p3, L)
  mu <- out$mu
  s  <- out$sigma
  V  <- out$vcov
  IndexVal <- out$Value

  z1 <- qnorm(p1)
  z3 <- qnorm(p3)

  Q1   <- exp(mu + s*z1)
  Qmed <- exp(mu)
  Q3   <- exp(mu + s*z3)

  Range <- Q3 - Q1

  dQ_dmu    <- function(Q) Q
  dQ_dsigma <- function(Q, z) Q * z

  dI_dQmed <- 1/Range
  dI_dQ1   <- (Qmed - L)/Range^2
  dI_dQ3   <- -(Qmed - L)/Range^2

  dI_dmu <- dI_dQ1*dQ_dmu(Q1) +
    dI_dQmed*dQ_dmu(Qmed) +
    dI_dQ3*dQ_dmu(Q3)

  dI_ds  <- dI_dQ1*dQ_dsigma(Q1,z1) +
    dI_dQ3*dQ_dsigma(Q3,z3)

  grad <- c(dI_dmu, dI_ds)

  varI <- t(grad) %*% V %*% grad
  seI  <- sqrt(varI)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# Bootstrap
###############################################################

LPI2LogNormal_bootstrap <- function(x, R=2000, p1=0.25, p3=0.75, L, alpha=0.05){

  boot_stat <- function(data, indices){
    xx <- data[indices]
    LPI2LogNormal_point_est(xx, p1, p3, L)$Value
  }

  b <- boot::boot(data=x, statistic=boot_stat, R=R)
  ci <- boot::boot.ci(b, conf=(1-alpha),
                      type=c("norm","basic","perc","bca"))

  list(
    Value = b$t0,
    norm  = ci$normal[2:3],
    basic = ci$basic[4:5],
    perc  = ci$percent[4:5],
    bca   = ci$bca[4:5]
  )
}

###############################################################
# Nonparametric CI
###############################################################

LPI2LogNormal_nonparam <- function(x, p1=0.25, p3=0.75, L, alpha=0.05){

  n <- length(x)

  q1   <- quantile(x, p1, type=8)
  qmed <- quantile(x, 0.5, type=8)
  q3   <- quantile(x, p3, type=8)

  Range <- q3 - q1
  IndexVal <- (qmed - L)/Range

  seI <- sd(x)/sqrt(n)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# MAIN FUNCTION
###############################################################

#' LPI2 Log-Normal Index
#'
#' Computes the LPI2 index under the Log-Normal distribution.
#'
#' @param x numeric vector
#' @param L lower bound
#' @param p1 first percentile
#' @param p3 second percentile
#' @param R number of bootstrap replications
#' @param alpha significance level
#' @return A data frame of estimation results
#' @export
LPI2LogNormal <- function(x, L, p1=0.25, p3=0.75,
                          R=1000, alpha=0.05){

  mle    <- LPI2LogNormal_point_est(x, p1, p3, L)
  asym   <- LPI2LogNormal_asymptotic(x, p1, p3, L, alpha)
  bootci <- LPI2LogNormal_bootstrap(x, R, p1, p3, L, alpha)
  nonpar <- LPI2LogNormal_nonparam(x, p1, p3, L, alpha)

  TAB <- data.frame(
    Item = c("MLE_Value",
             "Asymp_LCL","Asymp_UCL",
             "Boot_norm_LCL","Boot_norm_UCL",
             "Boot_basic_LCL","Boot_basic_UCL",
             "Boot_perc_LCL","Boot_perc_UCL",
             "Boot_bca_LCL","Boot_bca_UCL",
             "Nonpar_Value","Nonpar_LCL","Nonpar_UCL"),
    Value = c(
      mle$Value,
      asym$LCL, asym$UCL,
      bootci$norm[1], bootci$norm[2],
      bootci$basic[1], bootci$basic[2],
      bootci$perc[1], bootci$perc[2],
      bootci$bca[1], bootci$bca[2],
      nonpar$Value, nonpar$LCL, nonpar$UCL
    )
  )

  invisible(TAB)
}
