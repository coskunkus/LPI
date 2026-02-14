##############################################################
# LPI2Gamma
# Distribution: Gamma (shape=k, scale=theta)
# Formula: (Q_median - L) / (Q3 - Q1)
##############################################################

###############################################################
# 1. Internal: Gamma MLE Helper
###############################################################

LPI2Gamma_mle_helper <- function(x){

  nll <- function(par){
    k <- par[1]
    lam <- par[2]
    if(k <= 0 || lam <= 0) return(1e10)
    -sum(dgamma(x, shape=k, scale=lam, log=TRUE))
  }

  mean_x <- mean(x)
  var_x  <- var(x)

  start_scale <- var_x / mean_x
  start_shape <- mean_x^2 / var_x

  fit <- nlm(nll, p=c(start_shape, start_scale), hessian=TRUE)

  k_hat  <- fit$estimate[1]
  lam_hat <- fit$estimate[2]

  vcov_mat <- solve(fit$hessian)

  list(shape=k_hat, scale=lam_hat, vcov=vcov_mat)
}

###############################################################
# 2. Point Estimation
###############################################################

LPI2Gamma_point_est <- function(x, p1=0.25, p3=0.75, L){

  out <- LPI2Gamma_mle_helper(x)
  k   <- out$shape
  lam <- out$scale

  Q1   <- qgamma(p1, shape=k, scale=lam)
  Qmed <- qgamma(0.5, shape=k, scale=lam)
  Q3   <- qgamma(p3, shape=k, scale=lam)

  Val <- (Qmed - L)/(Q3 - Q1)

  list(Value=Val, shape=k, scale=lam, vcov=out$vcov)
}

###############################################################
# 3. Asymptotic CI (Delta Method)
###############################################################

LPI2Gamma_asymptotic <- function(x, p1=0.25, p3=0.75, L, alpha=0.05){

  out <- LPI2Gamma_point_est(x, p1, p3, L)

  k   <- out$shape
  lam <- out$scale
  V   <- out$vcov
  IndexVal <- out$Value

  Q1   <- qgamma(p1, shape=k, scale=lam)
  Qmed <- qgamma(0.5, shape=k, scale=lam)
  Q3   <- qgamma(p3, shape=k, scale=lam)

  Range <- Q3 - Q1

  dQ_dlambda <- function(Q, lam) Q/lam

  dQ_dk_num <- function(p, k, lam){
    eps <- 1e-5
    (qgamma(p, k+eps, lam) - qgamma(p, k-eps, lam))/(2*eps)
  }

  dQ1_dlam   <- dQ_dlambda(Q1, lam)
  dQmed_dlam <- dQ_dlambda(Qmed, lam)
  dQ3_dlam   <- dQ_dlambda(Q3, lam)

  dQ1_dk     <- dQ_dk_num(p1, k, lam)
  dQmed_dk   <- dQ_dk_num(0.5, k, lam)
  dQ3_dk     <- dQ_dk_num(p3, k, lam)

  dI_dQmed <- 1 / Range
  dI_dQ1   <- (Qmed - L) / Range^2
  dI_dQ3   <- -(Qmed - L) / Range^2

  dI_dk   <- dI_dQ1*dQ1_dk + dI_dQmed*dQmed_dk + dI_dQ3*dQ3_dk
  dI_dlam <- dI_dQ1*dQ1_dlam + dI_dQmed*dQmed_dlam + dI_dQ3*dQ3_dlam

  grad <- c(dI_dk, dI_dlam)

  varI <- t(grad) %*% V %*% grad
  seI  <- sqrt(varI)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 4. Bootstrap
###############################################################

LPI2Gamma_bootstrap <- function(x, R=2000, p1=0.25, p3=0.75, L, alpha=0.05){

  boot_stat <- function(data, indices){
    xx <- data[indices]
    LPI2Gamma_point_est(xx, p1, p3, L)$Value
  }

  b <- boot::boot(data=x, statistic=boot_stat, R=R)

  ci <- boot::boot.ci(b, conf=1-alpha,
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
# 5. Nonparametric CI
###############################################################

LPI2Gamma_nonparam <- function(x, p1=0.25, p3=0.75, L, alpha=0.05){

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
# 6. MAIN FUNCTION
###############################################################

#' LPI2 Gamma Index
#'
#' Computes the LPI2 index under Gamma distribution.
#'
#' @param x numeric vector
#' @param L lower bound
#' @param p1 first percentile
#' @param p3 second percentile
#' @param R number of bootstrap replications
#' @param alpha significance level
#' @return A data frame of estimation results
#' @export
LPI2Gamma <- function(x, L, p1=0.25, p3=0.75,
                      R=1000, alpha=0.05){

  mle    <- LPI2Gamma_point_est(x, p1, p3, L)
  asym   <- LPI2Gamma_asymptotic(x, p1, p3, L, alpha)
  bootci <- LPI2Gamma_bootstrap(x, R, p1, p3, L, alpha)
  nonpar <- LPI2Gamma_nonparam(x, p1, p3, L, alpha)

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
