##############################################################
# LPIgamma
# Classic CL Index
# Formula: CL = (mu - L) / sigma
##############################################################

###############################################################
# 1. Gamma MLE (shape alpha, scale theta)
###############################################################

LPIgamma_mle_helper <- function(x){

  fit <- MASS::fitdistr(x, densfun="gamma")

  alpha_hat <- fit$estimate["shape"]
  theta_hat <- 1 / fit$estimate["rate"]

  list(
    shape = alpha_hat,
    scale = theta_hat,
    vcov  = fit$vcov
  )
}

###############################################################
# 2. Point Estimation
###############################################################

LPIgamma_point_est <- function(x, L){

  out <- LPIgamma_mle_helper(x)

  a  <- out$shape
  th <- out$scale

  mu  <- a * th
  sig <- th * sqrt(a)

  CL <- (mu - L) / sig

  list(Value=CL, mu=mu, sigma=sig,
       shape=a, scale=th, vcov=out$vcov)
}

###############################################################
# 3. Asymptotic CI (Delta Method)
###############################################################

LPIgamma_asymptotic <- function(x, L, alpha=0.05){

  out <- LPIgamma_point_est(x, L)

  a   <- out$shape
  th  <- out$scale
  mu  <- out$mu
  sig <- out$sigma
  V   <- out$vcov

  IndexVal <- out$Value

  dmu_da  <- th
  dmu_dth <- a

  dsig_da  <- th/(2*sqrt(a))
  dsig_dth <- sqrt(a)

  dCL_dmu  <- 1/sig
  dCL_dsig <- -(mu-L)/sig^2

  dCL_da  <- dCL_dmu*dmu_da  + dCL_dsig*dsig_da
  dCL_dth <- dCL_dmu*dmu_dth + dCL_dsig*dsig_dth

  grad <- c(dCL_da, dCL_dth)

  varCL <- as.numeric(t(grad) %*% V %*% grad)
  seCL  <- sqrt(varCL)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seCL,
          IndexVal + qnorm(1-alpha/2)*seCL)

  list(Value=IndexVal, SE=seCL,
       LCL=CI[1], UCL=CI[2])
}

###############################################################
# 4. Bootstrap
###############################################################

LPIgamma_bootstrap <- function(x, L, R=2000, alpha=0.05){

  boot_stat <- function(data, indices){
    xx <- data[indices]
    LPIgamma_point_est(xx, L)$Value
  }

  b <- boot::boot(data=x, statistic=boot_stat, R=R)

  ci <- tryCatch(
    boot::boot.ci(b, conf=1-alpha,
                  type=c("norm","basic","perc","bca")),
    error=function(e) NULL
  )

  norm_bounds  <- c(NA, NA)
  basic_bounds <- c(NA, NA)
  perc_bounds  <- c(NA, NA)
  bca_bounds   <- c(NA, NA)

  if(!is.null(ci)){
    if(!is.null(ci$normal))  norm_bounds  <- ci$normal[2:3]
    if(!is.null(ci$basic))   basic_bounds <- ci$basic[4:5]
    if(!is.null(ci$percent)) perc_bounds  <- ci$percent[4:5]
    if(!is.null(ci$bca))     bca_bounds   <- ci$bca[4:5]
  }

  list(
    Value = b$t0,
    norm  = norm_bounds,
    basic = basic_bounds,
    perc  = perc_bounds,
    bca   = bca_bounds
  )
}

###############################################################
# 5. Nonparametric
###############################################################

LPIgamma_nonparam <- function(x, L, alpha=0.05){

  mu_hat  <- mean(x)
  sig_hat <- sd(x)

  CL <- (mu_hat - L)/sig_hat

  n <- length(x)

  se_mu  <- sig_hat/sqrt(n)
  se_sig <- sig_hat*sqrt(1/(2*(n-1)))

  dCL_dmu  <- 1/sig_hat
  dCL_dsig <- -(mu_hat-L)/sig_hat^2

  varCL <- dCL_dmu^2*se_mu^2 +
    dCL_dsig^2*se_sig^2

  seCL <- sqrt(varCL)

  CI <- c(CL - qnorm(1-alpha/2)*seCL,
          CL + qnorm(1-alpha/2)*seCL)

  list(Value=CL, SE=seCL,
       LCL=CI[1], UCL=CI[2])
}

###############################################################
# 6. MAIN FUNCTION
###############################################################

#' LPI Gamma Classic CL Index
#'
#' Computes the classic capability index under the Gamma model.
#'
#' @param x numeric vector
#' @param L lower bound
#' @param R number of bootstrap replications
#' @param alpha significance level
#' @return data frame of results
#' @export
LPIgamma <- function(x, L, R=1000, alpha=0.05){

  mle    <- LPIgamma_point_est(x, L)
  asym   <- LPIgamma_asymptotic(x, L, alpha)
  bootci <- LPIgamma_bootstrap(x, L, R, alpha)
  nonpar <- LPIgamma_nonparam(x, L, alpha)

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
      bootci$norm[1],  bootci$norm[2],
      bootci$basic[1], bootci$basic[2],
      bootci$perc[1],  bootci$perc[2],
      bootci$bca[1],   bootci$bca[2],
      nonpar$Value, nonpar$LCL, nonpar$UCL
    )
  )

  invisible(TAB)
}
