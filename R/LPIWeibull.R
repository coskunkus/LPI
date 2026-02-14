##############################################################
# LPIWeibull
# Classic CL Index under Weibull model
# Formula: CL = (mu - L) / sigma
##############################################################

###############################################################
# 1. Weibull MLE
###############################################################

LPIWeibull_mle_helper <- function(x){

  nll <- function(par){
    k  <- par[1]
    lam <- par[2]
    if(k <= 0 || lam <= 0) return(1e10)
    -sum(dweibull(x, shape=k, scale=lam, log=TRUE))
  }

  fit <- nlm(nll,
             p=c(1.5, mean(x)),
             hessian=TRUE)

  k_hat  <- fit$estimate[1]
  lam_hat <- fit$estimate[2]

  V <- solve(fit$hessian)

  list(shape=k_hat,
       scale=lam_hat,
       vcov=V)
}

###############################################################
# 2. Point Estimation
###############################################################

LPIWeibull_point_est <- function(x, L){

  out <- LPIWeibull_mle_helper(x)

  k   <- out$shape
  lam <- out$scale

  mu    <- lam * gamma(1 + 1/k)
  sigma <- lam * sqrt(gamma(1 + 2/k) -
                        gamma(1 + 1/k)^2)

  CL <- (mu - L)/sigma

  list(Value=CL,
       mu=mu,
       sigma=sigma,
       shape=k,
       scale=lam,
       vcov=out$vcov)
}

###############################################################
# 3. Asymptotic CI (Delta Method)
###############################################################

LPIWeibull_asymptotic <- function(x, L, alpha=0.05){

  out <- LPIWeibull_point_est(x, L)

  k   <- out$shape
  lam <- out$scale
  mu  <- out$mu
  sigma <- out$sigma
  V   <- out$vcov

  IndexVal <- out$Value

  g1 <- gamma(1 + 1/k)
  g2 <- gamma(1 + 2/k)

  dmu_dlam <- g1
  dmu_dk   <- lam * g1 *
    digamma(1 + 1/k) * (-1/k^2)

  A <- g2 - g1^2

  dA_dk <- g2 * digamma(1+2/k)*(-2/k^2) -
    2*g1*digamma(1+1/k)*(-1/k^2)

  dsig_dlam <- sqrt(A)
  dsig_dk   <- lam*(dA_dk/(2*sqrt(A)))

  dCL_dmu  <- 1/sigma
  dCL_dsig <- -(mu-L)/sigma^2

  dCL_dlam <- dCL_dmu*dmu_dlam +
    dCL_dsig*dsig_dlam

  dCL_dk   <- dCL_dmu*dmu_dk +
    dCL_dsig*dsig_dk

  grad <- c(dCL_dk, dCL_dlam)

  varCL <- as.numeric(t(grad) %*% V %*% grad)
  seCL  <- sqrt(varCL)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seCL,
          IndexVal + qnorm(1-alpha/2)*seCL)

  list(Value=IndexVal,
       SE=seCL,
       LCL=CI[1],
       UCL=CI[2])
}

###############################################################
# 4. Bootstrap
###############################################################

LPIWeibull_bootstrap <- function(x, L, R=2000, alpha=0.05){

  boot_stat <- function(data, idx){
    xx <- data[idx]
    LPIWeibull_point_est(xx, L)$Value
  }

  b <- boot::boot(data=x,
                  statistic=boot_stat,
                  R=R)

  ci <- tryCatch(
    boot::boot.ci(b,
                  conf=1-alpha,
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

LPIWeibull_nonparam <- function(x, L, alpha=0.05){

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

  list(Value=CL,
       SE=seCL,
       LCL=CI[1],
       UCL=CI[2])
}

###############################################################
# 6. MAIN FUNCTION
###############################################################

#' LPI Weibull Classic CL Index
#'
#' Computes the classic capability index under the Weibull model.
#'
#' @param x numeric vector
#' @param L lower bound
#' @param R number of bootstrap replications
#' @param alpha significance level
#' @return data frame of results
#' @export
LPIWeibull <- function(x, L, R=1000, alpha=0.05){

  mle    <- LPIWeibull_point_est(x, L)
  asym   <- LPIWeibull_asymptotic(x, L, alpha)
  bootci <- LPIWeibull_bootstrap(x, L, R, alpha)
  nonpar <- LPIWeibull_nonparam(x, L, alpha)

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
      nonpar$Value,
      nonpar$LCL,
      nonpar$UCL
    )
  )

  invisible(TAB)
}
