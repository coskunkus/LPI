##############################################################
# LPIlognormal
# Classic CL Index
# Formula: CL = (mu - L) / sigma  (Lognormal moments)
##############################################################

###############################################################
# 1. Lognormal MLE (Closed-form)
###############################################################

LPIlognormal_mle_helper <- function(x){

  lx <- log(x)

  mu_hat  <- mean(lx)
  sig_hat <- sd(lx)

  n <- length(x)

  V <- matrix(
    c(sig_hat^2/n,        0,
      0,                  sig_hat^2/(2*n)),
    2, 2, byrow=TRUE
  )

  list(mu=mu_hat,
       sigma=sig_hat,
       vcov=V)
}

###############################################################
# 2. Moments & CL
###############################################################

LPIlognormal_point_est <- function(x, L){

  out <- LPIlognormal_mle_helper(x)

  mu <- out$mu
  s  <- out$sigma

  muX  <- exp(mu + s^2/2)
  sigX <- sqrt((exp(s^2)-1) * exp(2*mu + s^2))

  CL <- (muX - L)/sigX

  list(Value=CL,
       mu=mu,
       sigma=s,
       muX=muX,
       sigX=sigX,
       vcov=out$vcov)
}

###############################################################
# 3. Asymptotic CI (Delta Method)
###############################################################

LPIlognormal_asymptotic <- function(x, L, alpha=0.05){

  out <- LPIlognormal_point_est(x, L)

  mu   <- out$mu
  s    <- out$sigma
  muX  <- out$muX
  sigX <- out$sigX
  V    <- out$vcov

  IndexVal <- out$Value

  dmuX_dmu <- muX
  dmuX_ds  <- muX * s

  dsigX_dmu <- sigX
  dsigX_ds  <- sigX * ( (s*exp(s^2))/(exp(s^2)-1) + s )

  dCL_dmuX  <- 1/sigX
  dCL_dsigX <- -(muX - L)/sigX^2

  dCL_dmu <- dCL_dmuX*dmuX_dmu +
    dCL_dsigX*dsigX_dmu

  dCL_ds  <- dCL_dmuX*dmuX_ds  +
    dCL_dsigX*dsigX_ds

  grad <- c(dCL_dmu, dCL_ds)

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

LPIlognormal_bootstrap <- function(x, L, R=2000, alpha=0.05){

  boot_stat <- function(data, indices){
    xx <- data[indices]
    LPIlognormal_point_est(xx, L)$Value
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

LPIlognormal_nonparam <- function(x, L, alpha=0.05){

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

#' LPI Lognormal Classic CL Index
#'
#' Computes the classic capability index under the Lognormal model.
#'
#' @param x numeric vector
#' @param L lower bound
#' @param R number of bootstrap replications
#' @param alpha significance level
#' @return data frame of results
#' @export
LPIlognormal <- function(x, L, R=1000, alpha=0.05){

  mle    <- LPIlognormal_point_est(x, L)
  asym   <- LPIlognormal_asymptotic(x, L, alpha)
  bootci <- LPIlognormal_bootstrap(x, L, R, alpha)
  nonpar <- LPIlognormal_nonparam(x, L, alpha)

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
