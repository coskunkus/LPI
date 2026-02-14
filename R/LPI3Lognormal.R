##############################################################
# LPI3LogNormal
# Distribution: Log-Normal (meanlog, sdlog)
# Formula: (Median - L) / MAD
##############################################################

###############################################################
# 1. Internal: MLE Helper
###############################################################

LPI3LogNormal_mle_helper <- function(x){

  nll <- function(par){
    ml <- par[1]
    sl <- par[2]
    if(sl <= 0) return(1e10)
    -sum(dlnorm(x, meanlog=ml, sdlog=sl, log=TRUE))
  }

  start_mean <- mean(log(x))
  start_sd   <- sd(log(x))

  fit <- nlm(nll, p=c(start_mean, start_sd), hessian=TRUE)

  list(
    meanlog = fit$estimate[1],
    sdlog   = fit$estimate[2],
    vcov    = solve(fit$hessian)
  )
}

###############################################################
# 2. Point Estimation
###############################################################

LPI3LogNormal_point_est <- function(x, L){

  out <- LPI3LogNormal_mle_helper(x)
  ml  <- out$meanlog
  sl  <- out$sdlog

  Qmed <- exp(ml)

  mad_func <- function(m_val){
    lower <- Qmed - m_val
    if(lower < 0) lower <- 0

    plnorm(Qmed + m_val, meanlog=ml, sdlog=sl) -
      plnorm(lower, meanlog=ml, sdlog=sl) - 0.5
  }

  root_try <- try(uniroot(mad_func, c(0, Qmed*10)), silent=TRUE)

  if(inherits(root_try, "try-error")){
    return(list(Value=NA, meanlog=ml, sdlog=sl, vcov=out$vcov))
  }

  MAD_val <- root_try$root

  Val <- (Qmed - L) / MAD_val

  list(Value=Val, meanlog=ml, sdlog=sl, vcov=out$vcov)
}

###############################################################
# 3. Asymptotic CI
###############################################################

LPI3LogNormal_asymptotic <- function(x, L, alpha=0.05){

  out <- LPI3LogNormal_point_est(x, L)
  if(is.na(out$Value))
    return(list(Value=NA, SE=NA, LCL=NA, UCL=NA))

  ml_hat <- out$meanlog
  sl_hat <- out$sdlog
  V      <- out$vcov
  IndexVal <- out$Value

  index_fun <- function(par){

    ml <- par[1]
    sl <- par[2]
    if(sl <= 0) return(NA)

    Qmed <- exp(ml)

    mad_f <- function(m){
      lower <- Qmed - m
      if(lower < 0) lower <- 0
      plnorm(Qmed + m, meanlog=ml, sdlog=sl) -
        plnorm(lower, meanlog=ml, sdlog=sl) - 0.5
    }

    root_try <- try(uniroot(mad_f, c(0, Qmed*10)), silent=TRUE)
    if(inherits(root_try, "try-error")) return(NA)

    (Qmed - L) / root_try$root
  }

  eps <- 1e-5

  d_ml <- (index_fun(c(ml_hat+eps, sl_hat)) -
             index_fun(c(ml_hat-eps, sl_hat))) / (2*eps)

  d_sl <- (index_fun(c(ml_hat, sl_hat+eps)) -
             index_fun(c(ml_hat, sl_hat-eps))) / (2*eps)

  grad <- c(d_ml, d_sl)

  varI <- t(grad) %*% V %*% grad
  seI  <- sqrt(varI)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 4. Bootstrap
###############################################################

LPI3LogNormal_bootstrap <- function(x, R=2000, L, alpha=0.05){

  boot_stat <- function(data, indices){
    xx <- data[indices]
    LPI3LogNormal_point_est(xx, L)$Value
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

LPI3LogNormal_nonparam <- function(x, L, alpha=0.05){

  samp_med <- median(x)
  samp_mad <- mad(x, center=samp_med, constant=1)

  IndexVal <- (samp_med - L) / samp_mad

  boots <- replicate(1000, {
    xx <- sample(x, replace=TRUE)
    m  <- median(xx)
    md <- mad(xx, center=m, constant=1)
    if(md==0) return(NA)
    (m - L)/md
  })

  seI <- sd(boots, na.rm=TRUE)
  CI  <- quantile(boots,
                  probs=c(alpha/2, 1-alpha/2),
                  na.rm=TRUE)

  list(Value=IndexVal, SE=seI,
       LCL=CI[1], UCL=CI[2])
}

###############################################################
# 6. MAIN FUNCTION
###############################################################

#' LPI3 Log-Normal Index
#'
#' Computes the LPI3 index under the Log-Normal distribution.
#'
#' @param x numeric vector
#' @param L lower bound
#' @param R number of bootstrap replications
#' @param alpha significance level
#' @return data frame of results
#' @export
LPI3LogNormal <- function(x, L, R=1000, alpha=0.05){

  mle    <- LPI3LogNormal_point_est(x, L)
  asym   <- LPI3LogNormal_asymptotic(x, L, alpha)
  bootci <- LPI3LogNormal_bootstrap(x, R, L, alpha)
  nonpar <- LPI3LogNormal_nonparam(x, L, alpha)

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
