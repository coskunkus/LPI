##############################################################
# LPI3Gamma PACKAGE CODE (Application Mode)
# Formula: Index = (Median - L) / MAD
# Distribution: Gamma (shape, scale)
###############################################################

library(boot)
library(DT)
library(htmltools)

###############################################################
# 1. Internal: Gamma MLE Helper
###############################################################
#'LPI3Gamma
#' @export
#' @param x random sample
#' @param L lower bound
#' @param R the number of boostrap sample
#' @param alpha confidence level
#'  LPI3Gamma(x, L, R, alpha)

LPI3Gamma_mle_helper <- function(x){

  # Gamma Log-Likelihood
  nll <- function(par){
    k <- par[1]; lam <- par[2] # shape, scale
    if(k <= 0 || lam <= 0) return(1e10)
    -sum(dgamma(x, shape=k, scale=lam, log=TRUE))
  }

  # Ba??lang???? de??erleri (Method of Moments)
  mean_x <- mean(x)
  var_x  <- var(x)
  start_scale <- var_x / mean_x
  start_shape <- mean_x^2 / var_x

  fit <- nlm(nll, p=c(start_shape, start_scale), hessian=TRUE)
  k_hat <- fit$estimate[1]
  lam_hat <- fit$estimate[2]
  vcov_mat <- solve(fit$hessian)

  list(shape=k_hat, scale=lam_hat, vcov=vcov_mat)
}

###############################################################
# 2. Internal: Point Estimation (MLE)
###############################################################

LPI3Gamma_point_est <- function(x, L){

  out <- LPI3Gamma_mle_helper(x)
  k <- out$shape; lam <- out$scale

  # 1. Medyan (Q2)
  Qmed <- qgamma(0.5, shape=k, scale=lam)

  # 2. Teorik MAD (Gamma i??in N??merik ????z??m)
  mad_func <- function(m_val, k, lam, med_val){
    lower <- med_val - m_val
    if(lower < 0) lower <- 0

    p_upper <- pgamma(med_val + m_val, shape=k, scale=lam)
    p_lower <- pgamma(lower, shape=k, scale=lam)

    return((p_upper - p_lower) - 0.5)
  }

  # K??k bulma (Hata korumal??)
  try_root <- try(uniroot(mad_func, interval=c(0, Qmed*10), k=k, lam=lam, med_val=Qmed), silent=TRUE)

  if(inherits(try_root, "try-error")){
    return(list(Value=NA, shape=k, scale=lam, vcov=out$vcov))
  }

  MAD_theoretical <- try_root$root

  # FORM??L: (Median - L) / MAD
  Val <- (Qmed - L) / MAD_theoretical

  list(Value=Val, shape=k, scale=lam, vcov=out$vcov)
}

###############################################################
# 3. Internal: Asymptotic CI (Numerical Gradient)
###############################################################

LPI3Gamma_asymptotic <- function(x, L, alpha=0.05){

  out <- LPI3Gamma_point_est(x, L)
  if(is.na(out$Value)) return(list(Value=NA, SE=NA, LCL=NA, UCL=NA))

  k_hat <- out$shape; lam_hat <- out$scale; V <- out$vcov
  IndexVal <- out$Value

  # Endeks Fonksiyonu (T??rev i??in)
  index_fun <- function(pars){
    k <- pars[1]; lam <- pars[2]
    if(k<=0 || lam<=0) return(NA)

    qm <- qgamma(0.5, shape=k, scale=lam)

    mad_f <- function(m, k, lam, med){
      low <- med - m; if(low<0) low <- 0
      pgamma(med+m, shape=k, scale=lam) - pgamma(low, shape=k, scale=lam) - 0.5
    }

    res <- try(uniroot(mad_f, c(0, qm*10), k=k, lam=lam, med=qm)$root, silent=TRUE)
    if(inherits(res, "try-error")) return(NA)

    return( (qm - L) / res )
  }

  # N??merik T??rev
  my_grad <- function(func, p){
    eps <- 1e-5
    v1 <- func(c(p[1]+eps, p[2])); v2 <- func(c(p[1]-eps, p[2]))
    dk <- (v1-v2)/(2*eps)
    v3 <- func(c(p[1], p[2]+eps)); v4 <- func(c(p[1], p[2]-eps))
    dlam <- (v3-v4)/(2*eps)
    return(c(dk, dlam))
  }

  grad_vec <- my_grad(index_fun, c(k_hat, lam_hat))

  varI <- t(grad_vec) %*% V %*% grad_vec
  seI <- sqrt(varI)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 4. Internal: Bootstrap CI (Robust BCa)
###############################################################

LPI3Gamma_boot_worker <- function(data, indices, L){
  x <- data[indices]
  val <- try(LPI3Gamma_point_est(x, L)$Value, silent=TRUE)
  if(inherits(val, "try-error")) return(NA)
  return(val)
}

LPI3Gamma_bootstrap <- function(x, R=2000, L, alpha=0.05){

  boot_stat <- function(data, indices){
    LPI3Gamma_boot_worker(data, indices, L)
  }

  b <- boot(data=x, statistic=boot_stat, R=R)

  # BCa Hatas??na Kar???? Koruma (Robust) ve Alpha Eklendi
  ci <- tryCatch({
    boot.ci(b, conf = 1-alpha, type=c("norm","basic","perc","bca"))
  }, error = function(e) {
    warning("BCa aral?????? hesaplanamad??. Di??er aral??klar sunuluyor.")
    tryCatch({
      boot.ci(b, conf = 1-alpha, type=c("norm","basic","perc"))
    }, error = function(e2) NULL)
  })

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
# 5. Internal: Nonparametric Method (Sample MAD)
###############################################################

LPI3Gamma_nonparam <- function(x, L, alpha=0.05){
  samp_med <- median(x)
  samp_mad <- mad(x, center = samp_med, constant = 1)

  IndexVal <- (samp_med - L) / samp_mad

  # Non-param Bootstrap
  boots <- replicate(1000, {
    xx <- sample(x, replace=TRUE)
    m <- median(xx)
    md <- mad(xx, center=m, constant=1)
    if(md==0) return(NA)
    (m - L)/md
  })

  seI <- sd(boots, na.rm=TRUE)
  CI <- quantile(boots, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 6. OUTPUT HELPER FUNCTIONS
###############################################################

LPI3Gamma_print <- function(summary_table, x=NULL, alpha=0.05){

  conf_pct <- (1-alpha)*100

  mleV   <- summary_table$Value[summary_table$Item=="MLE_Value"]
  asymL  <- summary_table$Value[summary_table$Item=="Asymp_LCL"]
  asymU  <- summary_table$Value[summary_table$Item=="Asymp_UCL"]

  normL  <- summary_table$Value[summary_table$Item=="Boot_norm_LCL"]
  normU  <- summary_table$Value[summary_table$Item=="Boot_norm_UCL"]
  basicL <- summary_table$Value[summary_table$Item=="Boot_basic_LCL"]
  basicU <- summary_table$Value[summary_table$Item=="Boot_basic_UCL"]
  percL  <- summary_table$Value[summary_table$Item=="Boot_perc_LCL"]
  percU  <- summary_table$Value[summary_table$Item=="Boot_perc_UCL"]
  bcaL   <- summary_table$Value[summary_table$Item=="Boot_bca_LCL"]
  bcaU   <- summary_table$Value[summary_table$Item=="Boot_bca_UCL"]

  npV    <- summary_table$Value[summary_table$Item=="Nonpar_Value"]
  npL    <- summary_table$Value[summary_table$Item=="Nonpar_LCL"]
  npU    <- summary_table$Value[summary_table$Item=="Nonpar_UCL"]

  cat("=====================================\n")
  cat("     LPI3Gamma RESULTS\n")
  cat(sprintf("   Confidence Level: %.0f%%\n", conf_pct))
  cat("     Formula: (Median - L) / MAD\n")
  cat("=====================================\n\n")

  if(!is.null(x)){
    mle_par <- LPI3Gamma_mle_helper(x)
    cat(sprintf("   Shape (k_hat)  : %.4f\n", mle_par$shape))
    cat(sprintf("   Scale (lam_hat): %.4f\n\n", mle_par$scale))
  }

  cat(sprintf("MLE Value              : %.4f\n\n", mleV))
  cat(sprintf("Asymptotic CI          : ( %.4f , %.4f )\n", asymL, asymU))
  cat(sprintf("Bootstrap Normal       : ( %.4f , %.4f )\n", normL, normU))
  cat(sprintf("Bootstrap Basic        : ( %.4f , %.4f )\n", basicL, basicU))
  cat(sprintf("Bootstrap Percent      : ( %.4f , %.4f )\n", percL, percU))
  cat(sprintf("Bootstrap BCa          : ( %.4f , %.4f )\n", bcaL, bcaU))

  cat(sprintf("\nNonparametric Value    : %.4f\n", npV))
  cat(sprintf("Nonparametric CI       : ( %.4f , %.4f )\n", npL, npU))
  cat("\n=====================================\n")
}

LPI3Gamma_viewer <- function(summary_table, x, alpha, explanation="LPI3Gamma Analysis") {

  mle_par <- LPI3Gamma_mle_helper(x)
  conf_level <- (1 - alpha) * 100

  items <- c("Shape (k_hat)", "Scale (lam_hat)",
             "MLE Value", "Asymptotic CI",
             "Bootstrap Normal", "Bootstrap Basic", "Bootstrap Percent", "Bootstrap BCa",
             "Nonparametric Value", "Nonparametric CI")

  vals <- c(
    sprintf("%.4f", mle_par$shape), sprintf("%.4f", mle_par$scale),
    sprintf("%.4f", summary_table$Value[summary_table$Item=="MLE_Value"]),
    sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Asymp_LCL"], summary_table$Value[summary_table$Item=="Asymp_UCL"]),
    sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Boot_norm_LCL"], summary_table$Value[summary_table$Item=="Boot_norm_UCL"]),
    sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Boot_basic_LCL"], summary_table$Value[summary_table$Item=="Boot_basic_UCL"]),
    sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Boot_perc_LCL"], summary_table$Value[summary_table$Item=="Boot_perc_UCL"]),
    sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Boot_bca_LCL"], summary_table$Value[summary_table$Item=="Boot_bca_UCL"]),
    sprintf("%.4f", summary_table$Value[summary_table$Item=="Nonpar_Value"]),
    sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Nonpar_LCL"], summary_table$Value[summary_table$Item=="Nonpar_UCL"])
  )

  df <- data.frame(Item = items, Value = vals, stringsAsFactors = FALSE)

  info_text <- paste0(explanation, " (Confidence Level: ", conf_level, "%)")

  explanation_html <- tags$div(
    tags$h3("LPI3Gamma Summary", style="margin-bottom:5px;"),
    tags$p(style="font-size:16px; color:#0056b3; margin-top:0px; font-weight:bold;", info_text),
    tags$hr()
  )

  css_fix <- tags$style(HTML(".dataTables_wrapper { height: auto !important; overflow-y: hidden !important; } table.dataTable { width: 100% !important; }"))

  return(browsable(tagList(css_fix, explanation_html, DT::datatable(df, rownames = FALSE, options = list(pageLength = nrow(df), scrollY = FALSE, paging = FALSE, dom = 't')))))
}

###############################################################
# 7. MAIN FUNCTION: LPI3Gamma (Application Mode)
###############################################################

LPI3Gamma <- function(x, L, R=1000, alpha=0.05){


  mle      <- LPI3Gamma_point_est(x, L)
  asym     <- LPI3Gamma_asymptotic(x, L, alpha=alpha)
  bootci   <- LPI3Gamma_bootstrap(x, R=R, L, alpha=alpha)
  nonpar   <- LPI3Gamma_nonparam(x, L, alpha=alpha)

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

  # Output
  LPI3Gamma_print(TAB, x, alpha)
  print(LPI3Gamma_viewer(TAB, x, alpha))

  invisible(TAB)
}
