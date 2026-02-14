##############################################################
# LPI3LogNormal PACKAGE CODE (Application Mode)
# Formula: Index = (Median - L) / MAD
# Distribution: Log-Normal (meanlog, sdlog)
###############################################################

library(boot)
library(DT)
library(htmltools)

###############################################################
# 1. Internal: Log-Normal MLE Helper
###############################################################
#' LPI3Lognormal
#' @export
#' @param x random sample
#' @param L lower bound
#' @param R the number of boostrap sample
#' @param alpha confidence level
#'  LPI3LogNormal(x, L, R, alpha)
LPI3LogNormal_mle_helper <- function(x){

  # Log-Normal Log-Likelihood
  nll <- function(par){
    ml <- par[1]; sl <- par[2] # meanlog, sdlog
    if(sl <= 0) return(1e10)
    -sum(dlnorm(x, meanlog=ml, sdlog=sl, log=TRUE))
  }

  # Başlangıç değerleri
  start_mean <- mean(log(x))
  start_sd   <- sd(log(x))

  fit <- nlm(nll, p=c(start_mean, start_sd), hessian=TRUE)
  ml_hat <- fit$estimate[1]
  sl_hat <- fit$estimate[2]
  vcov_mat <- solve(fit$hessian)

  list(meanlog=ml_hat, sdlog=sl_hat, vcov=vcov_mat)
}

###############################################################
# 2. Internal: Point Estimation (MLE)
###############################################################

LPI3LogNormal_point_est <- function(x, L){

  out <- LPI3LogNormal_mle_helper(x)
  ml <- out$meanlog; sl <- out$sdlog

  # 1. Medyan (Q2) = exp(meanlog)
  Qmed <- exp(ml)

  # 2. Teorik MAD (Log-Normal için Numerik Çözüm)
  mad_func <- function(m_val, ml, sl, med_val){
    lower <- med_val - m_val
    if(lower < 0) lower <- 0

    p_upper <- plnorm(med_val + m_val, meanlog=ml, sdlog=sl)
    p_lower <- plnorm(lower, meanlog=ml, sdlog=sl)

    return((p_upper - p_lower) - 0.5)
  }

  try_root <- try(uniroot(mad_func, interval=c(0, Qmed*10), ml=ml, sl=sl, med_val=Qmed), silent=TRUE)

  if(inherits(try_root, "try-error")){
    return(list(Value=NA, meanlog=ml, sdlog=sl, vcov=out$vcov))
  }

  MAD_theoretical <- try_root$root

  Val <- (Qmed - L) / MAD_theoretical

  list(Value=Val, meanlog=ml, sdlog=sl, vcov=out$vcov)
}

###############################################################
# 3. Internal: Asymptotic CI (Numerical Gradient)
###############################################################

LPI3LogNormal_asymptotic <- function(x, L, alpha=0.05){

  out <- LPI3LogNormal_point_est(x, L)
  if(is.na(out$Value)) return(list(Value=NA, SE=NA, LCL=NA, UCL=NA))

  ml_hat <- out$meanlog; sl_hat <- out$sdlog; V <- out$vcov
  IndexVal <- out$Value

  index_fun <- function(pars){
    ml <- pars[1]; sl <- pars[2]
    if(sl<=0) return(NA)

    qm <- exp(ml)

    mad_f <- function(m, ml, sl, med){
      low <- med - m; if(low<0) low <- 0
      plnorm(med+m, ml, sl) - plnorm(low, ml, sl) - 0.5
    }

    res <- try(uniroot(mad_f, c(0, qm*10), ml=ml, sl=sl, med=qm)$root, silent=TRUE)
    if(inherits(res, "try-error")) return(NA)

    return( (qm - L) / res )
  }

  my_grad <- function(func, p){
    eps <- 1e-5
    v1 <- func(c(p[1]+eps, p[2])); v2 <- func(c(p[1]-eps, p[2]))
    dml <- (v1-v2)/(2*eps)
    v3 <- func(c(p[1], p[2]+eps)); v4 <- func(c(p[1], p[2]-eps))
    dsl <- (v3-v4)/(2*eps)
    return(c(dml, dsl))
  }

  grad_vec <- my_grad(index_fun, c(ml_hat, sl_hat))

  varI <- t(grad_vec) %*% V %*% grad_vec
  seI <- sqrt(varI)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 4. Internal: Bootstrap CI
###############################################################

LPI3LogNormal_boot_worker <- function(data, indices, L){
  x <- data[indices]
  val <- try(LPI3LogNormal_point_est(x, L)$Value, silent=TRUE)
  if(inherits(val, "try-error")) return(NA)
  return(val)
}

LPI3LogNormal_bootstrap <- function(x, R=2000, L, alpha=0.05){

  boot_stat <- function(data, indices){
    LPI3LogNormal_boot_worker(data, indices, L)
  }

  b <- boot(data=x, statistic=boot_stat, R=R)

  ci <- tryCatch({
    boot.ci(b, conf = (1 - alpha), type=c("norm","basic","perc","bca"))
  }, error = function(e) {
    warning("BCa aralığı hesaplanamadı. Diğer aralıklar sunuluyor.")
    tryCatch({
      boot.ci(b, conf = (1 - alpha), type=c("norm","basic","perc"))
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
# 5. Internal: Nonparametric Method
###############################################################

LPI3LogNormal_nonparam <- function(x, L, alpha=0.05){
  samp_med <- median(x)
  samp_mad <- mad(x, center = samp_med, constant = 1)

  IndexVal <- (samp_med - L) / samp_mad

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

LPI3LogNormal_print <- function(summary_table, alpha, x = NULL){

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
  cat("     LPI3LogNormal RESULTS\n")
  cat(sprintf("     Confidence Level: %.0f%%\n", (1-alpha)*100))
  cat("     Formula: (Median - L) / MAD\n")
  cat("=====================================\n\n")

  if(!is.null(x)){
    mle_par <- LPI3LogNormal_mle_helper(x)
    cat(sprintf("   meanlog (mu_hat): %.4f\n", mle_par$meanlog))
    cat(sprintf("   sdlog (sig_hat) : %.4f\n\n", mle_par$sdlog))
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

LPI3LogNormal_viewer <- function(summary_table, x, alpha, explanation="Log-Normal Distribution Analysis") {

  mle_par <- LPI3LogNormal_mle_helper(x)
  conf_pct <- (1 - alpha) * 100

  items <- c("Location (mu_hat)", "Scale (sigma_hat)",
             "MLE Value", "Asymptotic CI",
             "Bootstrap Normal", "Bootstrap Basic", "Bootstrap Percent", "Bootstrap BCa",
             "Nonparametric Value", "Nonparametric CI")

  vals <- c(
    sprintf("%.4f", mle_par$meanlog), sprintf("%.4f", mle_par$sdlog),
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

  explanation_html <- tags$div(
    tags$h3("LPI3LogNormal Summary", style="margin-bottom:5px; font-weight:bold;"),
    tags$p(style="font-size:16px; color:#2c5aa0; font-weight:bold; margin-top:0px;",
           sprintf("LogNormal Distribution Analysis (Confidence Level: %g%%)", conf_pct)),
    tags$hr()
  )

  css_fix <- tags$style(HTML(".dataTables_wrapper { height: auto !important; overflow-y: hidden !important; } table.dataTable { width: 100% !important; }"))

  return(browsable(tagList(
    css_fix,
    explanation_html,
    DT::datatable(df, rownames = FALSE,
                  options = list(
                    pageLength = nrow(df),
                    scrollY = FALSE,
                    paging = FALSE,
                    dom = 't'
                  ))
  )))
}

###############################################################
# 7. MAIN FUNCTION: LPI3LogNormal (Application Mode)
###############################################################

LPI3LogNormal <- function(x, L, R=1000, alpha=0.05){

  mle     <- LPI3LogNormal_point_est(x, L)
  asym    <- LPI3LogNormal_asymptotic(x, L, alpha=alpha)
  bootci  <- LPI3LogNormal_bootstrap(x, R=R, L, alpha=alpha)
  nonpar  <- LPI3LogNormal_nonparam(x, L, alpha=alpha)

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

  LPI3LogNormal_print(TAB, alpha, x)
  print(LPI3LogNormal_viewer(TAB, x, alpha))

  invisible(TAB)
}

