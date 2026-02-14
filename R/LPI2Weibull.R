##############################################################
# LPI2Weibull PACKAGE CODE (Auto-Output & No True Params)
# Formula: Index = (Q_median - L) / (Q3 - Q1)
###############################################################

library(boot)
library(DT)
library(htmltools)

###############################################################
# 1. Internal: Weibull MLE Helper
###############################################################

LPI2Weibull_mle_helper <- function(x){

  nll <- function(par){
    k <- par[1]; lam <- par[2]
    if(k <= 0 || lam <= 0) return(1e10)
    -sum(dweibull(x, shape=k, scale=lam, log=TRUE))
  }

  # Ba??lang???? de??erleri
  fit <- nlm(nll, p=c(1.2, mean(x)), hessian=TRUE)
  k_hat <- fit$estimate[1]
  lam_hat <- fit$estimate[2]
  vcov_mat <- solve(fit$hessian)

  list(shape=k_hat, scale=lam_hat, vcov=vcov_mat)
}

###############################################################
# 2. Internal: Point Estimation (MLE)
###############################################################

LPI2Weibull_point_est <- function(x, p1=0.25, p3=0.75, L){

  out <- LPI2Weibull_mle_helper(x)
  k <- out$shape; lam <- out$scale

  # Q1 (p1), Q2 (Median=0.5), Q3 (p3)
  Q1   <- lam * (-log(1-p1))^(1/k)
  Qmed <- lam * (-log(1-0.5))^(1/k)
  Q3   <- lam * (-log(1-p3))^(1/k)

  # FORM??L: (Q_median - L) / (Q3 - Q1)
  Val <- (Qmed - L)/(Q3 - Q1)

  list(Value=Val, shape=k, scale=lam, vcov=out$vcov)
}

###############################################################
# 3. Internal: Asymptotic CI (Delta Method)
###############################################################

LPI2Weibull_asymptotic <- function(x, p1=0.25, p3=0.75, L, alpha=0.05){

  out <- LPI2Weibull_point_est(x, p1, p3, L)
  k <- out$shape; lam <- out$scale; V <- out$vcov
  IndexVal <- out$Value

  # Quantile Hesaplar??
  p_med <- 0.5
  Q1   <- lam * (-log(1-p1))^(1/k)
  Qmed <- lam * (-log(1-p_med))^(1/k)
  Q3   <- lam * (-log(1-p3))^(1/k)

  Range <- Q3 - Q1

  # --- T??REVLER ---
  dQ_dlambda <- function(Q, lam) Q/lam
  dQ_dk      <- function(p, Q, k) Q * (-log(-log(1-p)))/k^2

  dQ1_dlam   <- dQ_dlambda(Q1, lam)
  dQmed_dlam <- dQ_dlambda(Qmed, lam)
  dQ3_dlam   <- dQ_dlambda(Q3, lam)

  dQ1_dk     <- dQ_dk(p1, Q1, k)
  dQmed_dk   <- dQ_dk(p_med, Qmed, k)
  dQ3_dk     <- dQ_dk(p3, Q3, k)

  # Endeksin Quantile'lara g??re t??revleri: I = (Qmed - L) / (Q3 - Q1)
  dI_dQmed <- 1 / Range
  dI_dQ1   <- (Qmed - L) / Range^2
  dI_dQ3   <- -(Qmed - L) / Range^2

  # Zincir Kural??
  dI_dk   <- dI_dQ1*dQ1_dk + dI_dQmed*dQmed_dk + dI_dQ3*dQ3_dk
  dI_dlam <- dI_dQ1*dQ1_dlam + dI_dQmed*dQmed_dlam + dI_dQ3*dQ3_dlam

  grad <- c(dI_dk, dI_dlam)

  varI <- t(grad) %*% V %*% grad
  seI <- sqrt(varI)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 4. Internal: Bootstrap CI
###############################################################
#' LPI2Weibull
#' @export
#' @param x random sample
#' @param L lower bound
#' @param p1 first percentile
#' @param p3 second percentile
#' @param R the number of boostrap sample
#' @param alpha confidence level
#' LPI2LogNormal(x, L, p1, p3, R, alpha)
LPI2Weibull_boot_worker <- function(data, indices, p1, p3, L){
  x <- data[indices]
  LPI2Weibull_point_est(x, p1, p3, L)$Value
}

LPI2Weibull_bootstrap <- function(x, R=2000, p1=0.25, p3=0.75, L, alpha=0.05){

  boot_stat <- function(data, indices){
    LPI2Weibull_boot_worker(data, indices, p1, p3, L)
  }

  b <- boot(data=x, statistic=boot_stat, R=R)
  ci <- boot.ci(b, conf = (1 - alpha), type=c("norm","basic","perc","bca"))

  list(
    Value = b$t0,
    norm = ci$normal[2:3],
    basic = ci$basic[4:5],
    perc  = ci$percent[4:5],
    bca   = ci$bca[4:5]
  )
}

###############################################################
# 5. Internal: Nonparametric Method
###############################################################

LPI2Weibull_nonparam <- function(x, p1=0.25, p3=0.75, L, alpha=0.05){

  n <- length(x)
  p_med <- 0.5

  # Sample Quantiles (Type 8)
  q1   <- quantile(x, p1, type=8)
  qmed <- quantile(x, p_med, type=8)
  q3   <- quantile(x, p3, type=8)

  Range <- q3 - q1
  IndexVal <- (qmed - L) / Range

  # Density estimation
  f1   <- density(x, from=q1, to=q1)$y[1]
  fmed <- density(x, from=qmed, to=qmed)$y[1]
  f3   <- density(x, from=q3, to=q3)$y[1]

  # Covariance terms
  v11 <- p1*(1-p1)/(f1^2)
  v22 <- p_med*(1-p_med)/(fmed^2)
  v33 <- p3*(1-p3)/(f3^2)

  v12 <- (p1*(1-p_med))/(f1*fmed)
  v13 <- (p1*(1-p3))/(f1*f3)
  v23 <- (p_med*(1-p3))/(fmed*f3)

  Sigma <- matrix(c(v11, v12, v13,
                    v12, v22, v23,
                    v13, v23, v33), 3, 3)

  # Gradients
  dI_dqmed <- 1 / Range
  dI_dq1   <- (qmed - L) / Range^2
  dI_dq3   <- -(qmed - L) / Range^2

  grad <- c(dI_dq1, dI_dqmed, dI_dq3)

  varI <- t(grad) %*% Sigma %*% grad / n
  seI  <- sqrt(varI)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 6. OUTPUT HELPER FUNCTIONS (Print & Viewer)
###############################################################

LPI2Weibull_print <- function(summary_table, x=NULL, alpha=0.05){

  # --- Extract All Values ---
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


  # --- HEADER ---
  cat("=====================================\n")
  cat("          LPI2Weibull RESULTS\n")
  cat(sprintf("     Confidence Level: %d%%\n", (1-alpha)*100))
  cat("      Dist: Weibull (Shape, Scale)\n")
  cat("      Formula: (Qmed - L) / (Q3 - Q1)\n")
  cat("=====================================\n\n")


  # --- MLE PARAMS ---
  if(!is.null(x)){
    mle_par <- LPI2Weibull_mle_helper(x)
    cat(sprintf("   Shape (k_hat)  : %.4f\n", mle_par$shape))
    cat(sprintf("   Scale (lam_hat): %.4f\n\n", mle_par$scale))
  }


  # --- MLE VALUE ---
  cat(sprintf("MLE Value                : %.4f\n\n", mleV))


  # --- ASYMPTOTIC CI ---
  cat(sprintf("Asymptotic CI            : ( %.4f , %.4f )\n\n", asymL, asymU))


  # --- BOOTSTRAP CIs (ALL FOUR TYPES) ---
  cat("Bootstrap CIs:\n")

  cat(sprintf("  Normal CI              : ( %.4f , %.4f )\n", normL, normU))
  cat(sprintf("  Basic CI               : ( %.4f , %.4f )\n", basicL, basicU))
  cat(sprintf("  Percentile CI          : ( %.4f , %.4f )\n", percL, percU))
  cat(sprintf("  BCa CI                 : ( %.4f , %.4f )\n\n", bcaL, bcaU))


  # --- NONPARAMETRIC ---
  cat(sprintf("Nonparametric Value      : %.4f\n", npV))
  cat(sprintf("Nonparametric CI         : ( %.4f , %.4f )\n", npL, npU))


  cat("=====================================\n")
}

LPI2Weibull_viewer <- function(summary_table, x, alpha=0.05, explanation="LPI2Weibull Analysis") {

  mle_par <- LPI2Weibull_mle_helper(x)
  conf_text <- sprintf("Confidence Level: %d%%", (1-alpha)*100)

  df <- data.frame(
    Item = c(
      "Shape (k_hat)",
      "Scale (lam_hat)",
      "MLE Value",
      "Asymptotic CI",
      "Bootstrap Normal",
      "Bootstrap Basic",
      "Bootstrap Percent",
      "Bootstrap BCa",
      "Nonparametric Value",
      "Nonparametric CI"
    ),
    Value = c(
      sprintf("%.4f", mle_par$shape),
      sprintf("%.4f", mle_par$scale),
      sprintf("%.4f", summary_table$Value[summary_table$Item=="MLE_Value"]),
      sprintf("( %.4f , %.4f )",
              summary_table$Value[summary_table$Item=="Asymp_LCL"],
              summary_table$Value[summary_table$Item=="Asymp_UCL"]),
      sprintf("( %.4f , %.4f )",
              summary_table$Value[summary_table$Item=="Boot_norm_LCL"],
              summary_table$Value[summary_table$Item=="Boot_norm_UCL"]),
      sprintf("( %.4f , %.4f )",
              summary_table$Value[summary_table$Item=="Boot_basic_LCL"],
              summary_table$Value[summary_table$Item=="Boot_basic_UCL"]),
      sprintf("( %.4f , %.4f )",
              summary_table$Value[summary_table$Item=="Boot_perc_LCL"],
              summary_table$Value[summary_table$Item=="Boot_perc_UCL"]),
      sprintf("( %.4f , %.4f )",
              summary_table$Value[summary_table$Item=="Boot_bca_LCL"],
              summary_table$Value[summary_table$Item=="Boot_bca_UCL"]),
      sprintf("%.4f",
              summary_table$Value[summary_table$Item=="Nonpar_Value"]),
      sprintf("( %.4f , %.4f )",
              summary_table$Value[summary_table$Item=="Nonpar_LCL"],
              summary_table$Value[summary_table$Item=="Nonpar_UCL"])
    ),
    stringsAsFactors = FALSE
  )

  explanation_html <- tags$div(
    tags$h3("LPI2Weibull Summary", style="margin-bottom:5px;"),
    tags$p(style="font-size:16px; color:#0056b3; font-weight:bold; margin-top:0px;",
           explanation, paste0("(", conf_text, ")")),
    tags$hr()
  )

  css_fix <- tags$style(HTML("
    .dataTables_wrapper { height: auto !important; overflow-y: hidden !important; }
    table.dataTable { width: 100% !important; }
  "))

  return(browsable(tagList(css_fix, explanation_html, DT::datatable(df, rownames = FALSE, options = list(pageLength = nrow(df), scrollY = FALSE, paging = FALSE, dom = 't')))))
}

###############################################################
# 7. MAIN FUNCTION: LPI2Weibull (Auto Executes Output)
###############################################################

LPI2Weibull <- function(x, L, p1=0.25, p3=0.75, R=1000, alpha=0.05){

  # 1. Calculate Results
  mle     <- LPI2Weibull_point_est(x, p1, p3, L)
  asym    <- LPI2Weibull_asymptotic(x, p1, p3, L, alpha=alpha)
  bootci  <- LPI2Weibull_bootstrap(x, R=R, p1, p3, L, alpha=alpha)
  nonpar  <- LPI2Weibull_nonparam(x, p1, p3, L, alpha=alpha)

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
      bootci$bca[1],   bootci$bca[2],
      nonpar$Value, nonpar$LCL, nonpar$UCL
    )
  )

  # 2. AUTO OUTPUT: Print to Console
  LPI2Weibull_print(TAB, x, alpha=alpha)

  # 3. AUTO OUTPUT: Show in Viewer
  print(LPI2Weibull_viewer(TAB, x, alpha=alpha))

  # 4. Return data silently
  invisible(TAB)
}

