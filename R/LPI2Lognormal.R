##############################################################
# LPI2LogNormal PACKAGE CODE
# Formula: Index = (Q_median - L) / (Q3 - Q1)
##############################################################

library(boot)
library(DT)
library(htmltools)

###############################################################
# 1. Internal: Lognormal MLE Helper
###############################################################
#' LPI2Lognormal
#' @export
#' @param x random sample
#' @param L lower bound
#' @param p1 first percentile
#' @param p3 second percentile
#' @param R the number of boostrap sample
#' @param alpha confidence level
#' LPI2LogNormal(x, L, p1, p3, R, alpha)




#'
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @return the output from list
LPI2LogNormal_mle_helper <- function(x){

  lx <- log(x)
  mu_hat <- mean(lx)
  sigma_hat <- sqrt(mean((lx - mu_hat)^2))

  n <- length(x)
  V <- matrix(c(
    sigma_hat^2/n, 0,
    0,             sigma_hat^2/(2*n)
  ), 2, 2, byrow=TRUE)

  list(mu=mu_hat, sigma=sigma_hat, vcov=V)
}


###############################################################
# 2. Internal: Point Estimation (MLE)
###############################################################

LPI2LogNormal_point_est <- function(x, p1=0.25, p3=0.75, L){

  out <- LPI2LogNormal_mle_helper(x)
  mu <- out$mu; s <- out$sigma

  Q1   <- exp(mu + s*qnorm(p1))
  Qmed <- exp(mu + s*qnorm(0.5))
  Q3   <- exp(mu + s*qnorm(p3))

  Val <- (Qmed - L)/(Q3 - Q1)

  list(Value=Val, mu=mu, sigma=s, vcov=out$vcov)
}

###############################################################
# 3. Internal: Asymptotic CI (Delta Method)
###############################################################

LPI2LogNormal_asymptotic <- function(x, p1=0.25, p3=0.75, L, alpha=0.05){

  out <- LPI2LogNormal_point_est(x, p1, p3, L)
  mu <- out$mu; s <- out$sigma; V <- out$vcov
  IndexVal <- out$Value

  z1  <- qnorm(p1)
  zmed <- qnorm(0.5)
  z3  <- qnorm(p3)

  Q1   <- exp(mu + s*z1)
  Qmed <- exp(mu + s*zmed)
  Q3   <- exp(mu + s*z3)

  Range <- Q3 - Q1

  dQ_dmu    <- function(Q) Q
  dQ_dsigma <- function(Q, z) Q * z

  dQ1_dmu    <- dQ_dmu(Q1)
  dQmed_dmu  <- dQ_dmu(Qmed)
  dQ3_dmu    <- dQ_dmu(Q3)

  dQ1_ds     <- dQ_dsigma(Q1, z1)
  dQmed_ds   <- dQ_dsigma(Qmed, zmed)
  dQ3_ds     <- dQ_dsigma(Q3, z3)

  dI_dQmed <- 1/Range
  dI_dQ1   <- (Qmed - L)/Range^2
  dI_dQ3   <- -(Qmed - L)/Range^2

  dI_dmu <- dI_dQ1*dQ1_dmu + dI_dQmed*dQmed_dmu + dI_dQ3*dQ3_dmu
  dI_ds  <- dI_dQ1*dQ1_ds  + dI_dQmed*dQmed_ds  + dI_dQ3*dQ3_ds

  grad <- c(dI_dmu, dI_ds)

  varI <- t(grad) %*% V %*% grad
  seI  <- sqrt(varI)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}


###############################################################
# 4. Internal: Bootstrap CI
###############################################################

LPI2LogNormal_boot_worker <- function(data, indices, p1, p3, L){
  x <- data[indices]
  LPI2LogNormal_point_est(x, p1, p3, L)$Value
}

LPI2LogNormal_bootstrap <- function(x, R=2000, p1=0.25, p3=0.75, L, alpha=0.05){

  boot_stat <- function(data, indices){
    LPI2LogNormal_boot_worker(data, indices, p1, p3, L)
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
# 5. Internal: Nonparametric CI
###############################################################

LPI2LogNormal_nonparam <- function(x, p1=0.25, p3=0.75, L, alpha=0.05){

  n <- length(x)

  q1   <- quantile(x, p1, type=8)
  qmed <- quantile(x, 0.5, type=8)
  q3   <- quantile(x, p3, type=8)

  Range <- q3 - q1
  IndexVal <- (qmed - L)/Range

  f1   <- density(x, from=q1, to=q1)$y[1]
  fmed <- density(x, from=qmed, to=qmed)$y[1]
  f3   <- density(x, from=q3, to=q3)$y[1]

  v11 <- p1*(1-p1)/(f1^2)
  v22 <- 0.5*(1-0.5)/(fmed^2)
  v33 <- p3*(1-p3)/(f3^2)

  v12 <- p1*(1-0.5)/(f1*fmed)
  v13 <- p1*(1-p3)/(f1*f3)
  v23 <- 0.5*(1-p3)/(fmed*f3)

  Sigma <- matrix(c(v11, v12, v13,
                    v12, v22, v23,
                    v13, v23, v33), 3, 3)

  dI_dQmed <- 1/Range
  dI_dQ1   <- (qmed - L)/Range^2
  dI_dQ3   <- -(qmed - L)/Range^2

  grad <- c(dI_dQ1, dI_dQmed, dI_dQ3)

  varI <- t(grad) %*% Sigma %*% grad / n
  seI  <- sqrt(varI)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)

  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}


###############################################################
# 6. OUTPUT (Print & Viewer)
###############################################################

LPI2LogNormal_print <- function(summary_table, alpha, x = NULL){


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
  cat("         LPI2LogNormal RESULTS\n")
  cat(sprintf("         Confidence Level: %.0f%%\n", (1-alpha)*100))
  cat("         Dist: Log-Normal (mu, sigma)\n")
  cat("         Formula: (Qmed - L) / (Q3 - Q1)\n")
  cat("=====================================\n\n")

  if(!is.null(x)){
    mle_par <- LPI2LogNormal_mle_helper(x)
    cat(sprintf("   mu_hat       : %.4f\n", mle_par$mu))
    cat(sprintf("   sigma_hat    : %.4f\n\n", mle_par$sigma))
  }

  cat(sprintf("MLE Value                : %.4f\n\n", mleV))
  cat(sprintf("Asymptotic CI            : ( %.4f , %.4f )\n\n", asymL, asymU))
  cat("Bootstrap CIs:\n")
  cat(sprintf("  Normal CI              : ( %.4f , %.4f )\n", normL, normU))
  cat(sprintf("  Basic CI               : ( %.4f , %.4f )\n", basicL, basicU))
  cat(sprintf("  Percentile CI          : ( %.4f , %.4f )\n", percL, percU))
  cat(sprintf("  BCa CI                 : ( %.4f , %.4f )\n\n", bcaL, bcaU))
  cat(sprintf("Nonparametric Value      : %.4f\n", npV))
  cat(sprintf("Nonparametric CI         : ( %.4f , %.4f )\n", npL, npU))
  cat("=====================================\n")
}

LPI2LogNormal_viewer <- function(summary_table, x, alpha, explanation="Log-Normal Distribution Analysis") {

  mle_par <- LPI2LogNormal_mle_helper(x)
  conf_pct <- (1 - alpha) * 100


  df <- data.frame(
    Item = c(
      "Location (mu_hat)",
      "Scale (sigma_hat)",
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
      sprintf("%.4f", mle_par$mu),
      sprintf("%.4f", mle_par$sigma),
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
    tags$h3("LPI2LogNormal Summary", style="margin-bottom:5px;"),
    tags$p(style="font-size:16px; color:#2c5aa0; font-weight:bold; margin-top:0px;",
           sprintf("LogNormal Distribution Analysis (Confidence Level: %g%%)", conf_pct)),
    tags$hr()
  )

  css_fix <- tags$style(HTML("
    .dataTables_wrapper { height: auto !important; overflow-y: visible !important; }
    table.dataTable { width: 100% !important; }
  "))

  return(browsable(tagList(
    css_fix,
    explanation_html,
    DT::datatable(
      df,
      rownames = FALSE,
      options = list(
        pageLength = nrow(df),
        scrollY = FALSE,
        paging = FALSE,
        dom = 't'
      )
    )
  )))
}
###############################################################
# 7. MAIN FUNCTION
###############################################################

LPI2LogNormal <- function(x, L, p1=0.25, p3=0.75, R=1000, alpha=0.05){

  mle     <- LPI2LogNormal_point_est(x, p1, p3, L)
  asym    <- LPI2LogNormal_asymptotic(x, p1, p3, L, alpha=alpha)
  bootci  <- LPI2LogNormal_bootstrap(x, R=R, p1, p3, L, alpha=alpha)
  nonpar  <- LPI2LogNormal_nonparam(x, p1, p3, L, alpha=alpha)

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

  LPI2LogNormal_print(TAB, alpha, x)
  print(LPI2LogNormal_viewer(TAB, x, alpha))

  invisible(TAB)
}

