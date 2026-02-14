##############################################################
# LPIlognormal — Classic CL Index with Full Bootstrap CI
# Formula:  CL = (mu - L) / sigma   (Lognormal moments)
##############################################################

library(boot)
library(DT)
library(htmltools)

###############################################################
# 1. Lognormal MLE (closed-form)
###############################################################
#'Lognormal
#' @export
#' @param x random sample
#' @param L lower bound
#' @param R the number of boostrap sample
#' @param alpha confidence level
#' LPIlognormal(x, L, R, alpha)
LPIlognormal_mle_helper <- function(x){
  lx <- log(x)
  mu_hat <- mean(lx)
  sig_hat <- sd(lx)

  n <- length(x)
  V <- matrix(c(
    sig_hat^2/n,       0,
    0,                 sig_hat^2/(2*n)
  ), 2, 2, byrow=TRUE)

  list(mu=mu_hat, sigma=sig_hat, vcov=V)
}

###############################################################
# 2. Moments & CL
###############################################################

LPIlognormal_point_est <- function(x, L){
  out <- LPIlognormal_mle_helper(x)
  mu <- out$mu
  s  <- out$sigma

  mu_X <- exp(mu + s^2/2)
  sig_X <- sqrt((exp(s^2)-1) * exp(2*mu + s^2))

  CL <- (mu_X - L)/sig_X

  list(Value=CL, mu=mu, sigma=s, muX=mu_X, sigX=sig_X, vcov=out$vcov)
}

###############################################################
# 3. Asymptotic CI (Delta Method)
###############################################################

LPIlognormal_asymptotic <- function(x, L, alpha=0.05){
  out <- LPIlognormal_point_est(x, L)
  mu <- out$mu; s <- out$sigma
  muX <- out$muX; sigX <- out$sigX
  V <- out$vcov
  IndexVal <- out$Value

  dmuX_dmu <- muX
  dmuX_ds  <- muX * s

  dsigX_dmu <- sigX * (1)
  dsigX_ds  <- sigX * ( (s*exp(s^2))/(exp(s^2)-1) + s )

  dCL_dmuX  <- 1/sigX
  dCL_dsigX <- -(muX - L)/sigX^2

  dCL_dmu <- dCL_dmuX*dmuX_dmu + dCL_dsigX*dsigX_dmu
  dCL_ds  <- dCL_dmuX*dmuX_ds  + dCL_dsigX*dsigX_ds

  grad <- c(dCL_dmu, dCL_ds)

  varCL <- as.numeric(t(grad) %*% V %*% grad)
  seCL <- sqrt(varCL)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seCL,
          IndexVal + qnorm(1-alpha/2)*seCL)

  list(Value=IndexVal, SE=seCL, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 4. Bootstrap Worker
###############################################################

LPIlognormal_boot_worker <- function(data, indices, L){
  x <- data[indices]
  LPIlognormal_point_est(x, L)$Value
}

###############################################################
# 5. Bootstrap CI
###############################################################

LPIlognormal_bootstrap <- function(x, L, R, alpha=0.05){

  b <- boot(x, statistic=function(data,idx)
    LPIlognormal_boot_worker(data, idx, L),
    R=R)


  ci <- boot.ci(b, conf = (1 - alpha), type=c("norm","basic","perc","bca"))

  list(
    Value = b$t0,
    norm  = ci$normal[2:3],
    basic = ci$basic[4:5],
    perc  = ci$percent[4:5],
    bca   = ci$bca[4:5]
  )
}

###############################################################
# 6. Nonparametric method
###############################################################

LPIlognormal_nonparam <- function(x, L, alpha=0.05){

  mu_hat <- mean(x)
  sig_hat <- sd(x)
  CL <- (mu_hat - L)/sig_hat

  n <- length(x)
  se_mu <- sig_hat/sqrt(n)
  se_sig <- sig_hat*sqrt(1/(2*(n-1)))

  dCL_dmu <- 1/sig_hat
  dCL_dsig <- -(mu_hat-L)/sig_hat^2

  varCL <- dCL_dmu^2*se_mu^2 + dCL_dsig^2*se_sig^2
  seCL <- sqrt(varCL)

  CI <- c(CL - qnorm(1-alpha/2)*seCL,
          CL + qnorm(1-alpha/2)*seCL)

  list(Value=CL, SE=seCL, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 7. PRINT
###############################################################

LPIlognormal_print <- function(TAB, x, alpha){

  cat("=====================================\n")
  cat("     LPI-lognormal (Classic CL Index)\n")
  cat(sprintf("     Confidence Level: %.0f%%\n", (1-alpha)*100))
  cat("     Formula: CL = (mu - L) / sigma\n")
  cat("=====================================\n\n")

  mle_par <- LPIlognormal_mle_helper(x)
  cat(sprintf("   mu_hat     : %.4f\n", mle_par$mu))
  cat(sprintf("   sigma_hat  : %.4f\n\n", mle_par$sigma))

  cat(sprintf("MLE Value           : %.4f\n\n",
              TAB$Value[TAB$Item=="MLE_Value"]))

  cat(sprintf("Asymptotic CI       : ( %.4f , %.4f )\n\n",
              TAB$Value[TAB$Item=="Asymp_LCL"],
              TAB$Value[TAB$Item=="Asymp_UCL"]))

  cat("Bootstrap CIs:\n")
  cat(sprintf("  Normal            : ( %.4f , %.4f )\n",
              TAB$Value[TAB$Item=="Boot_norm_LCL"],
              TAB$Value[TAB$Item=="Boot_norm_UCL"]))
  cat(sprintf("  Basic             : ( %.4f , %.4f )\n",
              TAB$Value[TAB$Item=="Boot_basic_LCL"],
              TAB$Value[TAB$Item=="Boot_basic_UCL"]))
  cat(sprintf("  Percentile        : ( %.4f , %.4f )\n",
              TAB$Value[TAB$Item=="Boot_perc_LCL"],
              TAB$Value[TAB$Item=="Boot_perc_UCL"]))
  cat(sprintf("  BCa               : ( %.4f , %.4f )\n\n",
              TAB$Value[TAB$Item=="Boot_bca_LCL"],
              TAB$Value[TAB$Item=="Boot_bca_UCL"]))

  cat(sprintf("Nonparametric Value : %.4f\n",
              TAB$Value[TAB$Item=="Nonpar_Value"]))
  cat(sprintf("Nonparametric CI    : ( %.4f , %.4f )\n",
              TAB$Value[TAB$Item=="Nonpar_LCL"],
              TAB$Value[TAB$Item=="Nonpar_UCL"]))
  cat("=====================================\n")
}

###############################################################
# 8. VIEWER
###############################################################

LPIlognormal_viewer <- function(summary_table, x, alpha, explanation="Lognormal Distribution Analysis") {

  mle_par <- LPIlognormal_mle_helper(x)
  conf_pct <- (1 - alpha) * 100


  items <- c("Location (mu_hat)", "Scale (sigma_hat)",
             "MLE Value", "Asymptotic CI",
             "Bootstrap Normal", "Bootstrap Basic", "Bootstrap Percentile", "Bootstrap BCa",
             "Nonparametric Value", "Nonparametric CI")

  vals <- c(
    sprintf("%.4f", mle_par$mu),
    sprintf("%.4f", mle_par$sigma),
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
    tags$h3("LPILognormal Summary", style="margin-bottom:5px;"),
    tags$p(style="font-size:16px; color:#2c5aa0; font-weight:bold; margin-top:0px;",
           sprintf("Lognormal Distribution Analysis (Confidence Level: %g%%)", conf_pct)),
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
# 9. MAIN FUNCTION
###############################################################

LPIlognormal <- function(x, L, R=500, alpha=0.05){

  mle     <- LPIlognormal_point_est(x, L)
  asym    <- LPIlognormal_asymptotic(x, L, alpha=alpha)
  bootci  <- LPIlognormal_bootstrap(x, L, R, alpha=alpha)
  nonpar  <- LPIlognormal_nonparam(x, L, alpha=alpha)

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

  LPIlognormal_print(TAB, x, alpha)
  print(LPIlognormal_viewer(TAB, x, alpha))

  invisible(TAB)
}

