##############################################################
# LPIgamma — Classic CL Index with Full Bootstrap CI
# Formula:  CL = (mu - L) / sigma   (Gamma Moments)
##############################################################

library(boot)
library(DT)
library(htmltools)
library(MASS)

###############################################################
# 1. Gamma MLE (shape alpha, scale theta)
###############################################################
#' LPIGamma
#' @export
#' @param x random sample
#' @param L lower bound
#' @param R the number of boostrap sample
#'LPIgamma(x, L, R, alpha)

LPIgamma_mle_helper <- function(x){

  fit <- fitdistr(x, densfun="gamma")
  alpha_hat <- fit$estimate["shape"]
  theta_hat <- 1/fit$estimate["rate"]   # MASS returns rate = 1/scale

  V <- fit$vcov

  list(shape=alpha_hat, scale=theta_hat, vcov=V)
}

###############################################################
# 2. Theoretical mu, sigma, CL
###############################################################

LPIgamma_point_est <- function(x, L){

  out <- LPIgamma_mle_helper(x)
  a <- out$shape
  th <- out$scale

  mu  <- a*th
  sig <- th*sqrt(a)

  CL <- (mu - L)/sig

  list(Value=CL, mu=mu, sigma=sig, shape=a, scale=th, vcov=out$vcov)
}

###############################################################
# 3. Asymptotic CI (Delta Method)
###############################################################

LPIgamma_asymptotic <- function(x, L, alpha=0.05){

  out <- LPIgamma_point_est(x, L)
  a <- out$shape
  th <- out$scale
  mu <- out$mu
  sig <- out$sigma
  V <- out$vcov

  IndexVal <- out$Value

  # mu = a*th ---> derivatives
  dmu_da  <- th
  dmu_dth <- a

  # sigma = th*sqrt(a)
  dsig_da  <- th/(2*sqrt(a))
  dsig_dth <- sqrt(a)

  # CL = (mu - L)/sig
  dCL_dmu  <- 1/sig
  dCL_dsig <- -(mu-L)/sig^2

  # chain rule
  dCL_da  <- dCL_dmu*dmu_da  + dCL_dsig*dsig_da
  dCL_dth <- dCL_dmu*dmu_dth + dCL_dsig*dsig_dth

  grad <- c(dCL_da, dCL_dth)

  varCL <- as.numeric(t(grad) %*% V %*% grad)
  seCL <- sqrt(varCL)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seCL,
          IndexVal + qnorm(1-alpha/2)*seCL)

  list(Value=IndexVal, SE=seCL, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 4. Bootstrap Worker
###############################################################

LPIgamma_boot_worker <- function(data, indices, L){
  x <- data[indices]
  LPIgamma_point_est(x, L)$Value
}

###############################################################
# 5. Bootstrap CI (4 types)
###############################################################

LPIgamma_bootstrap <- function(x, L, R, alpha=0.05){

  boot_stat <- function(data, idx){
    LPIgamma_boot_worker(data, idx, L)
  }

  b <- boot(x, statistic=boot_stat, R=R)


  ci <- boot.ci(b, conf = 1-alpha, type=c("norm","basic","perc","bca"))

  list(
    Value = b$t0,
    norm  = ci$normal[2:3],
    basic = ci$basic[4:5],
    perc  = ci$percent[4:5],
    bca   = ci$bca[4:5]
  )
}

###############################################################
# 6. Nonparametric Estimation
###############################################################

LPIgamma_nonparam <- function(x, L, alpha=0.05){

  mu_hat <- mean(x)
  sig_hat <- sd(x)
  CL <- (mu_hat - L)/sig_hat

  n <- length(x)

  se_mu <- sig_hat/sqrt(n)
  se_sig <- sig_hat*sqrt(1/(2*(n-1)))

  dCL_dmu  <- 1/sig_hat
  dCL_dsig <- -(mu_hat-L)/sig_hat^2

  varCL <- dCL_dmu^2*se_mu^2 + dCL_dsig^2*se_sig^2
  seCL <- sqrt(varCL)

  CI <- c(CL - qnorm(1-alpha/2)*seCL,
          CL + qnorm(1-alpha/2)*seCL)

  list(Value=CL, SE=seCL, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 7. PRINT OUTPUT
###############################################################

LPIgamma_print <- function(TAB, x, alpha=0.05){

  conf_pct <- (1-alpha)*100

  cat("=====================================\n")
  cat("        LPIgamma (Classic CL Index)\n")
  cat(sprintf("   Confidence Level: %.0f%%\n", conf_pct))
  cat("       Formula: CL = (mu - L)/sigma\n")
  cat("=====================================\n\n")

  mle <- LPIgamma_mle_helper(x)
  cat(sprintf("   Shape (alpha_hat) : %.4f\n", mle$shape))
  cat(sprintf("   Scale (theta_hat) : %.4f\n\n", mle$scale))

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
# 8. VIEWER — LPIGamma Summary
###############################################################

LPIgamma_viewer <- function(summary_table, x, alpha, explanation="Gamma Distribution Analysis") {

  mle <- LPIgamma_mle_helper(x)
  conf_level <- (1 - alpha) * 100


  items <- c("Shape (k_hat)", "Scale (lam_hat)",
             "MLE Value", "Asymptotic CI",
             "Bootstrap Normal", "Bootstrap Basic", "Bootstrap Percentile", "Bootstrap BCa",
             "Nonparametric Value", "Nonparametric CI")


  vals <- c(
    sprintf("%.4f", mle$shape),
    sprintf("%.4f", mle$scale),
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
    tags$h3("LPIGamma Summary", style="margin-bottom:5px;"),
    tags$p(style="font-size:16px; color:#0056b3; margin-top:0px; font-weight:bold;", info_text),
    tags$hr()
  )

  css_fix <- tags$style(HTML("
    .dataTables_wrapper { height: auto !important; overflow-y: hidden !important; }
    table.dataTable { width: 100% !important; }
  "))

  return(browsable(tagList(
    css_fix,
    explanation_html,
    DT::datatable(df, rownames = FALSE,
                  options = list(
                    pageLength = nrow(df),
                    scrollY = FALSE,
                    paging = FALSE,
                    dom = 't',
                    ordering = FALSE
                  ))
  )))
}

###############################################################
# 9. MAIN FUNCTION
###############################################################

LPIgamma <- function(x, L, R=500, alpha=0.05){

  mle      <- LPIgamma_point_est(x, L)
  asym     <- LPIgamma_asymptotic(x, L, alpha=alpha)
  bootci   <- LPIgamma_bootstrap(x, L, R, alpha=alpha)
  nonpar   <- LPIgamma_nonparam(x, L, alpha=alpha)

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

  LPIgamma_print(TAB, x, alpha)
  print(LPIgamma_viewer(TAB, x, alpha))

  invisible(TAB)
}
