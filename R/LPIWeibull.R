##############################################################
# LPIWeibull ??? Classic CL Index with Full Bootstrap CI
# Formula:  CL = (mu - L) / sigma
##############################################################

library(boot)
library(DT)
library(htmltools)

###############################################################
# 1. Weibull MLE
###############################################################

LPIWeibull_mle_helper <- function(x){

  nll <- function(par){
    k <- par[1]; lam <- par[2]
    if(k <= 0 || lam <= 0) return(1e10)
    -sum(dweibull(x, shape=k, scale=lam, log=TRUE))
  }

  fit <- nlm(nll, p=c(1.5, mean(x)), hessian=TRUE)
  k_hat <- fit$estimate[1]
  lam_hat <- fit$estimate[2]
  V <- solve(fit$hessian)

  list(shape=k_hat, scale=lam_hat, vcov=V)
}

###############################################################
# 2. Point Estimation: mu, sigma, CL
###############################################################

LPIWeibull_point_est <- function(x, L){

  out <- LPIWeibull_mle_helper(x)
  k <- out$shape
  lam <- out$scale

  mu    <- lam * gamma(1 + 1/k)
  sigma <- lam * sqrt( gamma(1 + 2/k) - gamma(1 + 1/k)^2 )

  CL <- (mu - L) / sigma

  list(Value=CL, mu=mu, sigma=sigma, shape=k, scale=lam, vcov=out$vcov)
}

###############################################################
# 3. Asymptotic CI (Delta Method)
###############################################################

LPIWeibull_asymptotic <- function(x, L, alpha=0.05){

  out <- LPIWeibull_point_est(x, L)
  k <- out$shape; lam <- out$scale
  mu <- out$mu;   sigma <- out$sigma
  V <- out$vcov

  IndexVal <- out$Value

  g1 <- gamma(1 + 1/k)
  g2 <- gamma(1 + 2/k)

  dmu_dlam <- g1
  dmu_dk   <- lam * g1 * digamma(1 + 1/k) * (-1/k^2)

  A  <- g2 - g1^2
  dA_dk <- g2 * digamma(1+2/k)*(-2/k^2) -
    2*g1*digamma(1+1/k)*(-1/k^2)

  dsig_dlam <- sqrt(A)
  dsig_dk   <- lam*(dA_dk/(2*sqrt(A)))

  dCL_dmu  <- 1/sigma
  dCL_dsig <- -(mu-L)/sigma^2

  dCL_dlam <- dCL_dmu*dmu_dlam + dCL_dsig*dsig_dlam
  dCL_dk   <- dCL_dmu*dmu_dk   + dCL_dsig*dsig_dk

  grad <- c(dCL_dk, dCL_dlam)

  varCL <- as.numeric(t(grad) %*% V %*% grad)
  seCL  <- sqrt(varCL)

  CI <- c(IndexVal - qnorm(1-alpha/2)*seCL,
          IndexVal + qnorm(1-alpha/2)*seCL)

  list(Value=IndexVal, SE=seCL, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 4. Bootstrap CL Worker
###############################################################
#' LPIWeibull
#' @export
#' @param x random sample
#' @param L lower bound
#' @param R the number of boostrap sample
#' @param alpha confidence level
#'  LPIWeibull(x, L, R, alpha)
LPIWeibull_boot_worker <- function(data, indices, L){
  x <- data[indices]
  LPIWeibull_point_est(x, L)$Value
}

###############################################################
# 5. Bootstrap CI (4 Types)
###############################################################

LPIWeibull_bootstrap <- function(x, L, R, alpha=0.05){

  boot_stat <- function(data, idx){
    LPIWeibull_boot_worker(data, idx, L)
  }

  b <- boot(x, statistic=boot_stat, R=R)
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
# 6. Nonparametric Estimation
###############################################################

LPIWeibull_nonparam <- function(x, L, alpha=0.05){

  mu_hat  <- mean(x)
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
# 7. PRINT OUTPUT
###############################################################

LPIWeibull_print <- function(TAB, x, alpha){

  cat("=====================================\n")
  cat("     LPIWeibull (Classic CL Index)\n")
  cat(sprintf("     Confidence Level: %d%%\n", (1-alpha)*100))
  cat("      Formula: CL = (mu - L) / sigma\n")
  cat("=====================================\n\n")

  mle_par <- LPIWeibull_mle_helper(x)
  cat(sprintf("   Shape (k_hat)  : %.4f\n", mle_par$shape))
  cat(sprintf("   Scale (lam_hat): %.4f\n\n", mle_par$scale))

  cat(sprintf("MLE Value             : %.4f\n\n",
              TAB$Value[TAB$Item=="MLE_Value"]))

  cat(sprintf("Asymptotic CI         : ( %.4f , %.4f )\n\n",
              TAB$Value[TAB$Item=="Asymp_LCL"],
              TAB$Value[TAB$Item=="Asymp_UCL"]))

  cat("Bootstrap CIs:\n")
  cat(sprintf("  Normal              : ( %.4f , %.4f )\n",
              TAB$Value[TAB$Item=="Boot_norm_LCL"],
              TAB$Value[TAB$Item=="Boot_norm_UCL"]))
  cat(sprintf("  Basic               : ( %.4f , %.4f )\n",
              TAB$Value[TAB$Item=="Boot_basic_LCL"],
              TAB$Value[TAB$Item=="Boot_basic_UCL"]))
  cat(sprintf("  Percentile          : ( %.4f , %.4f )\n",
              TAB$Value[TAB$Item=="Boot_perc_LCL"],
              TAB$Value[TAB$Item=="Boot_perc_UCL"]))
  cat(sprintf("  BCa                 : ( %.4f , %.4f )\n\n",
              TAB$Value[TAB$Item=="Boot_bca_LCL"],
              TAB$Value[TAB$Item=="Boot_bca_UCL"]))

  cat(sprintf("Nonparametric Value   : %.4f\n",
              TAB$Value[TAB$Item=="Nonpar_Value"]))
  cat(sprintf("Nonparametric CI      : ( %.4f , %.4f )\n",
              TAB$Value[TAB$Item=="Nonpar_LCL"],
              TAB$Value[TAB$Item=="Nonpar_UCL"]))
  cat("=====================================\n")
}

###############################################################
# 8. VIEWER ??? NO SCROLL
###############################################################

LPIWeibull_viewer <- function(TAB, x, alpha, explanation="Weibull Distribution Analysis") {

  mle_par <- LPIWeibull_mle_helper(x)
  conf_text <- sprintf("Confidence Level: %d%%", (1-alpha)*100)

  df <- data.frame(
    Item = c(
      "Shape (k_hat)",
      "Scale (lam_hat)",
      "MLE Value",
      "Asymptotic CI",
      "Bootstrap Normal",
      "Bootstrap Basic",
      "Bootstrap Percentile",
      "Bootstrap BCa",
      "Nonparametric Value",
      "Nonparametric CI"
    ),
    Value = c(
      sprintf("%.4f", mle_par$shape),
      sprintf("%.4f", mle_par$scale),
      sprintf("%.4f", TAB$Value[TAB$Item=="MLE_Value"]),
      sprintf("( %.4f , %.4f )",
              TAB$Value[TAB$Item=="Asymp_LCL"],
              TAB$Value[TAB$Item=="Asymp_UCL"]),
      sprintf("( %.4f , %.4f )",
              TAB$Value[TAB$Item=="Boot_norm_LCL"],
              TAB$Value[TAB$Item=="Boot_norm_UCL"]),
      sprintf("( %.4f , %.4f )",
              TAB$Value[TAB$Item=="Boot_basic_LCL"],
              TAB$Value[TAB$Item=="Boot_basic_UCL"]),
      sprintf("( %.4f , %.4f )",
              TAB$Value[TAB$Item=="Boot_perc_LCL"],
              TAB$Value[TAB$Item=="Boot_perc_UCL"]),
      sprintf("( %.4f , %.4f )",
              TAB$Value[TAB$Item=="Boot_bca_LCL"],
              TAB$Value[TAB$Item=="Boot_bca_UCL"]),
      sprintf("%.4f", TAB$Value[TAB$Item=="Nonpar_Value"]),
      sprintf("( %.4f , %.4f )",
              TAB$Value[TAB$Item=="Nonpar_LCL"],
              TAB$Value[TAB$Item=="Nonpar_UCL"])
    )
  )

  explanation_html <- tags$div(
    tags$h3("LPIWeibull Summary", style="margin-bottom:5px;"),
    tags$p(style="font-size:16px; color:#0056b3; font-weight:bold; margin-top:0px;", explanation, paste0("(", conf_text, ")")),
    tags$hr()
  )

  css_fix <- tags$style(HTML(".dataTables_wrapper { height: auto !important; overflow-y: hidden !important; } table.dataTable { width: 100% !important; }"))

  browsable(
    tagList(
      css_fix,
      explanation_html,
      datatable(df, rownames=FALSE,
                options=list(
                  paging=FALSE,
                  searching=FALSE,
                  ordering=FALSE,
                  scrollY=FALSE,
                  dom='t'
                ))
    )
  )
}
###############################################################
# 9. MAIN FUNCTION
###############################################################

LPIWeibull <- function(x, L, R=500, alpha=0.05){

  mle     <- LPIWeibull_point_est(x, L)
  asym    <- LPIWeibull_asymptotic(x, L, alpha=alpha)
  bootci  <- LPIWeibull_bootstrap(x, L, R, alpha=alpha)
  nonpar  <- LPIWeibull_nonparam(x, L, alpha=alpha)

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

  LPIWeibull_print(TAB, x, alpha=alpha)
  print(LPIWeibull_viewer(TAB, x, alpha=alpha))

  invisible(TAB)
}
