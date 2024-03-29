# Test third and fourth moment for normality
jarque_bera <- function(x, k=1)
{
  d <- x - mean(x)
  n <- length(x)
  S1 <- sum(d**3) / n;  S2 <- (sum(d**2) / n)**(3/2)
  C1 <- sum(d**4) / n;  C2 <- (sum(d**2) / n)**2
  JB <- ((n - k + 1) / 6) * ((S1/S2)**2 + (C1/C2 - 3)**2 / 4)
  pvalue <- 1 - pchisq(JB, df = 2)
  return(pvalue)
}

# Weighted Ljung-Box tests for correct ARMA specification
weighted_ljung_box_test <- function(x, moment=1, df)
{
  pvalues <- rugarch:::.weightedBoxTest(x, p=moment, df=df)[, "p-value"]
  out <- matrix(pvalues, dimnames=list(names(pvalues), sprintf("moment=%d", moment)))
  structure(out, moment=moment, df=df)
}

# 
weighted_arch_lm_test <- function(obj, df)
{
  # Weighted ARCH-LM tests for well fitted ARCH specification
  L2 <- rugarch:::.weightedarchlmtest(residuals(obj), sigma(obj),  lags=df+1,  fitdf=df)
  L5 <- rugarch:::.weightedarchlmtest(residuals(obj), sigma(obj),  lags=df+3,  fitdf=df)
  L10 <- rugarch:::.weightedarchlmtest(residuals(obj), sigma(obj), lags=df+5,  fitdf=df)
  alm <- matrix(0, ncol = 4, nrow = 3)
  alm[1,1:4] = as.numeric(c(L2$statistic, L2$parameter, L2$p.value))
  alm[2,1:4] = as.numeric(c(L5$statistic, L5$parameter, L5$p.value))
  alm[3,1:4] = as.numeric(c(L10$statistic, L10$parameter, L10$p.value))
  colnames(alm) = c("Statistic", "Shape", "Scale", "P-Value")
  rownames(alm) = sapply(c(1,3,5), function(x) sprintf("ARCH Lag[%d]", x))
  return(alm)
}


# Extract distribution parameters from uGARCHfit object
extract_ugarchfit_dist <- function(obj)
{
  if (!inherits(obj, "uGARCHfit"))
    stop("Object must inherit from uGARCHfit model from rugarch package.")
  model_coefs <- rugarch::coef(obj)
  dentype <- obj@model$modeldesc$distribution
  shape <- model_coefs[grep("shape", names(model_coefs))]
  skew <- model_coefs[grep("skew", names(model_coefs))]
  lambda <- model_coefs[grep("lambda", names(model_coefs))]
  return(list(dentype=dentype, shape=shape, skew=skew, lambda=lambda))
}


# provide the PIT of the residuals to the marginal models
# obj here is the @fit slot from a uGARCHfit object from the rugarch package
# 
# DEPRECATED in favor of rugarch::pit for internal consistency
PIT <- function(obj)
{
  dpars <- extract_ugarchfit_dist(obj)
  z <- residuals(obj, standardize=TRUE)
  return(pdist(distribution=dpars$dentype,
               as.numeric(z),
               lambda=dpars$lambda,
               skew=dpars$skew,
               shape=dpars$shape))
}


# Given a model produced by rugarch, pull out the residuals, perform PIT and run
# gauntlet of marginal tests
marginal_tests <- function(obj, PRINT=FALSE, PLOT=FALSE)
{
  isuGARCHfit <- inherits(obj, "uGARCHfit")
  if (isuGARCHfit) {
    # Calculate 1. standard residuals: z = (x - mu_t ) / sd_t
    #           2. probability integral transform: u
    #           3. density distribution: d
    dpars <- extract_ugarchfit_dist(obj)
    z <- residuals(obj, standardize=TRUE) # standardize refers to condition standard deviation and means
    u <- pdist(distribution=dpars$dentype, z, lambda=dpars$lambda, skew=dpars$skew, shape=dpars$shape)
    d <- ddist(distribution=dpars$dentype, z, lambda=dpars$lambda, skew=dpars$skew, shape=dpars$shape)
    armadf <- sum(obj@model$modelinc[2:3])
    gdf <- sum(obj@model$modelinc[8:9])
  } else {
    stop(sprintf('obj does not inherit from uGARCHfit class'))
  }
  
  # visual PIT transform
  if (PLOT) {
    on.exit(par(mfrow=c(1,1)))
    par(mfrow = c(1,2))
    support <- seq(-max(z), max(z), length=500)
    plot(density(z), col = "darkgoldenrod2", lwd = 2, xlab="", ylab="", main="Empirical vs theoretical")
    lines(support, ddist(dpars$dentype, support, lambda=dpars$lambda, skew=dpars$skew, shape=dpars$shape))
  }

  ####### Specification Tests
  ## Weighted LM tests of the standard residuals
  # H0: The data are independently distributed
  lm_matrix <- do.call(cbind, lapply(1:2, function(x) weighted_ljung_box_test(z, x, df=armadf)))

  # Normality Tests
  nybl <- nyblom(obj) # Parameter Stability - H0: tests the null hypothesis that all parameters are constant against the alternative of that some of the parameters are unstable
  sgnb <- signbias(obj) # Test for Asymmetric Impact of Errors - H0: null of no asymmetric effects present in the standard residuals
  gofm <- gof(obj, groups=2:5 * 10) # Distribution Check
  wtd_arch_table <- weighted_arch_lm_test(obj, df=gdf) # Test ARCH Specification
  
  ks_test <- ks.test(u, punif)     # Kolmogorov-Smirnov test (stats package)
  sw_test <- shapiro.test(as.numeric(qnorm(u))) # Shapiro Test for Normality (stats package)
  jb_test <- jarque_bera(as.numeric(qnorm(u)), length(coef(obj)))  # Jarque-Bera test for joint normality of Skew and Kurtosis
  
  if (FALSE) {
    cat("\n-------------------------------------------------------------------------")
    cat("\nWeighted Ljung-Box test on standard residuals")
    cat("\nH0: Independent Moment Condition\n")
    print(lm_matrix)
    cat("\n-------------------------------------------------------------------------")
    cat("\nWeighted ARCH-LM tested on standard residuals\n")
    print(wtd_arch_table)
    cat("\n-------------------------------------------------------------------------")
    cat("\nShapiro-Wilks normality test on the standard residuals")
    cat("\nH0: Distribution is normally distributed\n")
    msg <- sprintf("W = %.5f, p-value = %.4f", sw_test[["statistic"]], sw_test[["p.value"]])
    cat(msg,"\n")
    cat("\n-------------------------------------------------------------------------\n")
    cat("\n-------------------------------------------------------------------------\n")
    cat("Test that PIT are U(0,1)\n")
    print(gofm)
  }
  return(invisible(list(
    lm_tests=lm_matrix,
    # k_s=ks_test,
    # shapiro_wilks=sw_test,
    # jarque_bera=jb_test,
    gof=gofm,
    wtd_arch_table=wtd_arch_table,
    signbias=sgnb,
    nyblom=nybl
  )))
}


# Verify if all the distributional checks have been passed
verify_marginal_test <- function(mt, alpha=c("0.01", "0.05", "0.10"), ignore_nyblom=FALSE)
{
  alpha <- match.arg(alpha)
  alpha <- as.numeric(alpha)
  alpha.name <- sprintf("%d%%", as.integer(alpha * 100))
  
  tmpnms <- c(outer(row.names(mt$lm_tests), colnames(mt$lm_tests), FUN = function(x,y) paste(x, y)))
  p_arma_lm_test <- data.table("P-Value"=c(mt$lm_tests),
                               "Spec"=tmpnms,
                               "Test"="ARMA-LM")
  
  p_arch_lm_test <- data.table("P-Value"=mt$wtd_arch_table[, "P-Value"],
                               "Spec"=row.names(mt$wtd_arch_table),
                               "Test"="ARCH-LM")

  p_gof_test <- data.table("P-Value"=mt$gof[, "p-value(g-1)"],
                           "Spec"=paste0(mt$gof[, "group"], "bins"),
                           "Test"="GOF")
  
  p_signbias_test <- data.table("P-Value"=mt$signbias$prob,
                                "Spec"=row.names(mt$signbias),
                                "Test"="SIGN-BIAS")
  
  p_join_par_stab <- data.table("Stat"=mt$nyblom$JointStat,
                                "CV"=mt$nyblom$JointCritical[alpha.name],
                                "Spec"=alpha.name,
                                "Test"="NYBLOM-J")
  
  p_indv_par_stab <- data.table("Stat"=c(mt$nyblom$IndividualStat),
                                "CV"=mt$nyblom$IndividualCritical[alpha.name],
                                "Spec"=alpha.name,
                                "Test"=paste0("NYBLOM-I-", row.names(mt$nyblom$IndividualStat)))
  
  out <- rbindlist(list(p_arch_lm_test, p_arma_lm_test, p_gof_test,
                        p_indv_par_stab, p_join_par_stab, p_signbias_test),
                   use.names=TRUE, fill=TRUE)
  out[!is.na(`P-Value`), pass_test := `P-Value` > alpha]
  out[is.na(`P-Value`), pass_test := Stat < CV]
  if (ignore_nyblom) {
    out[grepl("NYBLOM", Test), pass_test := NA]
  }
  return(out)
}
