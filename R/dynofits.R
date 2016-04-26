#' Sequential break point copula model.
#' 
#' \code{BPfit} performs a recursive search over a bivariate time series of
#'   uniform marginal distributions.
#' 
#' @param x A numeric vector of uniform marginal values.
#' 
#' @param y A numeric vector of uniform marginal values.
#' 
#' @param fam1 An integer representing the family of the copula to use.
#'   If \code{fam2} is not \code{NULL}, this value indicates the copula family
#'   of the first copula component.
#'
#' @param fam2 A second (optional) integer indicating the use of a 
#'   mixture copula and the family of the second copula component. Defaults to
#'   \code{NULL}.
#'
#' @param parallel Logical switch to run the breakpoint search in parallel.
#'
#' @param date_names Character vector of (optional) date/timestamp names for the
#'   marginal distributions.
#' 
#' @param ncores Integer specifying the number cores to use in the parallelization.
#'   If the user specifies more cores (real or logical) than the CPU can
#'   support, the max number of supported cores will be used. 
#' 
#' @references Dias, Alexandra and Paul Embrechts, 2004, Change-point analysis
#'   for dependence structures in finance and insurance. In G. Szegoe (ed.),
#'   \emph{Risk Measures for the 21st Century}, Wiley Finance Series, 321-335
#' 
#' @references Dias, Alexandra and Paul Embrechts, 2009, Testing for structural
#' changes in exchange rates dependence beyound linear correlation,
#' \emph{European Journal of Finance}, 15(7), 619-637
#' 
#' @details The sequential break point model can be more computationally
#'   expensive than \code{\link{MSfit}} and \code{\link{STfit}}, especially
#'   as the length of the time series grows. This is due to the sequential hunt
#'   for break points.
#'   
#'   In the model's favor, robust standard errors can be calculated for each
#'   respective regime. In addition to the standard errors obtained by inverting
#'   the hessian, outer product of gradients (OPG) and sandwich estimates are
#'   also available. Standard errors based on the sandwich estimator are used
#'   in the summary.
#' 
#' @return \code{BPfit} returns an S3 object of \code{\link[base]{class}}
#'   \code{seqBreakPoint}.
#'   
#'   The summary, plot, coef, and logLik methods will give a decent snapshot of
#'   the model.
#'   
#'   An object of class \code{seqBreakPoint} is a list of lists, one for each
#'   regime. Each individual list contains the following components:
#'   \tabular{ll}{
#'     \code{pars} \tab a vector of coefficients for the copula \cr
#'     \code{log.likelihood} \tab log-likelihood value for the regime \cr
#'     \code{dep.measures} \tab a list tail dependence and Kendall's tau \cr
#'     \code{emp_hess} \tab the emprical hessian \cr
#'     \code{opg} \tab the outer product of the gradient \cr
#'     \code{sandwich} \tab the sandwich estimator \cr
#'     \code{family} \tab integers recorded which copula family to used \cr
#'     \code{points} \tab start and ending index values that subset the marginal
#'       series for the regime \cr
#'   }
#'   
#'   In addition, the following three attributes are included:
#'   
#'   \tabular{ll}{
#'     \code{marinal_names} \tab names of the marginal series \cr
#'     \code{initial_bp} \tab information on the initial break points before a
#'       re-partition method is applied \cr
#'     \code{class} \tab class of the model \cr
#'   }
#' 
#' @export
#' 

BPfit <- function(x, y, fam1, fam2 = NULL, parallel = FALSE, date_names = NULL,
                  ncores = 2) {
  
  xnme <- deparse(substitute(x))
  ynme <- deparse(substitute(y))
  stopifnot(!missing(x),
            !missing(y),
            !is.null(fam1),
            length(y) == length(x))
  
  if (!fam1 %in% c(1, 2, 3, 4, 13, 14)) {
    stop("Copula family 1 outside scope of analysis")
  }
  if (!is.null(fam2)) {
    if (!fam2 %in% c(3, 4, 13, 14))
      stop("Copula family 2 outside scope of analysis") 
  }
  n <- length(x)
  
  # here is the place to set up the parallelization stuff
  clus <- NULL
  if (parallel) {
    dtc <- max(detectCores(logical = T), detectCores(logical = F))
    if (missing(ncores)) {
      ncores <- dtc
    } else if (ncores > dtc) {
      warning(sprintf("The number of cores has been re-assigned to %d", dtc))
      ncores <- dtc
    }
    on.exit(closeAllConnections())
    clus <- makeCluster(ncores)
    clusterExport(clus, c("BreakPointLL",".cop_pdf", ".cop_cdf", ".cop_llh",
                          "cop_static", "BiCopPDF", "BiCopCDF", "BiCopEst",
                          "dependence_measures", "numerical_Ktau", "BiCopPar2TailDep",
                          "BiCopPar2Tau"))
  }
  closure <- autoBPtest(x, y, fam1, fam2, parallel, date_names, clus)
  closure$BPtest(1:n)
  repart <- if (length(closure$value()) <= 1) {
    closure$value()
  } else {
    repartitionBP(closure$value(), x, y, fam1, fam2, parallel, date_names, clus)
  }
  if (parallel) {
    stopCluster(clus)
    rm(clus)
  }
  bps <- vapply(repart, `[[`, numeric(1), 'index')
  idx <- sort(c(0, bps, 1264))
  nn <- length(bps) + 1
  output <- list()
  for (ii in 1:nn) {
    t <- idx[ii] + 1; tt <- idx[ii + 1]
    xx <- x[t:tt]; yy <- y[t:tt]
    # estimate parameters of mixture model
    res <- cop_static(xx, yy, fam1, fam2)
    # n x k matrix of gradient contributions
    G <- cop_gradcontr(xx, yy, res$pars, fam1, fam2)
    # k x k hessian
    H <- cop_hessian(xx, yy, res$pars, fam1, fam2)
    # three diffent approaches for the var-cov (VCV) matrix
    GG <- try(crossprod(G))
    # 1) Empirical Hessian
    sH <- try(solve(H))
    # 2) Sandwich Estimator
    VCV <- try(sH %*% GG %*% sH)
    # 3) Outer-Product-of-the-Gradient
    OPG <- try(solve(GG))
    output[[length(output) + 1]] <- list(
      pars = res$pars,
      log.likelihood = res$log.likelihood,
      dep.measures = dependence_measures(res$pars, fam1, fam2),
      emp_hess = -sH,
      opg = OPG,
      sandwich = VCV,
      family = c(fam1, fam2),
      points = c(start = t, end = tt)
    )
  }
  attr(output, "marginal_names") <- c(xnme, ynme)
  attr(output, "date_names") <- date_names
  attr(output, "initial_bp") <- closure$value()
  class(output) <- "seqBreakPoint"
  return(output)
}



#' Summary of the sequential break point model.
#'
#' @param Obj An object of class seqBreakPoint 
#'
#' @importFrom VineCopula BiCopName
#'
#' @export
#' 

summary.seqBreakPoint <- function(obj) {
  summarize <- function(cmp) {
    pp <- if (cmp$family == 1) cmp$pars[[1]] else cmp$pars
    n <- diff(cmp$points) + 1
    se <- sqrt(diag(cmp$sandwich))
    zstat <- pp / se
    tbl <- cbind(pp, se, zstat, 2 * pnorm(-abs(zstat)))
    colnames(tbl) <- c("Estimate", "Std.Err", "z-stat", "Pr(>|z|)")
    printCoefmat(tbl, digits = 3)
  }
  for (ii in seq_along(obj)) {
    cat("\n")
    fms <- obj[[ii]]$family
    points <- obj[[ii]]$points
    dp <- obj[[ii]]$dep.measures
    print(sprintf("Regime %d runs from index %d to %d", ii, points[1], points[2]))
    cat("\n\n")
    summarize(obj[[ii]])
    cat("\nMixture Components: \n")
    if (length(fms) == 1) {
      print(BiCopName(fms, F))
    } else {
      print(paste(BiCopName(fms[1], F), BiCopName(fms[2], F), sep = ", "))
    }
    cat("\nDependence Measures: \n")
    print(matrix(round(c(dp[[1]], dp[[2]]), 4), nrow = 3,
                 dimnames = list(c("ltd", "utd", "ktau"), "..")))
    cat("\n")
    cat("-----------------------------------------------------------\n")
    cat("-----------------------------------------------------------\n")
  }
  # lapply(lst, summarize)
}



#' Plot of the sequential break point model.
#'
#' @param Obj An object of class seqBreakPoint 
#'
#' @export
#'

plot.seqBreakPoint <- function(obj) {
  
  nms <- attr(obj, "marginal_names")
  idx <- unique(as.vector(vapply(obj, `[[`, double(2), "points")))
  idx2 <- seq(length(diff(idx))) %% 2 != 0
  ktau <- vapply(obj, function(x) x$dep.measures$Ktau, double(1))
  tdep <- vapply(obj, function(x) unlist(x$dep.measures$tail_dep), double(2))
  taudif <- apply(tdep, 2, diff)
  
  on.exit(par(mfrow = c(1, 1)))
  layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
         widths = c(1, 1), heights = c(2, 3))
  # Top Graph
  plot(rep(ktau, diff(idx)[idx2]),
       type = "l", col = "tomato",
       ylim = range(ktau) + c(-0.15, 0.15),
       ylab = "",
       xlab = "", 
       main = paste(nms, collapse = ' '))
  legend("topright", "Kendall's Tau", bty = "n")
  abline(h = 0, lty = 2)
  grid(col = "grey65")
  # Bottom Left
  plot(rep(taudif, diff(idx)[idx2]), type = "l", col = "blue",
       ylim = range(taudif) + c(-0.05, 0.05),
       xlab = "", ylab = "",
       main = paste("family:", paste(obj[[1]]$family, collapse = ", ")))
  legend("topright", "UTD - LTD", bty = "n")
  abline(h = 0, lty = 2)
  grid(col = "grey65")
  # Bottom Right
  plot(rep(tdep["upper", ], diff(idx)[idx2]), type = "l", col = "olivedrab4",
       ylim = range(tdep["upper", ]) + c(-0.025, 0.025),
       xlab = "", ylab = "")
  legend("topright", "UTD", bty = "n")
  abline(h = 0, lty = 2)
  grid(col = "grey65")
}



#' Extract the coefficients from the sequential break point model.
#'
#' @param Obj An object of class seqBreakPoint 
#'
#' @export
#' 

coef.seqBreakPoint <- function(obj) {
  out <- unlist(lapply(obj, `[[`, "pars"), F, T)
  if (obj[[1]]$family == 1)
    return(out[grepl("^mu", names(out))])
  else
    return(out)
}



#' Extract the log-likelihoods from the sequential break point model
#'
#' @param Obj An object of class seqBreakPoint 
#'
#' @export
#'

logLik.seqBreakPoint <- function(obj) {
  unlist(lapply(obj, `[[`, "log.likelihood"), F, T)
}



#' Markov-switching copula model.
#'
#' \code{MSfit} estimates a markov switching copula on a bivariate time series
#'   of uniform marginal distributions.
#' 
#' @param x A numeric vector of uniform marginal values.
#' 
#' @param y A numeric vector of uniform marginal values.
#' 
#' @param family A list of integers specifying the family of the copula to use
#'   in each regime.
#'
#' @param initValues Sometimes optional numeric vector of starting values. See
#'   Details.
#'   
#' @param tol A numeric value specifying the convergence tolerance of the model.
#'   Specifically, the model reaches convergene if the difference in likelihood
#'   from successive models falls below \code{tol}. Defaults to \code{1e-5}
#'
#' @details For \code{initValues}, if the same copula family is used for each regime, no initial values
#'   need to be supplied. If the user wants different copula families estimated
#'   in different regimes, initValues need to be supplied. For a model with K
#'   regimes, the order of the values should be provided as follows:
#'   \enumerate{
#'     \item Copula parameters in the order they appear in \code{family}
#'     \item K * (K - 1) transition variables: \eqn{p_{1,1},...,p_{1,k-1},p_{2,1},...,p_{2,k-1},...,p_{k,k-1}}
#'     \item K - 1 initial state parameters: \eqn{p_{0,1},...,p_{0,k-1} }
#'   }
#'   
#' @return \code{MSfit} returns an S3 object of \code{\link[base]{class}}
#'   \code{markovCopula}.
#'   
#'   The summary, plot, coef, and logLik functions will, repectively, print a 
#'   summarization of the output, a plot of dependence measures, extract model
#'   parameters, and extract the log-likelihood values.
#'   
#'   An object of class \code{markovCopula} has the following components: 
#'   
#'   \tabular{ll}{
#'     \code{log.likelihood} \tab log-likelihood value for the regime \cr
#'     \code{pars} \tab a vector of coefficients for the copula \cr
#'     \code{N} \tab the length of the time-series \cr
#'     \code{solver} \tab the final output from \code{\link[stats]{optim}} \cr
#'     \code{regime.inference} \tab the model's condition density, conditional
#'       probability, conditional forecasts, and the smoothed probabilities \cr
#'     \code{copula} \tab details of the estimated copulas in each regime \cr
#'     \code{transition} \tab the transition matrix and initial regime vector \cr
#'     \code{nregimes} \tab the number of regimes in the model \cr
#'   }
#'   
#' @importFrom VineCopula BiCopPar2TailDep
#'
#' @export
#'

MSfit <- function(x, y, family = list(1, 1), initValues, tol = 1e-5) {
  
  # Capture x and y names
  xnme <- deparse(substitute(x))
  ynme <- deparse(substitute(y))
  # if (any(unlist(family) > 2)) require(cubature)
  
  # Input controls
  stopifnot(length(x) == length(y), is.numeric(x),
            is.numeric(y), is.list(family),
            all(unlist(family) %in% c(1, 2, 3, 4, 13, 14)))
  if (any(vapply(family, function(x) !identical(x, family[[1]]), TRUE)) && 
        (missing(initValues) || is.null(initValues))) {
    stop("Initial values need to be supplied if you want to use different copula
         families across regimes.")
  }
  
  # Given copula type, yield the number of copula parameters per state
  nc <- vapply(family, function(x) switch(as.character(x[1]), '1' = 1, '2' = 2, 3),
               numeric(1))
  # Number of regimes equals number of elements in 'family' list
  nr <- length(family)
    
  # Make sure initValues vector is of appropriate length given copula family
  # and number of regimes. Parameters should be in the following order:
  # cop par -> transition par -> initial state par
  if (missing(initValues) || is.null(initValues)) {
    initValues <- unlist(starting_guess(x, y, unique(unlist(family)), length(family)))
    initValues <- c(initValues, rep(1, nr * (nr - 1)), c(10, rep(1, nr - 2)))
  }
  
  Pnames <- paste0(unlist(mapply(ms_pnames, family, SIMPLIFY = F)), "_",
                   unlist(mapply(function(v, r) rep(v, r), 1:nr, nc, SIMPLIFY = F)))
  Mnames <- paste0(rep("p", nr * (nr - 1)), rep(1:nr, each = nr - 1))
  Mnames <- paste0(Mnames, rep(1:(nr - 1), nr))
  Inames <- paste0("p0", 1:(nr - 1))
  names(initValues) <- c(Pnames, Mnames, Inames)
  
  if (length(initValues) != nr**2 - 1 + sum(nc))
    stop("'initValues' is not the appropriate length. Double-check copula family
         and the number of regimes.")
  
  # Iterate back and forth between copula estimation and markov matrix estimation
  res_trn <- markov_optim(initValues, x, y, family)
  res_cop <- copula_optim(res_trn$par, x, y, family)
  llh_dif <- abs(res_cop$solver$value  - res_trn$solver$value)
  llh_new <- res_cop$solver$value
  while (llh_dif > tol) {
    llh_old <- llh_new
    res_trn <- markov_optim(res_cop$par, x, y, family)
    res_cop <- copula_optim(res_trn$par, x, y, family)
    llh_new <- res_cop$solver$value
    llh_dif <- abs(llh_new - llh_old)
  }
  
  # Take output form last EM-step and store results in a list
  parm <- res_cop$par
  ncnr <- sum(nc)
  
  cop_parm <- parm[1:ncnr]
  trns_parm <- parm[(ncnr + 1):(ncnr + nr^2 - nr)]
  init_parm <- parm[(ncnr + nr^2 - nr + 1):(ncnr + nr^2 - 1)]
  
  P <- matrix(trns_parm, ncol = nr)
  P <- rbind(P^2, matrix(1, ncol = nr))
  P <- apply(P, 2, function(cc) cc / sum(cc))
  F_0 <- c(init_parm, 1)
  F_0 <- as.matrix(F_0^2 / sum(F_0^2))
  
  idx <- c(0, cumsum(nc))
  FF <- matrix(ncol = length(x), nrow = nr)
  ktau <- taildep <- vector("list", nr)
  for (rr in 1:nr) {
    aa <- idx[rr] + 1
    zz <- idx[rr + 1]
    fam <- family[[rr]]
    pp <- cop_parm[aa:zz]
    FF[rr, ] <- ms_cop_vpdf(fam, x, y)(pp)
    if (length(fam) == 1) {
      ktau[[rr]] <- BiCopPar2Tau(fam, pp[1], if (length(pp) > 1) pp[2])
      taildep[[rr]] <- unlist(BiCopPar2TailDep(fam, pp[1], if (length(pp) > 1) pp[2]))
    } else {
      ktau[[rr]] <- numerical_Ktau(pp, fam[1], fam[2])
      TD1 <- BiCopPar2TailDep(fam[1], pp[1])
      TD2 <- BiCopPar2TailDep(fam[2], pp[2])
      taildep[[rr]] <- pp[3] * unlist(TD1) + (1 - pp[3]) * unlist(TD2)
    }
  }
  
  # Regime Inference
  regime_inf <- markov_llh(P, FF, F_0)
  smooth_prob <- smooth_inf(P, regime_inf$Xi_t_t, regime_inf$Xi_t1_t)
  
  # Obtain names of copula data
  U <- data.frame(x, y, row.names = 1:length(x))
  names(U) <- c(xnme, ynme)
  
  output <- list(
    # markov_llh returns the tru logL
    log.likelihood = regime_inf$logL,
    par = res_cop$par,
    N = length(x),
    solver = res_cop$solver,
    regime.inference = list(cond.density = FF,
                            cond.probability = regime_inf$Xi_t_t,
                            cond.forecast = regime_inf$Xi_t1_t,
                            smooth.prob = smooth_prob),
    copula = list(family = family, ktau = ktau,
                  tail_dep = taildep, data = U),
    transition = list(matrix = P, init = F_0),
    nregimes = nr,
    parm.idx = idx,
    ncxnr = ncnr,
    initValues = initValues)
  
  class(output) <- "markovCopula"
  output
}



#' Plot of the markov switching model.
#'
#' @param Obj An object of class markovCopula 
#'
#' @export
#'

plot.markovCopula <- function(obj) {
  nr <- obj$nregimes
  sp <- t(obj$regime.inference$smooth.prob)
  on.exit(par(mfrow = c(1,1)))
  par(mfrow = c(nr, 1))
  for (pp in 1:nr) {
    plot(sp[, pp], type = "l", xlab = paste("regime", pp), ylab = "smooth prob.")
    grid()
  }
}



#' Summary of the markov switching model.
#'
#' @param Obj An object of class markovCopula 
#'
#' @export
#'

summary.markovCopula <- function(obj) {
  coefs <- obj$par[1:obj$ncxnr]
  se <- sqrt(diag(solve(obj$solver$hessian)))
  tstat <- coefs / se
  sumtable <- cbind(round(coefs, 3), round(se, 3), round(tstat, 3),
                    round(pchisq(abs(tstat), df = 1, 0, F), 4))

  colnames(sumtable) <- c("Estimate", "Std. Err", "t-test", "Pr(>|t|)")
  printCoefmat(as.matrix(sumtable))
  
  fam <- obj$copula$family
  nr <- length(fam)
  names(fam) <- paste0("family.", 1:nr)
  
  cat("------------------------------------------------\n")
  tranM <- obj$transition$matrix
  initM <- obj$transition$init
  dimnames(tranM) <- list(paste0("p(", 1:nr, ",.)"), paste0("p(.,", 1:nr, ")"))
  dimnames(initM) <- list(paste0("p(", 1:nr, ")"), "...")
  cat("Transition Matrix\n\n")
  print(round(tranM, 3))
  cat("\nInitial State\n")
  print(round(initM, 3)) 
  cat("------------------------------------------------\n")
  cat("Copula Family\n\n")
  print(fam)
  cat("------------------------------------------------\n")
  cat("Log-Likelihood\n")
  print(logLik(obj))
  cat("------------------------------------------------\n")
  
  ktau <- obj$copula$ktau
  tail <- obj$copula$tail_dep
  mat <- matrix(NA_real_, nrow = nr, ncol = 3,
                dimnames = list(paste0("reg.", 1:nr), c("ktau", "ltd", "utd")))
  mat[, 1] <- round(unlist(ktau), 3)
  mat[, 2] <- round(vapply(tail, function(x) x[grep("^lower", names(x))],
                           numeric(1)), 3)
  mat[, 3] <- round(vapply(tail, function(x) x[grep("^upper", names(x))],
                           numeric(1)), 3)
  cat("Dependence Measures\n")
  print(mat)
  
}



#' Extract the coefficients from the markov switching model.
#'
#' @param Obj An object of class markovCopula 
#'
#' @export
#'

coef.markovCopula <- function(obj) {
  obj$par
}



#' Extract the log-likelihoods from the markov switching model.
#'
#' @param Obj An object of class markovCopula 
#'
#' @export
#'

logLik.markovCopula <- function(obj) {
  obj$log.likelihood
}


#' Smooth-transition copula model.
#' 
#' \code{STfit} estimates a smooth-transition copula on a bivariate time series
#'   of uniform marginal distributions.
#' 
#' @param x A numeric vector of uniform marginal values.
#' 
#' @param y A numeric vector of uniform marginal values.
#' 
#' @param family A vector of integers specifying the family of the copula to use.
#' 
#' @param regimes An integers specifying how many regimes to estimate.
#' 
#' @param initValues Optional starting values. See Details.
#' 
#' @details For \code{initValues}, if the user has specific insight into any
#'   distributional changes and wants to provide specific initial parameters,
#'   for a model with K regimes the order of the values should be provided as
#'   follows:
#'   \enumerate{
#'     \item K first parameters of each regime
#'     \item K second parameters of each regime
#'     \item (Optional) weight parameter for a mixture copula
#'     \item K - 1 parameters governing speed of transition
#'     \item K - 1 location parameters
#'   }
#'
#' @return \code{STfit} returns an S3 object of \code{\link[base]{class}}
#'   \code{smoothTransCopula}.
#'   
#'   The summary, plot, coef, and logLik functions will, repectively, print a 
#'   summarization of the output, a plot of dependence measures, extract model
#'   parameters, and extract the log-likelihood values.
#'   
#'   An object of class \code{smoothTransCopula} has the following components: 
#'   
#'   \tabular{ll}{
#'     \code{log.likelihood} \tab log-likelihood value for the regime \cr
#'     \code{pars} \tab a vector of coefficients for the copula \cr
#'     \code{dep.measures} \tab list containing tail dependence and Kendall's tau \cr
#'     \code{smooth.parameters} \tab matrix of smooth parameter paths \cr
#'     \code{N} \tab the length of the time-series \cr
#'     \code{solver} \tab the output from \code{\link[Rsolnp]{solnp}} \cr
#'     \code{copula} \tab details of the estimated copulas in each regime \cr
#'     \code{transition} \tab the transition matrix and initial regime vector \cr
#'     \code{nregimes} \tab the number of regimes in the model \cr
#'   }
#'
#'
#' @importFrom Rsolnp solnp
#'
#' @export
#'

STfit <- function(x, y, family = 1, regimes = 2, initValues = NULL) {
  
  # Capture x and y names
  xnme <- deparse(substitute(x))
  ynme <- deparse(substitute(y))
  k <- regimes
  
  # Input controls
  stopifnot(length(x) == length(y),
            is.numeric(x), 
            is.numeric(y),
            all(family %in% c(1, 2, 3, 4, 13, 14)),
            regimes > 1)
  nc <- switch(as.character(family)[1], '1' = 1,	'2' = 2, 3)
  nct <- if (length(family) == 2 && (all(family > 2))) {
    4 * k - 1
  } else {
    k * nc + 2 * (k - 1)
  }
  
  if (missing(initValues) || is.null(initValues)) {
    initValues <- unlist(starting_guess(x, y, family, k))
    inflects <- (c(1:k - 1) / k)[-1]
    deltas <- rep(5, k - 1)
    
    if (length(family) == 1 && family == 1) {
      initValues <- c(initValues, deltas, inflects)
      
    } else if (length(family) == 1 && family == 2) {
      initValues <- c(initValues[grepl("^mu", names(initValues))],
                      initValues[grepl("^nu", names(initValues))],
                      deltas, inflects)
      
    } else if (length(family) == 2) {
      initValues <- c(initValues[grepl("^th1", names(initValues))],
                      initValues[grepl("^th2", names(initValues))],
                      mean(initValues[grepl("^wt", names(initValues))]),
                      deltas, inflects)
    }
  }
  
  if (length(initValues) != nct)
    stop(sprintf("'initValues' is not the appropriate length. Double-check copula
                 family and the number of regimes. With the specified values for 'family'
                 and 'regimes', the function is expecting %d initial parameters.", nct))
  
  # Create transition variable.
  N <- length(x)
  St <- cumsum(rep(1, N)) / N
  
  # Function that produces the log-likelihood of the model
  stLLH <- function(parm, family, k, St) {
    TH <- transform_parm(parm, family, k, St, sd(St))
    smooth_llh(x, y, TH, family)
  }
  
  bounds <- st_boundary(family, k)
  # Constraint function for the transition locations
  ineqfn <- function(parm, family, k, St) { 
    T <- length(parm)
    x <- parm[(T - k + 2):T]
    diff(x)
  }
  
  # require(Rsolnp)
  if (k == 2) {
    # No inequality bounds
    Fit <- solnp(initValues, stLLH, control = list(trace = 1),
                 LB = bounds$LB, UB = bounds$UB, family = family, k = k,
                 St = St)
  } else {
    # Inequality bounds need for inflection points
    Fit <- solnp(initValues, stLLH, control = list(trace = 1),
                 LB = bounds$LB, UB = bounds$UB, ineqfun = ineqfn,
                 ineqLB = rep(1 / N, k - 2), ineqUB = rep(1 - 1 / N, k - 2),
                 family = family, k = k, St = St)
  }
  
  smoothTrans <- transform_parm(Fit$pars, family, k, St, sd(St))
  names(Fit$pars) <- bounds$pnames
  
  nfk <- length(family) * k
  if (nc %in% c(1, 3)) {
    pp <- split(Fit$pars[1:nfk],
                ceiling(seq_along(Fit$pars[1:nfk]) / ifelse(nc == 1, 1, 2)))
    if (nc == 3) {
      pp <- lapply(pp, function(x) c(x, Fit$pars['wt']))
    }
  } else if (nc == 2) {
    pp <- lapply(1:k, function(x) Fit$pars[c(x, x + k)])
  }
  
  f1 <- family[[1]]
  f2 <- if (length(family) == 2) family[[2]]
  depMeas <- lapply(pp, dependence_measures, fam1 = f1, fam2 = f2)
  U = data.frame(x, y, row.names = 1:length(x))
  names(U) <- c(xnme, ynme)
  output <- list(
    log.likelihood = Fit$values,
    pars = round(Fit$pars, 6),
    dep.measures = depMeas,
    smooth.parameters = smoothTrans,
    solver = Fit,
    family = family,
    nregimes = k,
    N = N,
    npars = nc,
    U = U,
    initValues = initValues)
  
  class(output) <- "smoothTransCopula"
  output
}



#' Plot of the smooth transition model.
#'
#' @param Obj An object of class smoothTransCopula 
#'
#' @export
#'

plot.smoothTransCopula <- function(obj) {
  smooth.pars <- obj$smooth.parameters
  nmes <- names(obj$U)
  fam <- obj$family
  nc <- ncol(obj$smooth.parameters)
  ylabs <- switch(as.character(nc), "1" = "mu", "2" = c("mu", "v"),
                  c("th1", "th2"))
  on.exit(par(mfrow = c(1, 1)))
  par(mfrow = c(ifelse(nc == 1, 1, 2), 1))
  for (cc in 1:nc) {
    plot(smooth.pars[ ,cc], xlab = "",
         ylab = ylabs[cc],
         main = c(paste(nmes, collapse = " "), "")[cc])
    grid()
  }
}



#' Summary of the smooth transition model.
#'
#' @param Obj An object of class smoothTransCopula 
#'
#' @export
#'

summary.smoothTransCopula <- function(obj) {
  coefs <- obj$pars
  fam <- obj$family
  nr <- obj$nregimes
  se <- sqrt(diag(solve(obj$solver$hessian)))
  if (nr > 2) {
    se <- se[-(1:(nr - 2))]
  }
  tstat <- coefs / se
  sumtable <- cbind(round(coefs, 3), round(se, 3), round(tstat, 3),
                    round(pchisq(abs(tstat), df = 1, 0, F), 4))
  colnames(sumtable) <- c("Estimate", "Std.Err", "z-stat", "Pr(>|t|)")
  row.names(sumtable) <- names(coefs)
  printCoefmat(sumtable)
  
  cat("------------------------------------------------\n")
  cat("Dependence Measures\n")
  depMeas <- t(sapply(obj$dep.measures, function(x) c(x$Ktau, x$tail_dep)))
  colnames(depMeas) <- c("Ktau", "ltd", "utd")
  row.names(depMeas) <- paste0("regime.", 1:nr)
  print(round(depMeas, 3))
  
  cat("------------------------------------------------\n")
  cat("Log-Likelihood\n")
  print(logLik(obj))
  cat("------------------------------------------------\n")
}



#' Extract the log-likelihoods from the smooth transition model.
#'
#' @param Obj An object of class smoothTransCopula 
#'
#' @export
#'

logLik.smoothTransCopula <- function(obj) {
  obj$log.likelihood[length(obj$log.likelihood)]
}



#' Extract the coefficients from the smooth transition model.
#'
#' @param Obj An object of class smoothTransCopula 
#'
#' @export
#'

coef.smoothTransCopula <- function(obj) {
  obj$pars
}


#' Compare different models using AIC/BIC
#' 
#' @param ... One or more dynamic copula models to compare. They can be comma
#'   seperated or in a list.
#'
#' @importFrom VineCopula BiCopName
#'
#' @export
#'

aic_bic <- function(...) {
  
  objsLst <- list(...)
  sIn <- mapply(inherits, objsLst, "smoothTransCopula")
  N <- vapply(objsLst, `[[`, integer(1), 'N')
  stopifnot(all(sIn), all(N == N[1]))
  # Information Criterion
  logL <- vapply(objsLst,
                 function(x) x$log.likelihood[length(x$log.likelihood)],
                 double(1))
  k <- vapply(objsLst, function(x) length(x$pars), integer(1))
  aic <- 2 * k - 2 * (-logL)
  bic <- log(unique(N)) * k - 2 * (-logL)
  
  nregimes <- vapply(objsLst, `[[`, double(1), "nregimes")
  fam_int <- vapply(objsLst, function(x) x$family[1], double(1))
  copula_names <- vapply(fam_int, BiCopName, character(1), FALSE)
  margs <- paste0(names(objsLst[[1]]$U), collapse = "")
  outnms <- paste0(paste0(copula_names, toupper(margs)), nregimes)
  
  compare <- function(crit) {
    relativeProb <- exp((min(crit) - crit) / 2)
    out <- matrix(c(logL, crit, relativeProb), nrow = 3, byrow = T)
    dimnames(out) <- list(c("LogL", "Criterion", "Relative Prob."),
                          outnms)
    round(out, 3)
  }
  
  out <- lapply(list(aic, bic), compare)
  names(out) <- c("AIC", "BIC")
  out
}
