# Calculate density of two-component mixture copula
# 
# @param theta A vector of length 3. The first item is the (single) parameter
#    for the first component copula, the second item is the (single) parameter
#    for the second component copula, and the third is the copula weight.
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param fam1 An integer representing the family of the component copula applied
#    to u. See Details.
#
# @param fam2 An integer representing the family of the component copula applied
#    to v. See Details.
#    
# @return A numeric vector of the same length as \code{u} and \code{v} whose
#    values are the density (pdf) of the given mixture copula
#
# @importFrom VineCopula BiCopPDF
#
 

cop_pdf <- function(theta, u, v, fam1, fam2) {
  th1 = theta[1]; th2 = theta[2]; th3 = theta[3]
  stopifnot(th3 >= 0 && th3 <= 1)
  
  th3 * BiCopPDF(u, v, fam1, th1) + (1 - th3) * BiCopPDF(u, v, fam2, th2)
}  



# Calculate distribution of two-component mixture copula
# 
# @param theta A vector of length 3. The first item is the (single) parameter
#    for the first component copula, the second item is the (single) parameter
#    for the second component copula, and the third is the copula weight.
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param fam1 An integer representing the family of the component copula applied
#    to u. See Details.
#
# @param fam2 An integer representing the family of the component copula applied
#    to v. See Details.
# 
# @return A numeric vector of the same length as \code{u} and \code{v} whose
#    values are the distribution (cdf) of the given mixture copula
#
# @importFrom VineCopula BiCopCDF
#

cop_cdf <- function(theta, u, v, fam1, fam2) {
  th1 = theta[1]; th2 = theta[2]; th3 = theta[3]
  stopifnot(th3 >= 0 && th3 <= 1,
            length(u) == length(v))
  
  th3 * BiCopCDF(u, v, fam1, th1) + (1 - th3) * BiCopCDF(u, v, fam2, th2)
}



# Convienence function that returns the sum of the negative log-likelihood of
# a mixture copula
# 
# @param theta A vector of length 3. The first item is the (single) parameter
#    for the first component copula, the second item is the (single) parameter
#    for the second component copula, and the third is the copula weight.
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param fam1 An integer representing the family of the component copula applied
#    to u. See Details.
#
# @param fam2 An integer representing the family of the component copula applied
#    to v. See Details.
#  
# @return The negative log-likelhood value
#    

cop_llh <- function(theta, u, v, fam1, fam2) {

  -sum(log(cop_pdf(theta, u, v, fam1, fam2)))
}



#' Estimates a static copula
#' 
#' @param u A vector of uniform marginal values.
#' 
#' @param v A vector of uniform marginal values.
#' 
#' @param fam1 An integer representing the family of the component copula applied
#'    to u. See Details.
#'
#' @param fam2 An integer representing the family of the component copula applied
#'    to v. See Details.
#'  
#' @return A list of: \enumerate{
#'   \item the log-likelihood value
#'   \item the estimated copula parameters
#'   }
#'
#' @importFrom VineCopula BiCopEst
#' @importFrom VineCopula BiCopPDF
#' @export
#' 

cop_static <- function(u, v, fam1, fam2 = NULL) {

  stopifnot(length(u) == length(v))
  isNullFam2 <- is.null(fam2)
  
  unpack_ <- function(fit)
  {
    return(list(par=fit[["par"]], par2=fit[["par2"]]))
  }
  
  # Non-mixture Copula
  if (isNullFam2) {
    Fit <- unpack_(BiCopEst(u, v, fam1, max.df = 50))
    Lc <- log(BiCopPDF(u, v, fam1, Fit$par, par2 = if (length(Fit) == 2L) Fit$par2 else 0))
    theta <- unlist(Fit)
    names(theta) <- if (fam1 < 3) {c("mu", "nu")} else {c("th1", "th2")}
    LL <- sum(Lc)
  } else if (fam1 > 2 && fam2 > 2) { # Mixture Copulas
    # Set lower bounds for parameters
    S = 1e-7
    LB <- vector("integer", 3)
    LB[1] <- switch(as.character(fam1), "3" = S, "4" = 1 + S, "13" = S, "14" = 1 + S)
    LB[2] <- switch(as.character(fam2), "3" = S, "4" = 1 + S, "13" = S, "14" = 1 + S)
    LB[3] <- S
    Fit <- optim(c(1.5, 1.5, 0.5), cop_llh, u = u, v = v, fam1 = fam1, fam2 = fam2,
                 method = "L-BFGS-B", lower = LB, upper = c(30, 30, 1 - S))
    # Return Log-like as max value
    LL <- -Fit$value
    theta <- Fit$par
    names(theta) <- c("th1", "th2", "wt")
  }
  out <- list(
    log.likelihood = LL,
    pars = theta,
    N = length(u),
    family = c(fam1, fam2),
    nregimes = 1L,
    dep.measures = dependence_measures(theta, fam1, fam2)
  )
  class(out) <- "cop_static"
  out
}



# For a mixture copula, estimate Kendall's tau using numerical integration
#  
# @param theta A vector of length 3. The first item is the (single) parameter
#    for the first component copula, the second item is the (single) parameter
#    for the second component copula, and the third is the copula weight.
#    
# @param fam1 An integer representing the family of the first component copula.
#    
# @param fam1 An integer representing the family of the second component copula.
# 
# @return Kendall's tau for the given mixture copula
#
# @references Nelson (2006) pg ___, eq ___
# 
# @importFrom cubature adaptIntegrate
# 

numerical_Ktau <- function(theta, fam1, fam2) {
  
  Q <- function(U, pars) {
    u <- U[1]; v <- U[2]
    cop_cdf(pars, u, v, fam1, fam2) * cop_pdf(pars, u, v, fam1, fam2)
  }
  4 * adaptIntegrate(Q, lowerLimit = c(0,0), upperLimit = c(1,1),
                     pars = theta)$integral - 1
}



#' Provide semi-smart starting values for dynamic copulas.
#'
#' @param u A vector of uniform marginal values.
#' 
#' @param v A vector of uniform marginal values.
#'
#' @param fam A vector of length one or two indicating which copula family to
#'   find initial values. A vector of length two indicates a mixture copula
#'
#' @param k Integer indicating how many segments to disjoint sets to break (u, v)
#'   into. 
#' 
#' @return A list of length \code{k} with estimated copula parameters
#' 
#' @details This function breaks the given fixed marginals into \code{k} disjoint
#'   but adjacent sets of equal size. It then returns the estimated parameters of
#'   a static copula governed by copula family \code{fam} for each of the disjoint
#'   sets.
#'
#' @export
#' 

starting_guess <- function(u, v, fam, k) {
  
  N <- length(u)
  stopifnot(length(v) == N,
            length(fam) <= 2)
  idx <- round((1:(k - 1) / k) * N)
  ed <- c(idx, N)
  st <- c(0, idx) + 1
  idx <- mapply(function(aa, zz) seq(aa, zz), st, ed, SIMPLIFY = FALSE)
  X <- lapply(idx, function(ii) u[ii])
  Y <- lapply(idx, function(ii) v[ii])
  f1 <- lapply(1:k, function(xx) fam[[1]])
  f2 <- lapply(1:k, function(xx) if(length(fam) == 2) fam[[2]])
  
  cop_est2 <- function(x, y, ff1, ff2) {
    theta <- cop_static(x, y, ff1, ff2)[["pars"]]
    if (ff1 == 1)
      theta <- theta[1]
    theta
  }
  
  out <- mapply(cop_est2, X, Y, f1, f2, SIMPLIFY = FALSE)
  out
}



#' Calculates Kendall's tau, upper tail dependence, and lower tail dependence
#' for a given copula.
#' 
#' @param pars A vector of relevant parameter values
#' 
#' @param fam1 An integer representing the family of the first component copula.
#' 
#' @param fam2 An integer representing the family of the second component copula.
#'   Defaults to \code{NULL}
#'
#' @return A list of: \enumerate{
#'   \item upper and lower tail-dependence values
#'   \item Kendall's tau
#'   } 
#' 
#' @details Providing only \code{fam1} indicates that the user wants the
#'   dependence values for a single coplula. Providing both \code{fam1} and
#'   \code{fam2} will return the dependence measures for a mixture copula.
#' 
#' @importFrom VineCopula BiCopPar2TailDep
#' @importFrom VineCopula BiCopPar2Tau
#' @importFrom VineCopula BiCopCDF
#' @importFrom VineCopula BiCopPDF
#' @importFrom cubature adaptIntegrate
#' 
#' @export
#' 

dependence_measures <- function(pars, fam1, fam2 = NULL) {
  
  if (is.null(fam2)) {
    tdep <- unlist(BiCopPar2TailDep(fam1, pars[[1]], if (length(pars) > 1) pars[[2]] else 0))
    ktau <- BiCopPar2Tau(fam1, pars[[1]], if (length(pars) > 1) pars[[2]] else 0)
    
  } else if (length(pars) == 3) {
    convex_weights <- function(x, y, z)  x * z + y * (1 - z)
    th1 <- pars[1]; th2 <- pars[2]; wt <- pars[3]
    if (wt <= 0 || wt >= 1) {
      warning(sprintf("wt parameter takes a value outside the unit interval: %", wt))
    }
    tdep <- mapply(convex_weights, BiCopPar2TailDep(fam1, th1), 
                   BiCopPar2TailDep(fam2, th2), wt) 
    ktau <- numerical_Ktau(c(th1, th2, wt), fam1, fam2)
    
  } else {
    stop("Function doesn't know how to handle this particular combination of 
         family and parameter values.")
  }
  names(ktau) <- "Ktau"
  list(tail_dep = tdep, Ktau = ktau)
}



# For the smooth transition model, create boundaries with respect to copula
# model to pass on to the optimization function
# 
# @param fam Numeric vector of length one or two indicating the copula type and
#   family.
# 
# @param k Number of regimes
#

st_boundary <- function(fam, k) {
  
  S <- 1e-3
  if (length(fam) == 1) {
    if (fam == 1) {
      UB <- rep(1 - S, k)
      LB <- rep(S - 1, k)
      names(UB) <- names(LB) <- paste0("mu", 1:k)
    }
    else if (fam == 2) {
      UB <- rep(c(1 - S, 100), each = k)
      LB <- rep(c(S - 1, 2 + S), each = k)
      names(UB) <- names(LB) <- c(paste0("mu", 1:k), paste0("v", 1:k))  
    }
  }
  else if (length(fam) == 2) {
    l1 <- switch(as.character(fam[[1]]), "3" = 0, "4" = 1, "13" = 0, "14" = 1)
    l2 <- switch(as.character(fam[[2]]), "3" = 0, "4" = 1, "13" = 0, "14" = 1)
    UB <- c(rep(30, 2 * k), 1 - S)
    LB <- c(rep(c(l1 + S, l2 + S), each = k), S)
    names(UB) <- names(LB) <- c(paste0("th1", 1:k), paste0("th2", 1:k), "wt")
  } else {
    stop("More than two copulas specified. Currently, only a mixture
      of two Archimedean copulas is supported.")
  }
  gUB <- rep(300, k - 1); gLB <- rep(S, k - 1)
  cUB <- rep(0.95, k - 1); cLB <- rep(0.05, k - 1)
  names(gUB) <- names(gLB) <- paste0("g", 1:(k - 1))
  names(cUB) <- names(cLB) <- paste0("c", 1:(k - 1))
  UB <- c(UB, gUB, cUB); LB <- c(LB, gLB, cLB)
  list(UB = UB, LB = LB, pnames = names(UB))
}



# For the markov switching model, create boundaries with respect to copula
# model to pass on to the optimization function
# 
# @param fam Numeric vector of length one or two indicating the copula type and
#   family.
# 
# @param k Number of regimes
#

ms_boundary <- function(fam, k) {
  
  S <- 1e-2
  if (length(fam) == 1) {
    if (fam == 1)
      return(list(UB = 1 - S, LB = -1 + S))
    else if (fam == 2)
      return(list(UB = c(1 - S, 100), LB = c(-1, 2) + S))
  }
  else if (length(fam) == 2) {
    l1 <- switch(as.character(fam[[1]]), "3" = 0, "4" = 1, "13" = 0, "14" = 1)
    l2 <- switch(as.character(fam[[2]]), "3" = 0, "4" = 1, "13" = 0, "14" = 1)
    return(list(UB = c(30, 30, 1 - S), LB = c(l1 + S, l2 + S, S)))
  } else {
    stop("Input variable 'fam' needs to be of length two.")
  }
}



# For the markov switching model, create function factory that returns a
# copula density function conditioning on given copula paramaters
#
# @param fam Numeric vector of length one or two indicating the copula type and
#   family.
#
# @param x A vector of uniform marginal values.
# 
# @param y A vector of uniform marginal values.
#
# @importFrom VineCopula BiCopPDF
#

ms_cop_vpdf <- function(fam, x, y) {
  ff1 <- fam[[1]]
  ff2 <- if (length(fam) > 1) fam[[2]]
  if (ff1 %in% c(1, 2) && length(fam) == 1) {
    return(function(pp) BiCopPDF(u1 = x, u2 = y, family = ff1, par = pp[[1]],
                                 par2 = if (length(pp) == 2) pp[[2]] else 0))
  } else {
    return(function(pp) cop_pdf(theta = pp, u = x, v = y, fam1 = ff1, fam2 = ff2))
  }
}



# For the smooth transition model, create function factory that returns a
# copula density function conditioning on given copula paramaters
# 
# @param family Numeric vector of length one or two indicating the copula type
#   and family.
#
# @param x A vector of uniform marginal values.
# 
# @param y A vector of uniform marginal values.
# 
# @return A density function for a copula for arbitrary parameter(s)
# 

st_cop_vpdf <- function(family, x, y) {
  
  stopifnot(is.numeric(family), is.numeric(x), is.numeric(y))
  if (family == 1) {
    return(function(pp) gaussian_pdf(x, y, pp[, 1]))
  } else if (family == 2) {
    return(function(pp) tcop_pdf(x, y, pp[, 1], pp[, 2]))
  } else if (family == 3) {
    return(function(pp) clayton_pdf(x, y, pp[, 1]))
  } else if (family == 4) {
    return(function(pp) gumbel_pdf(x, y, pp[, 1]))
  } else if (family == 13) {
    return(function(pp) Sclayton_pdf(x, y, pp[, 1]))
  } else if (family == 14) {
    return(function(pp) Sgumbel_pdf(x, y, pp[, 1]))
  } else {
    stop("Copula PDF not chosen. Check family variable.")
  }
}


# For the markov switching model, return copula names
# 
# @param x Integer indicating copula family
# 
# @return A vector of parameter names
# 

ms_pnames <- function(x) {
  
  if (length(x) == 1) {
    if (x == 1) {
      return("mu")
    } else if (x == 2) {
      return(c("mu", "nu"))
    }
  } else if (length(x) == 2) {
    return(c("th1", "th2", "wt")) 
  }
}


# For the markov switching model, run the algorithm that produces the model's
# log-likelihood.
#
# @param P A markov transition matrix for nr regimes
# 
# @param FF An nr x n matrix where each row represents the densities for each
#   of the state dependent regimes
# 
# @param F_0 An nr x 1 regime of initial regime probabilities
# 
# @return A list of: \describe{
#   \item{logL}{the log-likelihood of the model}
#   \item{Xi_t_t}{Xi_t_t}
#   \item{Xi_t1_t}{Xi_t1_t}
#   }
#
# @references Hamilton (1994) pg ___, eq [22.4.5] and [22.4.5]
#

markov_llh <- function(P, FF, F_0) {
  
  stopifnot(nrow(P) == nrow(FF), nrow(P) == nrow(F_0))
  Xi_t_t <- matrix(nrow = nrow(FF), ncol = ncol(FF))
  Xi_t1_t <- matrix(nrow = nrow(FF), ncol = ncol(FF) + 1)
  Xi_t1_t[, 1] <- F_0
  logL <- 0
  for (jj in 2:(ncol(FF) + 1)) {
    fp <- Xi_t1_t[, jj - 1] * FF[, jj - 1]
    sfp <- sum(fp)
    logL <- logL + log(sfp)
    Xi_t_t[, jj - 1] <- fp / sfp
    Xi_t1_t[, jj] <- P %*% Xi_t_t[, jj - 1]
  }
  list(logL = logL, Xi_t_t = Xi_t_t, Xi_t1_t = Xi_t1_t)
}



# For the markov switching model, an implementation of Kim's smooth inference
# alogrithm utilizes the entire series of data to obtain improved inference for
# what regime the model occupies at each unit of time 
# 
# @param P A markov transition matrix for nr regimes
# 
# @param Xi_t_t An nr x 1 matrix
# 
# @param Xi_t1_t An nr x 1 matrix
# 
# @return An nr x n matrix whose columns are the units of time for the time
#   series and each row is the probability of being in that regime at that time
# 
# @references Hamilton (1994) pg ___ eq [22.4.14]
# 

smooth_inf <- function(P, Xi_t_t, Xi_t1_t) {
  
  nc <- ncol(Xi_t_t); nr <- nrow(Xi_t_t)
  Xi_t_T <- matrix(nrow = nr, ncol = nc)
  Xi_t_T[, nc] <- Xi_t_t[, nc]
  Xi_t1_t <- Xi_t1_t[, -1]
  
  for (tt in (nc - 1):1) {
    Xi_t_T[, tt] <- Xi_t_t[, tt] *
      (t(P) %*% (Xi_t_T[, tt + 1] / Xi_t1_t[, tt]))
  }
  Xi_t_T
}



# For the markov switching model, given the parameters that define the
# markov transition matrix, complete the maximization step of the E-M algorithm
# by optimizing the copula parameters
# 
# @param parm A vector of copula parameters
# 
# @param x A vector of uniform marginal values.
# 
# @param y A vector of uniform marginal values. 
#
# @param family Numeric vector of length one or two indicating the copula type
#   and family.
#
# @return A list of: \describe{
#   \item{par}{Estimated copula parameters}
#   \item{solver}{Output from the call to \code{\link[stats]{optim}}}
#   \item{EMstep}{Good to know for further development}
#   }
#

copula_optim <- function(parm, x, y, family) {
  
  # Basic constants
  nc <- vapply(family, function(x) switch(as.character(x[1]), "1" = 1, "2" = 2, 3),
               numeric(1))
  nr <- length(family)
  nx <- length(x)
  ncnr <- sum(nc)
  # Break up model paramaters
  cop_parm <- parm[1:ncnr]
  trns_parm <- parm[(ncnr + 1):(ncnr + nr^2 - nr)]
  init_parm <- parm[(ncnr + nr^2 - nr + 1):(ncnr + nr^2 - 1)]
  # Function to optimize log likelihood function
  rsLLH <- function(pp) {
    # Conditional Densities
    # Create paramater index to use in combination with the function factory
    idx <- c(0, cumsum(nc))
    FF <- matrix(ncol = length(x), nrow = nr)
    for (rr in 1:nr) {
      # Function factory here
      FF[rr, ] <- ms_cop_vpdf(family[[rr]], x, y)(pp[(idx[rr] + 1):idx[rr + 1]])
    }
    # Return the negative log likelihood value
    -markov_llh(P, FF, F_0)$logL
  }
  # Turn transition paramaters into markov matrix and initial state vector
  P <- matrix(trns_parm, ncol = nr)
  P <- rbind(P^2, matrix(1, ncol = nr))
  P <- apply(P, 2, function(cc) cc / sum(cc))
  F_0 <- c(init_parm, 1)
  F_0 <- as.matrix(F_0^2 / sum(F_0^2))   
  # Set Boundaries
  bound <- lapply(family, ms_boundary)
  # Fit regime copula paramaters given transition probabilities 
  out <- optim(cop_parm, rsLLH, method = "L-BFGS-B", hessian = TRUE,
               lower = unlist(lapply(bound, `[[`, "LB")),
               upper = unlist(lapply(bound, `[[`, "UB")))
  print(out)
  par <- c(out$par, trns_parm, init_parm)
  return(list(par = par, solver = out, EMstep = "copula"))
}



# For the markov switching model, given the copula parameters, start the
# expectation step of the E-M algorithm by optimizing the markov transition
# matrix
# 
# @param parm A vector of copula parameters.
# 
# @param x A vector of uniform marginal values.
# 
# @param y A vector of uniform marginal values. 
#
# @param family Numeric vector of length one or two indicating the copula type
#   and family.
#
# @return A list of: \describe{
#   \item{par}{Estimated transition parameters and initial state}
#   \item{solver}{Output from the call to \code{\link[stats]{optim}}}
#   \item{EMstep}{Good to know for further development}
#   }
#

markov_optim <- function(parm, x, y, family) {
  
  # Basic constants
  nc <- vapply(family, function(x) switch(as.character(x[1]), "1" = 1, "2" = 2, 3),
               numeric(1))
  nr <- length(family)
  nx <- length(x)
  ncnr <- sum(nc)
  # Break up model paramaters
  cop_parm <- parm[1:ncnr]
  trns_parm <- parm[(ncnr + 1):(ncnr + nr^2 - 1)]
  # Function to optimize
  rsLLH <- function(parm) {
    mark_parm <- parm[1:(nr^2 - nr)]
    init_parm <- parm[(nr^2 - nr + 1):(nr^2 - 1)]
    # Turn transition paramaters into markov matrix and initial state vector
    P <- matrix(mark_parm, ncol = nr)
    P <- rbind(P^2, matrix(1, ncol = nr))
    P <- apply(P, 2, function(cc) cc / sum(cc))
    F_0 <- c(init_parm, 1)
    F_0 <- as.matrix(F_0^2 / sum(F_0^2))  
    # Return the negative log likelihood value
    -markov_llh(P, FF, F_0)$logL
  }
  # Conditional Densities
  # Create paramater index to use in combination with the function factory
  idx <- c(0, cumsum(nc))
  FF <- matrix(ncol = length(x), nrow = nr)
  for (rr in 1:nr) {
    # Function factory here
    FF[rr, ] <- ms_cop_vpdf(family[[rr]], x, y)(cop_parm[(idx[rr] + 1):idx[rr + 1]])
  }
  # Fit regime copula paramaters given transition probabilities 
  out <- optim(trns_parm, rsLLH, method = "BFGS", hessian = TRUE)
  print(out)
  par <- c(cop_parm, out$par)
  return(list(par = par, solver = out, EMstep = "markov"))
}



# For the smooth transition model, calculate the negative log-likelihood value.
# 
# @param x A vector of uniform marginal values.
# 
# @param y A vector of uniform marginal values.
# 
# @param M An nr x np matrix whose columns contain the (np) parameters for each
#   of the nr regimes
# 
# @param family Numeric vector of length one or two indicating the copula type
#   and family.
#
# @return The sum of the negalitve log-likelhood
# 

smooth_llh <- function(x, y, M, family) {
  stopifnot(is.matrix(M))
  if (length(family) == 1) {
    den <- st_cop_vpdf(family, x, y)(M)
  } else {
    wt <- attr(M, "wt")
    den <- wt * st_cop_vpdf(family[[1]], x, y)(M[, 1, drop = FALSE]) +
      (1 - wt) * st_cop_vpdf(family[[2]], x, y)(M[, 2, drop = FALSE])
  }
  -sum(log(den))
}



# For the smooth transition model, produce the smooth transition values of a
# single parameter
# 
# @param parm A parameter vector containing a value for each regime
# 
# @param St A vector the same length as the time series whose value is defined
#   as \code{cumsum(rep(1, N)) \ N}.
#
# @param g A numeric value controlling the pace of transition between states
# 
# @param gsd The standard deviation of \code{St}. Combines with the parameter
#   \code{g} to control the pace of the transition between states
# 
# @param c A numeric value in the unit interval locating the inflection point
#   of the smooth transition
#
# @return A vector with the smooth transition values for a copula parameter
#

smooth_parm <- function(parm, St, g, c, gsd) {
  stopifnot(length(parm) == length(g) + 1,
            length(g) == length(c))
  # gsd <- sd(St)
  logit <- function(S, G, C, GSD) 1 / (1 + exp(-(G / GSD) * (S - C)))
  theta <- matrix(rep(parm[1], length(St)), nrow = 1)
  for (pp in seq_along(g)) {
    theta <- theta +
      (parm[pp + 1] - parm[pp]) * logit(St, g[pp], c[pp], gsd)
  }
  theta
}



# For the smooth transition model, produce the smooth transition values of a
# single parameter
# 
# @param parm A vector containing all the parameters of the smooth transition
#   model
# 
# @param family Numeric vector of length one or two indicating the copula type
#   and family.
#
# @param k An integer indicating the number of regimes
# 
# @param St A vector the same length as the time series whose value is defined
#   as \code{cumsum(rep(1, N)) \ N}.
# 
# @param gsd The standard deviation of \code{St}. Combines with the parameter
#   \code{g} to control the pace of the transition between states
#
# @return An n x np matrix whose columns are the smooth transition values of
#   the copula model. For a mixture copula, the (single) weight parameter is
#   included as an attribute.
#

transform_parm <- function(parm, family, k, St, gsd) {
  # gsd <- sd(St)
  if (length(family) == 1) {
    if (family != 2) {
      mu <- parm[1:k]
      gg <- parm[(k + 1):(2 * k - 1)]
      cc <- parm[(2 * k):(3 * k - 2)]
      TH <- matrix(smooth_parm(mu, St, gg, cc, gsd))
      
    } else if (family == 2) {
      mu <- parm[1:k]
      nu <- parm[(k + 1):(2 * k)]
      gg <- parm[(2 * k + 1):(3 * k - 1)]
      cc <- parm[(3 * k):(4 * k - 2)]
      TH <- matrix(c(smooth_parm(mu, St, gg, cc, gsd),
                     smooth_parm(nu, St, gg, cc, gsd)), ncol = 2)
    }
  }
  if (length(family) == 2) {
    th1 <- parm[1:k]
    th2 <- parm[(k + 1):(2 * k)]
    wt <- parm[2 * k + 1]
    gg <- parm[(2 * k + 2):(3 * k)]
    cc <- parm[(3 * k + 1):(4 * k - 1)]
    TH <- matrix(c(smooth_parm(th1, St, gg, cc, gsd), 
                   smooth_parm(th2, St, gg, cc, gsd)), ncol = 2)
    attr(TH, "wt") <- wt
  }
  TH
}



# Vectorized version of the Clayton copula density
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param alpha Clayton copula parameter
# 
# @return A vector the same length as \code{u} of density values
# 

clayton_pdf <- function(u, v, alpha){
  
  (1 + alpha) * u^-(alpha + 1) * v^-(alpha + 1) *
    (u^-alpha + v^-alpha - 1)^-(1 / alpha + 2)
}



# Vectorized version of the Survival Clayton copula density
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param alpha Clayton copula parameter
# 
# @return A vector the same length as \code{u} of density values
#

Sclayton_pdf <- function(u, v, alpha){
  
  clayton_pdf(1 - u, 1 - v, alpha)
}



# Vectorized version of the Frank copula density
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param delta Frank copula parameter
# 
# @return A vector the same length as \code{u} of density values
#

frank_pdf <- function(u, v, delta) {
  
  (delta * exp(-delta * (u + v)) * (1 - exp(-delta))) /
    ((exp(-delta) - 1) + (exp(-delta * u) - 1) * (exp(-delta * v) - 1))^2
}



# Vectorized version of the Survival Frank copula density
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param delta Frank copula parameter
# 
# @return A vector the same length as \code{u} of density values
#

Sfrank_pdf <- function(u, v, delta) {
 
  frank_pdf(1 - u, 1 - v, delta = delta)
}



# Vectorized version of the Gumbel copula density
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param gamma Gumbel copula parameter
# 
# @return A vector the same length as \code{u} of density values
#

gumbel_pdf <- function(u, v, gamma) {
  
  g1 <- (1 / u) * (1 / v) * (-log(u))^(gamma - 1) * (-log(v))^(gamma - 1) *
    exp(-((-log(u))^gamma + (-log(v))^gamma)^(1 / gamma))
  g2 <- (gamma - 1) * ((-log(u))^gamma + (-log(v))^gamma)^(1 / gamma - 2)
  g3 <- ((-log(u))^gamma + (-log(v))^gamma)^(2 / gamma - 2)
  g1 * (g2 + g3)
}



# Vectorized version of the Survival Gumbel copula density
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param gamma Gumbel copula parameter
# 
# @return A vector the same length as \code{u} of density values
#

Sgumbel_pdf <- function(u, v, gamma) {
  
  gumbel_pdf(1 - u, 1 - v, gamma = gamma)
}



# Vectorized version of the Gaussian copula density
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param rho Correlation parameter
# 
# @return A vector the same length as \code{u} of density values
#

gaussian_pdf <- function(u, v, rho){
  
  stopifnot(abs(rho) <= 1)
  X <- qnorm(u); Y <- qnorm(v)
  numerator <- (rho * X)^2 - (2 * rho * X * Y) + (rho * Y)^2
  sqrt(1 - rho^2) * exp(-numerator / (2 * (1 - rho^2)))
}



# Vectorized version of the Student-t copula density
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param rho Correlation parameter
# 
# @param nu Degrees of Freedom parameter
# 
# @return A vector the same length as \code{u} of density values
#

tcop_pdf <- function(u, v, rho, nu){

  stopifnot(abs(rho) <= 1, nu > 2)
  X <- qt(u, df = nu); Y <- qt(v, df = nu)
  K <- gamma(nu / 2) * gamma((nu + 1) / 2)^-2 * gamma((nu + 2) / 2)
  paren <- (1 + (X^2 + Y^2 - 2 * rho * X * Y) / (nu * (1 - rho^2)))^(-(nu + 2) / 2)
  K * (1 - rho^2)^(-1 / 2) * paren * ((1 + X^2 / nu) * (1 + Y^2 / nu))^((nu + 1) / 2)
}



# For the sequential breakpoint model, calculate a generalized likelihood ratio
# test
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param t An integer indicating where to seperate \code{u} and \code{v} into
#    pre and post regimes
#
# @param fam1 An integer representing the family of the component copula applied
#    to u. This is passed to \code{\link{cop_static}}
#
# @param fam2 An integer representing the family of the component copula applied
#    to v. This is passed to \code{\link{cop_static}}
#
# @return A list of: \enumerate{
#   \item Pre-regime log-likelihood value
#   \item Pre-regime copula parameters
#   \item Post-regime log-likeihood value
#   \item Post-regime copula parameters
#   }          
# 

BreakPointLL <- function(u, v, t, fam1, fam2) {
  
  TT <- length(u)
  stopifnot(TT == length(v))
  LLpre <- cop_static(u[1:t], v[1:t], fam1, fam2)
  LLpost <- cop_static(u[(t + 1):TT], v[(t + 1):TT], fam1, fam2)
  list(LLpre = LLpre$log.likelihood, theta.pre = LLpre$pars,
       LLpost = LLpost$log.likelihood, theta.post = LLpost$pars)
}



# For the sequential breakpoint model, calculate a generalized likelihood ratio
# test
# 
# @param t An integer indicating where to seperate \code{u} and \code{v} into
#    pre and post regimes
#
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param t An integer indexing where to seperate \code{u} and \code{v} into
#    pre and post regimes
#
# @param fam1 An integer representing the family of the component copula applied
#    to u. This is passed to \code{\link{cop_static}}
#
# @param fam2 An integer representing the family of the component copula applied
#    to v. This is passed to \code{\link{cop_static}}
#  
# @param Full A model from \code{\link{cop_static}} providing the full series'
#   log-likelihood value
#
# @return The test statistic for the generalized likelihood ratio test        
# 

BreakPointTestStatistic <- function(t, u, v, fam1, fam2, Full) {
  Break <- BreakPointLL(u, v, t, fam1, fam2)
  2 * (Break$LLpre + Break$LLpost - Full$log.likelihood)
}



# Approximate small sample p-value for the generalized likelihood test
# 
# @param x Test statistic from \code{\link{BreakPointTestStatistic}}
# 
# @param np The number of parameters for the chosen copula model
# 
# @param N The length of the time series
# 
# @return A p-value
# 

aprx <- function(x, np, N) {
  lh <- log(N)^(3/2) / N
  ((x^np * exp(-(x^2 / 2))) / (sqrt(2^np) * gamma(np / 2))) *
    (log(((1 - lh) / lh)^2) - (np / x^2) * log((1 - lh)^2 / lh^2) + 4 / x^2)
}



# For the sequential breakpoint model, search for a breakpoint in a subset of 
# values for two given marginal distributions#' 
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param series An index vector that subsets the full marginal values
#   \code{u, v}
# 
# @param fam1 An integer representing the family of the component copula applied
#    to u. This is passed to \code{\link{cop_static}}
#
# @param fam2 An integer representing the family of the component copula applied
#    to v. This is passed to \code{\link{cop_static}}
#
# @param date_names Optional vector of date names
# 
# @param parallel Logical switch to run the breakpoint search in parallel.
# 
# @param ncores Integer value specifying the number of cores to use
# 
# @param cluster If run in parallel, cluster is passed from a call to
#   \code{\link{BPfit}}
#
# @return A list of: \describe{
#   \item{index}{Index, with resepect to the \strong{full} marginal series, of
#     largest breakpoint test statistic}
#   \item{Date}{If \code{date_names} was provided, the date that corresponds
#     with \code{index}}
#   \item{sqrt(Z)}{The value of the largest breakpoint test statistic}
#   \item{p-value}{The approximated p-value of the breakpoint test statistic}
#   \item{N-size}{The number of observations in \code{series} where the function
#     searched for a breakpoint}
#   }
# 

BreakAnalysis <- function(u, v, series, fam1, fam2 = NULL, date_names = NULL,
                          parallel = FALSE, cluster) {

  stopifnot(!any(missing(u), missing(v), missing(series), missing(fam1)))
  N <- length(series)
  if (is.null(date_names)) {
    date_names <- rep(NA_character_, N)
  }
  # determining the number of parameters can be more robust
  p <- switch(as.character(fam1), "1" = 1, "2" = 2, 3)
  uu <- u[series]
  vv <- v[series]
  low <- ifelse(N < 140, 7, floor(N * 0.05))
  high <- ifelse(N < 140, N - 7, ceiling(N * 0.95))
  Full <- cop_static(uu, vv, fam1, fam2)
  if (parallel && !missing(cluster) & !is.null(cluster)) {
    lamK <- parallel::parSapply(cluster, low:high, FUN = BreakPointTestStatistic,
                                u = uu, v = vv, fam1 = fam1, fam2 = fam2, Full = Full) 
  } else {
    lamK <- sapply(low:high, FUN = BreakPointTestStatistic,
                   u = uu, v = vv, fam1 = fam1, fam2 = fam2, Full = Full)
  }
  plot(lamK)
  maxidx <- which.max(lamK)
  t <- low + maxidx - 1
  Z <- sqrt(lamK[maxidx])
  cc <- (-3:3) + maxidx
  cc <- cc[cc > 0]
  # print(cbind(dat[(min(series) + low:high - 1), c(1, 2)], lamK)[cc, ])
  Results <- list()
  Results[[1]] <- min(series) + t - 1
  Results[[2]] <- date_names[min(series) + t - 1]
  Results[[3]] <- Z
  Results[[4]] <- aprx(Z, p, N)
  Results[[5]] <- N
  names(Results) <- c("index", "Date", "sqrt(Z)", "p-value", "N-Size")
  print(Results)
}



# For the sequential breakpoint model, a closure that recursively searchs for
# all statistically significant breakpoints.
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param f1 An integer representing the family of the component copula applied
#    to u. This is passed to \code{\link{cop_static}}
#
# @param f2 An integer representing the family of the component copula applied
#    to v. This is passed to \code{\link{cop_static}}
# 
# @param parallel Logical switch to run the breakpoint search in parallel.
# 
# @param date_names Optional vector of date names
# 
# @param cluster If run in parallel, cluster is passed from a call to
#   \code{\link{BPfit}}
# 

autoBPtest <- function(x, y, f1, f2 = NULL, parallel = FALSE, date_names = NULL,
                       cluster = NULL) {
  output <- list()
  force(date_names)
  BPtest <- function(srs) {
    bptest <- BreakAnalysis(x, y, srs, f1, f2, date_names, parallel, cluster)
    if (bptest[['p-value']] < 0.05 && bptest[['p-value']] > 0.00 &&
          bptest[['N-Size']] > 20) {
      output[[length(output) + 1]] <<- bptest
      t <- bptest[[1]] # break point
      BPtest(srs = min(srs):t)
      BPtest(srs = (t + 1):max(srs))
    }
  }
  list(BPtest = BPtest,
       value = function() output,
       reset = function() output <<- list()
      )
}



# For the sequential breakpoint model, given a set of breakpoints, run a 
# repartition of the breakpoints for improved consistency
# 
# @param bpResultList The output form \code{\link{autoBPtest}}
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param f1 An integer representing the family of the component copula applied
#    to u. This is passed to \code{\link{cop_static}}
#
# @param f2 An integer representing the family of the component copula applied
#    to v. This is passed to \code{\link{cop_static}}
#    
# @param parallel Logical switch to run the breakpoint search in parallel.
# 
# @param date_names Optional vector of date names
# 
# @param cluster If run in parallel, cluster is passed from a call to
#   \code{\link{BPfit}}
# 
# @reference Repartition method from Bai (1997) and Dias & Embrechts (2009)
# 

repartitionBP <- function(bpResultList, u, v, f1, f2, parallel = F, date_names,
                          cluster = NULL) {
  bpList <- bpResultList
  if (length(bpList) < 2) 
    stop("Less than two break points. Repartition is not needed.")
  # Extract list of break dates and add end points
  idx <- vapply(bpList, `[[`, numeric(1), 1)
  idx <- sort(c(1, idx, length(u)))
  output <- list()
  # Implement repartition method and pass results to 'output'
  for(jj in 1:(length(idx) - 2)) {
    t <- idx[jj]; tt <- idx[jj + 2]
    res <- BreakAnalysis(u, v, t:tt, f1, f2, date_names, parallel, cluster)
    output[[length(output) + 1]] <- res
  }
  output
}



# For the sequential breakpoint model, calculate copula gradients for var-cov
# estimation
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param theta A vector of copula parameters
# 
# @param f1 An integer representing the family of the component copula applied
#    to u. This is passed to \code{\link{cop_static}}
#
# @param f2 An integer representing the family of the component copula applied
#    to v. This is passed to \code{\link{cop_static}}
# 
# @return An n x nr matrix of gradient contributions
# 
# @importFrom VineCopula BiCopPDF
# @importFrom VineCopula BiCopDeriv
# 

cop_gradcontr <- function(u, v, theta, f1, f2) {
  
  if (f1 == 1) {
    rho <- theta[1]
    g1 <- (1 / BiCopPDF(u, v, f1, rho)) * BiCopDeriv(u, v, f1, rho, deriv="par")
    return(matrix(g1, ncol = 1))
  } else if (f1 == 2) {
    rho <- theta[1]; dof <- theta[2]
    g1 <- (1 / BiCopPDF(u, v, f1, rho, dof)) * BiCopDeriv(u, v, f1, rho, dof, "par")
    g2 <- (1 / BiCopPDF(u, v, f1, rho, dof)) * BiCopDeriv(u, v, f1, rho, dof, "par2")
    return(matrix(c(g1, g2), ncol  = 2))
  } else {
    th1 <- theta[1]; th2 <- theta[2]; th3 <- theta[3]
    g1 <- (th3 / cop_pdf(theta, u, v, f1, f2)) * BiCopDeriv(u, v, f1, th1, "par")
    g2 <- ((1 - th3) / cop_pdf(theta, u, v, f1, f2)) * BiCopDeriv(u, v, f2, th2, "par")
    g3 <- (1 / cop_pdf(theta, u, v, f1, f2)) * (BiCopPDF(u, v, f1, th1) -
         BiCopPDF(u, v, f2, th2))
    return(matrix(c(g1, g2, g3), ncol = 3))
  }
}



# For the sequential breakpoint model, calculate copula hessians for var-cov
# estimation
# 
# @param u A vector of uniform marginal values.
# 
# @param v A vector of uniform marginal values.
# 
# @param theta A vector of copula parameters
# 
# @param f1 An integer representing the family of the component copula applied
#    to u. This is passed to \code{\link{cop_static}}
#
# @param f2 An integer representing the family of the component copula applied
#    to v. This is passed to \code{\link{cop_static}}
# 
# @return An nr x nr matrix of gradient contributions
# 
# @importFrom VineCopula BiCopPDF
# @importFrom VineCopula BiCopDeriv
# @importFrom VineCopula BiCopDeriv2
# 

cop_hessian <- function(u, v, theta, f1, f2) {
  
  if (f1 == 1) {
    rho <- theta[1]
    h11 <- (1 / BiCopPDF(u, v, f1, rho)) * BiCopDeriv2(u, v, f1, rho, deriv="par") - 
      (1 / BiCopPDF(u, v, f1, rho)^2) * BiCopDeriv(u, v, f1, rho)^2
    return(matrix(sum(h11)))
  } else if (f1 == 2) {
    rho <- theta[1]; dof <- theta[2]
    h11 <- (1 / BiCopPDF(u, v, f1, rho, dof)) * BiCopDeriv2(u, v, f1, rho, dof, "par") - 
      (1 / BiCopPDF(u, v, f1, rho, dof)^2) * BiCopDeriv(u, v, f1, rho, dof, "par")^2
    h12 <- (1 / BiCopPDF(u, v, f1, rho, dof)) * BiCopDeriv2(u, v, f1, rho, dof, "par1par2") - 
      (1 / BiCopPDF(u, v, f1, rho, dof)^2) * BiCopDeriv(u, v, f1, rho, dof, "par") *
      BiCopDeriv(u, v, f1, rho, dof, "par2")
    h22 <- (1 / BiCopPDF(u, v, f1, rho, dof)) * BiCopDeriv2(u, v, f1, rho, dof, "par2") - 
      (1 / BiCopPDF(u, v, f1, rho, dof)^2) * BiCopDeriv(u, v, f1, rho, dof, "par2")^2
    h11 <- sum(h11); h12 <- sum(h12); h22 <- sum(h22)
    H <- matrix(c(h11, h12, h12, h22), nrow = 2)
    return(H)
  } else {
    th1 <- theta[1]; th2 <- theta[2]; th3 <- theta[3]
    h11 <- (th3 / cop_pdf(theta, u, v, f1, f2)) * BiCopDeriv2(u, v, f1, th1, "par") - 
      ((th3 / cop_pdf(theta, u, v, f1, f2)) * BiCopDeriv(u, v, f1, th1, "par"))^2
    h12 <- -(1 / cop_pdf(theta, u, v, f1, f2)^2) * th3 * (1 - th3) * 
      BiCopDeriv(u, v, f1, th1, "par") * BiCopDeriv(u, v, f2, th2, "par")
    h13 <- BiCopDeriv(u, v, f1, th1) * ((1 / cop_pdf(theta, u, v, f1, f2)) -
          (th3 / cop_pdf(theta, u, v, f1, f2)^2) *
          (BiCopPDF(u, v, f1, th1) - BiCopPDF(u, v, f2, th2)))
    h22 <- ((1 - th3) / cop_pdf(theta, u, v, f1, f2)) *
      BiCopDeriv2(u, v, f2, th2, "par") - (((1 - th3) / cop_pdf(theta, u, v, f1, f2)) *
         BiCopDeriv(u, v, f2, th2, "par"))^2
    h23 <- BiCopDeriv(u, v, f2, th2) * ((1 / cop_pdf(theta, u, v, f1, f2)) -
          ((1 - th3) / cop_pdf(theta, u, v, f1, f2)^2) *
          (BiCopPDF(u, v, f1, th1) - BiCopPDF(u, v, f2, th2)))
    h33 <- -(1 / cop_pdf(theta, u, v, f1, f2)^2) * (BiCopPDF(u, v, f1, th1) -
         BiCopPDF(u, v, f2, th2))^2
    h11 <- sum(h11); h12 <- sum(h12); h13 <- sum(h13)
    h22 <- sum(h22); h23 <- sum(h23); h33 <- sum(h33)
    H <- matrix(c(h11, h12, h13, h12, h22, h23, h13, h23, h33), nrow = 3)
    return(H)
  }
}
