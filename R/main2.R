setwd("~/Git/dynocopula/")
library(VineCopula)
library(cubature)
library(data.table)
library(rugarch)
source("R/dynofits.R")
source("R/dynohelpers.R")
source("R/normality_tests.R")


# DTU is the data.frame that contains the iid uniform marginals
files_ <- list.files("./data/model_objects/")
files_ <- files_[-5]
collect_u <- vector("list", length(files_)); names(collect_u) <- files_
collect_e <- collect_u
for (f_ in files_) {
  tst <- readRDS(paste(c("./data/model_objects", f_), collapse = "/"))
  curr <- regmatches(f_, regexpr("[a-zA-Z]+", f_))
  collect_u[[f_]] <- data.table(DATE=tst$Dates, currency=curr, u=as.numeric(rugarch::pit(tst$model)))
  
  decade <- regmatches(f_, regexpr("[0-9]{2}", f_))
  MT <- verify_marginal_test(marginal_tests(tst$model), alpha = "0.05", ignore_nyblom=TRUE)
  MT[, `:=`(currency=curr, group=decade)]
  collect_e[[f_]] <- MT
}
dtfU <- dcast(rbindlist(collect_u), DATE ~ currency, value.var="u")
dtfEval <- rbindlist(collect_e)


# currencies
fxnames <- c("euro", "sterling", "yen")

year_list <- list(1980:1989, 1990:1999, 2000:2009, 2010:2018)


# Single Regime Data Estimation
# ------------------------------------------------------------------------------
uu <- dtfU[, sterling]
vv <- dtfU[, yen]
summary(BiCopEst(uu, vv, family=4))


optim(par=c(1.15, 1.06, 0.5), cop_llh, u=uu, v=vv, fam1=4, fam2=14,
      method="L-BFGS-B", lower=c(1, 1, 0), upper=c(17, 17, 0.999))


# BIC
3 * log(dtfU[, .N]) - 2 * 236.3893# 163.1549

# ------------------------------------------------------------------------------


# Break Point Optimization
# ------------------------------------------------------------------------------
# 2H T-cop
dt_names <- dtfU[year(DATE) %in% unlist(year_list[3:4])][["DATE"]]
for (pair in combn(fxnames2, 2, simplify=FALSE)[1:2]) {
  print(sprintf("Pair %s-%s started at %s", pair[1], pair[2], Sys.time()))
  aa <- try(BPfit(dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[1]]],
                  dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[2]]], 
                  date_names=as.character(dt_names), fam1=2, parallel = TRUE, ncores = 6))
  if (inherits(aa, "try-error"))
    next
  saveRDS(aa, file=sprintf("data/dynocop_objects/%s_%s_tcop_00_18.RDS", pair[1], pair[2]))
}

# 2H Clayton
dt_names <- dtfU[year(DATE) %in% unlist(year_list[3:4])][["DATE"]]
for (pair in combn(fxnames2, 2, simplify=FALSE)) {
  print(sprintf("Pair %s-%s started at %s", pair[1], pair[2], Sys.time()))
  aa <- try(BPfit(dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[1]]],
              dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[2]]], 
              date_names=as.character(dt_names), fam1=3, fam2=13, parallel = T, ncores = 6))
  if (inherits(aa, "try-error"))
    next
  saveRDS(aa, file=sprintf("data/dynocop_objects/%s_%s_clayton_00_18.RDS", pair[1], pair[2]))
}

# 2H Gumbel
dt_names <- dtfU[year(DATE) %in% unlist(year_list[3:4])][["DATE"]]
for (pair in combn(fxnames2, 2, simplify=FALSE)) {
  print(sprintf("Pair %s-%s started at %s", pair[1], pair[2], Sys.time()))
  aa <- try(BPfit(dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[1]]],
              dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[2]]], 
              date_names=as.character(dt_names), fam1=4, fam2=14, parallel = T, ncores = 6))
  if (inherits(aa, "try-error"))
    next
  saveRDS(aa, file=sprintf("data/dynocop_objects/%s_%s_gumbel_00_18.RDS", pair[1], pair[2]))
}

# ------------------------------------------------------------------------------


# Markov Switching
# ------------------------------------------------------------------------------
dt_names <- dtfU[, DATE]
fxnames <- c("euro", "sterling", "yen")
states <- 2:3
tol <- 1e-5

for (s in states) {
  for (pair in combn(fxnames, 2, simplify=FALSE)) {
    print(sprintf("Pair: %s-%s Regimes: %d started at %s", pair[1], pair[2], s, Sys.time()))
    aa <- MSfit(dtfU[[pair[1]]],
                dtfU[[pair[2]]], 
                family=as.list(rep(2, s)), tol=tol)
    saveRDS(aa, file=sprintf("data/dynocop_objects/MSfit_%s_%s_tcop_00_18_s%d.RDS", pair[1], pair[2], s))
  }
}


for (s in states) {
  for (pair in combn(fxnames, 2, simplify=FALSE)) {
    print(sprintf("Pair: %s-%s Regimes: %d started at %s", pair[1], pair[2], s, Sys.time()))
    aa <- MSfit(dtfU[[pair[1]]],
                dtfU[[pair[2]]], 
                family=split(rep(c(3, 13), s), floor((1:(2*s) / 2) - 0.5)), tol=tol)
    saveRDS(aa, file=sprintf("data/dynocop_objects/MSfit_%s_%s_clayton_00_18_s%d.RDS", pair[1], pair[2], s))
  }
}


for (s in states) {
  for (pair in combn(fxnames, 2, simplify=FALSE)) {
    print(sprintf("Pair: %s-%s Regimes: %d started at %s", pair[1], pair[2], s, Sys.time()))
    aa <- MSfit(dtfU[[pair[1]]],
                dtfU[[pair[2]]], 
                family=split(rep(c(4, 14), s), floor((1:(2*s) / 2) - 0.5)), tol=tol)
    saveRDS(aa, file=sprintf("data/dynocop_objects/MSfit_%s_%s_gumbel_00_18_s%d.RDS", pair[1], pair[2], s))
  }
}



aa <- MSfit(dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[1]]],
            dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[2]]], 
            family=list(2,2,2))
# ------------------------------------------------------------------------------


# Smooth Transition
# ------------------------------------------------------------------------------
dt_names <- dtfU[, DATE]
fxnames <- c("euro", "sterling", "yen")
states <- 4
tol <- 1e-5

for (s in states) {
  for (pair in combn(fxnames, 2, simplify=FALSE)) {
    print(sprintf("Pair: %s-%s Regimes: %d started at %s", pair[1], pair[2], s, Sys.time()))
    aa <- STfit(dtfU[[pair[1]]],
                dtfU[[pair[2]]], 
                family=2, regimes=s)
    saveRDS(aa, file=sprintf("data/dynocop_objects/STfit_%s_%s_tcop_00_18_s%d.RDS", pair[1], pair[2], s))
  }
}


for (s in states) {
  for (pair in combn(fxnames, 2, simplify=FALSE)) {
    print(sprintf("Pair: %s-%s Regimes: %d started at %s", pair[1], pair[2], s, Sys.time()))
    aa <- STfit(dtfU[[pair[1]]],
                dtfU[[pair[2]]], 
                family=c(3, 13), regimes=s)
    saveRDS(aa, file=sprintf("data/dynocop_objects/STfit_%s_%s_clayton_00_18_s%d.RDS", pair[1], pair[2], s))
  }
}


for (s in states) {
  for (pair in combn(fxnames, 2, simplify=FALSE)) {
    print(sprintf("Pair: %s-%s Regimes: %d started at %s", pair[1], pair[2], s, Sys.time()))
    aa <- STfit(dtfU[[pair[1]]],
                dtfU[[pair[2]]], 
                family=c(4, 14), regimes=s)
    saveRDS(aa, file=sprintf("data/dynocop_objects/STfit_%s_%s_gumbel_00_18_s%d.RDS", pair[1], pair[2], s))
  }
}




# ------------------------------------------------------------------------------


# Model Selection
# ----------------------------------------------------

files_ <- list.files("data/dynocop_objects/", pattern="BPfit_euro_sterling")
lst_mods <- lapply(files_, function(x) readRDS(sprintf("data/dynocop_objects/%s", x)))

lapply(lst_mods, length)
Ls <- lapply(lst_mods, function(x) sum(logLik(x)))
ks <- lapply(lst_mods, function(x) length(x) * length(x[[1]]$pars))
N <- dtfU[, .N]
data.table(file=files_, BIC=mapply(function(l, k) log(N) * k - 2 * (l), Ls, ks))

summary(lst_mods[[3]])




files_ <- list.files("data/dynocop_objects/", pattern="STfit_sterling_yen")
lst_mods <- lapply(files_, function(x) readRDS(sprintf("data/dynocop_objects/%s", x)))
aicbic <- aic_bic(lst_mods)
aicbic$BIC[, colnames(aicbic$BIC)[grepl("4", colnames(aicbic$BIC))]]



files_ <- list.files("data/dynocop_objects/", pattern="MSfit_sterling_yen")
lst_mods <- lapply(files_, function(x) readRDS(sprintf("data/dynocop_objects/%s", x)))

Ls <- lapply(lst_mods, function(x) logLik(x))
ks <- lapply(lst_mods, function(x) length(coef(x)))
N <- dtfU[, .N]
data.table(file=files_, BIC=mapply(function(l, k) log(N) * k - 2 * (l), Ls, ks))


summary(readRDS("data/dynocop_objects/STfit_sterling_yen_tcop_00_18_s3.RDS"))

# ----------------------------------------------------

# Date of Regime Switches
# ----------------------------------------------------
mod <- readRDS("data/dynocop_objects/MSfit_sterling_yen_tcop_00_18_s2.RDS")
P <- mod$regime.inference$smooth.prob
idx <- which(abs(diff(P[1, ] > 0.5)) == 1) + 1
dtfU[c(1, idx), .(DATE)]


mod <- readRDS("data/dynocop_objects/STfit_sterling_yen_tcop_00_18_s3.RDS")
p <- coef(mod)
dtfU[round(.N * p[grepl("c", names(p))])]


mod <- readRDS("data/dynocop_objects/BPfit_sterling_yen_tcop_00_18.RDS")
sapply(rev(attr(mod, "initial_bp")), `[[`, "Date")





# ----------------------------------------------------


# Plot Daily Parameter Values
# ----------------------------------------------------

prep_tcop_values <- function(obj)
{
  cls <- class(obj)
  if (cls == "smoothTransCopula") {
    rho <- obj$smooth.parameters[, 1]
    nu <- obj$smooth.parameters[, 2]
    ktaus <- BiCopPar2Tau(2, rho, nu)
    tdeps <- BiCopPar2TailDep(2, rho, nu)$upper # upper tail dep
  } else if (cls == "seqBreakPoint") {
    idx <- unique(as.vector(vapply(obj, `[[`, double(2), "points")))
    idx2 <- seq(length(diff(idx))) %% 2 != 0
    ktau <- vapply(obj, function(x) x$dep.measures$Ktau, double(1))
    tdep <- vapply(obj, function(x) unlist(x$dep.measures$tail_dep), double(2))
    ktaus <- rep(ktau, diff(idx)[idx2] + 1)
    tdeps <- rep(tdep["upper", ], diff(idx)[idx2] + 1) # upper tail dep
  } else if (cls == "markovCopula") {
    sp <- t(obj$regime.inference$smooth.prob)
    ktau <- unlist(obj$copula$ktau)
    tdep <- sapply(obj$copula$tail_dep, `[[`, 2) # upper tail dep
    ktaus <- rowSums(sweep(sp, 2, ktau, `*`))
    tdeps <- rowSums(sweep(sp, 2, tdep, `*`))
  } else {
    stop(sprintf("Function does not support obj of class '%s'", cls))
  }
  return(data.frame(ktau=ktaus, tdep=tdeps))
}

f1 <- c("BPfit_euro_sterling_tcop_00_18.RDS",
        "BPfit_euro_yen_tcop_00_18.RDS",
        "BPfit_sterling_yen_tcop_00_18.RDS",
        "MSfit_euro_sterling_tcop_00_18_s2.RDS",
        "MSfit_euro_yen_tcop_00_18_s2.RDS",
        "MSfit_sterling_yen_tcop_00_18_s2.RDS",
        "STfit_euro_sterling_tcop_00_18_s4.RDS",
        "STfit_euro_yen_tcop_00_18_s4.RDS",
        "STfit_sterling_yen_tcop_00_18_s3.RDS")

collect <- vector("list", length(f1)); names(collect) <- f1
for (ii in seq_along(collect)) {
  file_ <- f1[ii]
  file_loc <- sprintf("data/dynocop_objects/%s", file_)
  obj_ <- readRDS(file_loc)
  
  tmp_dtf <- prep_tcop_values(obj_)
  fnms <- strsplit(file_, split="_")[[1]]
  tmp_dtf$pair <- paste(fnms[2:3], collapse = "-")
  tmp_dtf$model <- fnms[1]
  tmp_dtf$Date <- DT[year(DATE) %in% 2000:2018, DATE]
  collect[[file_]] <- tmp_dtf
}
dtfPlot <- rbindlist(collect)
dtfPlot[, pair := factor(pair)]
dtfPlot[, model := factor(model)]

# Plots
# ---
par(mfrow=c(2, 1), mar=c(3, 3, 3, 2) + 0.1)
kylim <- c(-0.2, 0.7)
tylim <- c(-0.05, 0.6)
mcex <- 0.85
# Seq Break Point
dtfPlot[pair=="euro-sterling" & model=="BPfit",
        plot(Date, ktau, type="l",
             xlab="", ylab="", main="Sequential Break Point: Kendall's Tau",
             ylim=kylim,
             lty=1,
             cex.main=mcex)]
dtfPlot[pair=="euro-yen" & model=="BPfit", lines(Date, ktau, lty=2)]
dtfPlot[pair=="sterling-yen" & model=="BPfit", lines(Date, ktau, lty=3)]
grid()
legend("topright", legend = c("Euro-Pound", "Euro-Yen", "Pound-Yen"),
       cex=0.5, lty=c(1, 2, 3), bty="n")

dtfPlot[pair=="euro-sterling" & model=="BPfit",
        plot(Date, tdep, type="l",
             xlab="", ylab="", main="Sequential Break Point: Tail Dependence",
             ylim=tylim,
             lty=1,
             cex.main=mcex)]
dtfPlot[pair=="euro-yen" & model=="BPfit", lines(Date, tdep, lty=2)]
dtfPlot[pair=="sterling-yen" & model=="BPfit", lines(Date, tdep, lty=3)]
grid()
legend("topright", legend = c("Euro-Pound", "Euro-Yen", "Pound-Yen"),
       cex=0.5, lty=c(1, 2, 3), bty="n")



# Markov Switching
dtfPlot[pair=="euro-sterling" & model=="MSfit",
        plot(Date, ktau, type="l",
             xlab="", ylab="", main="Markov Transition: Kendall's Tau",
             ylim=kylim,
             lty=1,
             cex.main=mcex)]
dtfPlot[pair=="euro-yen" & model=="MSfit", lines(Date, ktau, lty=2)]
dtfPlot[pair=="sterling-yen" & model=="MSfit", lines(Date, ktau, lty=3)]
grid()
legend("topright", legend = c("Euro-Pound", "Euro-Yen", "Pound-Yen"),
       cex=0.5, lty=c(1, 2, 3), bty="n")

dtfPlot[pair=="euro-sterling" & model=="MSfit",
        plot(Date, tdep, type="l",
             xlab="", ylab="", main="Markov Transition: Tail Dependence",
             ylim=tylim,
             lty=1,
             cex.main=mcex)]
dtfPlot[pair=="euro-yen" & model=="MSfit", lines(Date, tdep, lty=2)]
dtfPlot[pair=="sterling-yen" & model=="MSfit", lines(Date, tdep, lty=3)]
grid()
legend("topright", legend = c("Euro-Pound", "Euro-Yen", "Pound-Yen"),
       cex=0.5, lty=c(1, 2, 3), bty="n")


# Smooth Transition
dtfPlot[pair=="euro-sterling" & model=="STfit",
        plot(Date, ktau, type="l",
             xlab="", ylab="", main="Smooth Transition: Kendall's Tau",
             ylim=kylim,
             lty=1,
             cex.main=mcex)]
dtfPlot[pair=="euro-yen" & model=="STfit", lines(Date, ktau, lty=2)]
dtfPlot[pair=="sterling-yen" & model=="STfit", lines(Date, ktau, lty=3)]
grid()
legend("topright", legend = c("Euro-Pound", "Euro-Yen", "Pound-Yen"),
       cex=0.5, lty=c(1, 2, 3), bty="n")

dtfPlot[pair=="euro-sterling" & model=="STfit",
        plot(Date, tdep, type="l",
             xlab="", ylab="", main="Smooth Transition: Tail Dependence",
             ylim=tylim,
             lty=1,
             cex.main=mcex)]
dtfPlot[pair=="euro-yen" & model=="STfit", lines(Date, tdep, lty=2)]
dtfPlot[pair=="sterling-yen" & model=="STfit", lines(Date, tdep, lty=3)]
grid()
legend("topright", legend = c("Euro-Pound", "Euro-Yen", "Pound-Yen"),
       cex=0.5, lty=c(1, 2, 3), bty="n")



# Rolling Kendall' Tau
# ------------------------------------------------------------------------------
taus <- c()
w <- 60
for (rr in (w + 1):(nrow(dtfU) - w)) {
  tau <- dtfU[rr:(rr + w - 1), cor(sterling, yen, method="kendall")]
  taus <- c(taus, tau)
}





# ------------------------------------------------------------------------------