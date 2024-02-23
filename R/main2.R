setwd("~/Git/dynocopula/")
library(VineCopula)
library(data.table)
library(rugarch)
source("R/dynofits.R")
source("R/dynohelpers.R")
source("R/normality_tests.R")


# DTU is the data.frame that contains the iid uniform marginals
files_ <- list.files("./data/model_objects/")
collect_u <- vector("list", length(files_)); names(collect_u) <- files_
collect_e <- collect_u
for (f_ in files_) {
  tst <- readRDS(paste(c("./data/model_objects", f_), collapse = "/"))
  curr <- regmatches(f_, regexpr("[a-zA-Z]+", f_))
  collect_u[[f_]] <- data.table(DATE=tst$Dates, currency=curr, u=as.numeric(rugarch::pit(tst$model)))
  
  decade <- regmatches(f_, regexpr("[0-9]{2}", f_))
  MT <- verify_marginal_test(marginal_tests(tst$model), alpha = "0.05")
  MT[, `:=`(currency=curr, group=decade)]
  collect_e[[f_]] <- MT
}
dtfU <- dcast(rbindlist(collect_u), DATE ~ currency, value.var="u")
dtfEval <- rbindlist(collect_e)


# currencies
fxnames1 <- c("austd", "deutsch", "sterling", "yen")
fxnames2 <- c("austd", "euro", "sterling", "yen")

year_list <- list(1980:1989, 1990:1999, 2000:2009, 2010:2018)
# Break Point Optimization
# ------------------------------------------------------------------------------
# 2H T-cop
dt_names <- dtfU[year(DATE) %in% unlist(year_list[3:4])][["DATE"]]
for (pair in combn(fxnames2, 2, simplify=FALSE)) {
  print(sprintf("Pair %s-%s started at %s", pair[1], pair[2], Sys.time()))
  aa <- BPfit(dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[1]]],
              dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[2]]], 
              date_names=as.character(dt_names), fam1=2, parallel = T, ncores = 4)
  saveRDS(aa, file=sprintf("data/dynocop_objects/%s_%s_tcop_00_18.RDS", pair[1], pair[2]))
}

# 2H Clayton
dt_names <- dtfU[year(DATE) %in% unlist(year_list[3:4])][["DATE"]]
for (pair in combn(fxnames2, 2, simplify=FALSE)) {
  print(sprintf("Pair %s-%s started at %s", pair[1], pair[2], Sys.time()))
  aa <- BPfit(dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[1]]],
              dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[2]]], 
              date_names=as.character(dt_names), fam1=3, fam2=13, parallel = T, ncores = 6)
  saveRDS(aa, file=sprintf("data/dynocop_objects/%s_%s_clayton_00_18.RDS", pair[1], pair[2]))
}

# 2H Gumbel
dt_names <- dtfU[year(DATE) %in% unlist(year_list[3:4])][["DATE"]]
for (pair in combn(fxnames2, 2, simplify=FALSE)) {
  print(sprintf("Pair %s-%s started at %s", pair[1], pair[2], Sys.time()))
  aa <- BPfit(dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[1]]],
              dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[2]]], 
              date_names=as.character(dt_names), fam1=4, fam2=14, parallel = T, ncores = 6)
  saveRDS(aa, file=sprintf("data/dynocop_objects/%s_%s_gumbel_00_18.RDS", pair[1], pair[2]))
}

# ------------------------------------------------------------------------------


# Markov Switching
# ------------------------------------------------------------------------------
dt_names <- dtfU[year(DATE) %in% unlist(year_list[3:4])][["DATE"]]
for (pair in combn(fxnames2, 2, simplify=FALSE)) {
  print(sprintf("Pair %s-%s started at %s", pair[1], pair[2], Sys.time()))
  aa <- MSfit(dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[1]]],
              dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[2]]], 
              date_names=as.character(dt_names), fam1=4, fam2=14, parallel = T, ncores = 6)
  saveRDS(aa, file=sprintf("data/dynocop_objects/MSfit_%s_%s_gumbel_00_18.RDS", pair[1], pair[2]))
}


aa <- MSfit(dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[1]]],
            dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[2]]], 
            family=list(2,2,2))
# ------------------------------------------------------------------------------


plot_country_pairs <- function(cty, modlist, Dt)
{
  require(wesanderson)
  
  modlist <- modlist[grepl(cty, names(modlist))]
  stopifnot(length(modlist) > 0)
  dts <- Reduce(range, lapply(modlist, attr, "date_names"))
  dts <- Dt[DATE >= dts[1] & DATE <= dts[2], DATE]
  
  bool <- !is.na(Dt[, eval(parse(text=cty))])
  
  pairctrys <- lapply(strsplit(names(modlist), "_"), setdiff, c(cty, "bp"))
  
  f_ <- function(obj)
  {
    idx <- unique(as.vector(vapply(obj, `[[`, double(2), "points")))
    idx2 <- seq(length(diff(idx))) %% 2 != 0
    ktau <- vapply(obj, function(x) x$dep.measures$Ktau, double(1))
    rep(ktau, diff(idx)[idx2] + 1)
  }
  
  colors_ <- wes_palette("Darjeeling1", length(modlist), "continuous")
  
  par(bg="gray40")
  plot(dts, rep(0, length(dts)), type="n", xlab="Date", ylab="Ktau")
  title(cty)
  grid()
  legend("bottomleft", legend=pairctrys, lty=1, lwd=2, col=colors_, bg="gray40")
  
  for (cc in seq_along(modlist)) {
    bool2 <- !is.na(Dt[, eval(parse(text=pairctrys[[cc]]))])
    dts2 <- Dt[bool & bool2, DATE]
    ktau <- f_(modlist[[cc]])
    lines(dts2, ktau, col=colors_[[cc]], lwd=2)
  }
}


for (fx in fxnames) {
  plot_country_pairs(fx, BP_gaussian, DTU)
}







lst <- fetch_union("sterling", "euro", DTU)

bpfit <- BPfit(lst$x, lst$y, fam1=2, parallel=T, ncores=4, date_names=seq.Date(lst$date_range[1], lst$date_range[2], by="day"))
msfit <- MSfit(lst$x, lst$y, family=list(1,1,1,1))
msfit3 <- MSfit(lst$x, lst$y, family=list(1,1,1))
stfit <- STfit(lst$x, lst$y, family=1, regimes=4)

