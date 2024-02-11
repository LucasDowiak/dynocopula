setwd("~/Git/dynocopula/")
library(VineCopula)
library(data.table)
source("R/dynofits.R")
source("R/dynohelpers.R")



# DTU is the data.frame that contains the iid uniform marginals
files_ <- list.files("./data/model_objects/")
collect_u <- vector("list", length(files_)); names(collect_u) <- files_
collect_e <- collect_u
for (f_ in files_) {
  tst <- readRDS(paste(c("./data/model_objects", f_), collapse = "/"))
  curr <- regmatches(f_, regexpr("[a-zA-Z]+", f_))
  collect_u[[f_]] <- data.table(DATE=tst$Dates, currency=curr, u=as.numeric(rugarch::pit(tst$mod)))
  
  decade <- regmatches(f_, regexpr("[0-9]{2}", f_))
  MT <- verify_marginal_test(marginal_tests(tst$model), alpha = "0.05")
  MT[, `:=`(currency=curr, group=decade)]
  collect_e[[f_]] <- MT
}
dtfU <- dcast(rbindlist(collect_u), DATE ~ currency, value.var="u")
dtfEval <- rbindlist(collect_e)


BPfit(dtfU[DATE %in% year_list[[4]], yen], dtfU[DATE %in% year_list[[4]], sterling],)

DTU_hist <- DTU
DTU <- DTU_hist
DTU <- data.table(DTU)
DTU[, index := 1:.N]




# currencies
fxnames1 <- c("austd", "deutsch", "sterling", "yen")
fxnames2 <- c("austd", "euro", "sterling", "yen")

fetch_union <- function(x, y, Dt)
{
  stopifnot(all(c(x, y, "DATE") %in% names(Dt),
                inherits(Dt, "data.table")))
  bool <- complete.cases(Dt[, .SD, .SDcols=c(x, y)])
  
  # return empty list if overlap of currencies has fewer than 2 years
  if (sum(bool) < 760) {
    return(list())
  } else {
    return(list(x=Dt[bool, eval(parse(text=x))],
                y=Dt[bool, eval(parse(text=y))],
                date_range=Dt[bool, range(DATE)]))
  }
}


# 2H
dt_names <- dtfU[year(DATE) %in% unlist(year_list[3:4])][["DATE"]]
for (pair in combn(fxnames2, 2, simplify=FALSE)) {
  print(sprintf("Pair %s-%s started at %s", pair[1], pair[2], Sys.time()))
  aa <- BPfit(dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[1]]],
              dtfU[year(DATE) %in% unlist(year_list[3:4])][[pair[2]]], 
              date_names=as.character(dt_names), fam1=2, parallel = T, ncores = 4)
  saveRDS(aa, file=sprintf("data/dynocop_objects/%s_%s_tcop_00_18.RDS", pair[1], pair[2]))
}



for (fx1 in fxnames) {

    ii <- which(fx1 == fxnames)
  
    if (ii < length(fxnames)) {
  
        for (fx2 in setdiff(fxnames[-c(1:ii)], fx1)) {
            lst <- fetch_union(fx1, fx2, DTU)
        
            if (length(lst) > 0) {
              
                print(sprintf("%s and %s started at %s.", fx1, fx2, Sys.time()))
                assign(paste(fx1, fx2, "bp", sep="_"),
                       BPfit(lst$x, lst$y, fam1=2, parallel=T, ncores=4,
                             date_names=seq.Date(lst$date_range[1], lst$date_range[2], by="day")))
          
            }
            # Take a moment to cool off after parallel process
            Sys.sleep(60)
        }
    }
}


ls(pattern="*_bp$")
cmd <- paste0("list(", paste(ls(pattern="*_bp$"), collapse=", "), ")")
assign("BP_tdist", eval(parse(text=cmd)))
names(BP_tdist) <- ls(pattern="*_bp$")
saveRDS(BP_tdist, file="data/BP_tdist.RDS")



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

