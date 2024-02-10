setwd("~/dynocopula/")
library(VineCopula)
library(data.table)
source("R/dynofits.R")
source("R/dynohelpers.R")
source("R/marginal_models2.R")


# DTU is the data.frame that contains the iid uniform marginals
DTU_hist <- DTU
DTU <- DTU_hist
DTU <- data.table(DTU)
DTU[, index := 1:.N]



# currencies
fxnames <- c("austd", "deutsch", "euro", "ffranc", "newz",
             "sterling", "swfranc", "won", "yen")
fxnames <- fxnames[-1]

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

