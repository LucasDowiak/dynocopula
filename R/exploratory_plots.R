FX <- c("sterling", "swfranc", "yen", "austd", "newz",
        "ffranc", "deutsch", "euro", "safrand", "candd")
# pars
oldpars <- par()
par(mar=c(4,3,2,1) + 0.1,
    mfcol=c(5, 2))

# levels
for (fx in toupper(FX)) {
  plot(DT[, "DATE"], DT[, fx], type="l", main=fx, xlab="", ylab="")
  grid()
}

# returns
for (fx in FX) {
  plot(DT[, "DATE"], DT[, fx], type="l", main=fx, xlab="", ylab="",
       ylim=c(-8,8))
  grid()
}




dtf <- read.csv("data/fx_ffranc_dollar.csv", stringsAsFactors = F)
dtf$X <- NULL
dts <- as.Date(dtf$Date, format="%Y-%m-%d")
dtf$Date <- as.character(format(dts, "%d %b %y"))
write.csv(dtf, file="data/fx_ffranc_dollar.csv", row.names=FALSE)


