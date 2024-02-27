setwd("~/Git/dynocopula")
library("quantmod")
source("R/auto_marginal.R")

# bring in data
# won comes from St Louis FRED site
# the rest comes from BOE site
files_ <- dir("data/", pattern="^fx_", full.names=T)
csvs <- vector("list", length(files_)); names(csvs) <- files_
for (f_ in files_) {
  fx <- strsplit(f_, "_")[[1]][[2]]
  DT <- read.csv(f_, colClasses="character", header=F, skip=1,
                 col.names=c("Date", fx), na.strings="\\.")
  DT[, fx] <- as.numeric(DT[, fx])
  DT$Date <- as.Date(DT$Date, format = "%d %b %y")
  csvs[[f_]] <- DT
}
DT <- Reduce(function(x,y) merge(x, y, by="Date", all=TRUE), csvs)
names(DT) <- toupper(names(DT))

# turn spot rates into returns log(X_t) - log(X_t-1)
fxreturns <- lapply(DT[, -1], function(x) c(NA_real_, diff(log(x))))
names(fxreturns) <- tolower(names(fxreturns))
DT <- as.data.table(cbind(DT, as.data.frame(fxreturns)))
rm(fxreturns)

# All Currencies
yr80s <- 1980:1989
yr90s <- 1990:1999
yr00s <- 2000:2009
yr10s <- 2010:2018
year_list <- list(yr80s=yr80s, yr90s=yr90s, yr00s=yr00s, yr10s=yr10s)
fxnames <- DT[1:10, names(.SD), .SDcols=which(names(DT) == tolower(names(DT)))]


# Full History Currencies
# -------------------------------------
for (yrs in names(year_list)) {
  x <- DT[year(DATE) %in% year_list[[yrs]],  sterling]
  dts <- DT[year(DATE) %in% year_list[[yrs]], DATE]
  mod <- auto_fit(x, max_arch=3, max_garch=3, ignore_nyblom=TRUE)
  mod[["Dates"]] <- dts
  saveRDS(mod, file=sprintf("data/model_objects/sterling_%s.RDS", yrs))
}

for (yrs in names(year_list)) {
  x <- DT[year(DATE) %in% year_list[[yrs]],  austd]
  dts <- DT[year(DATE) %in% year_list[[yrs]], DATE]
  mod <- auto_fit(x,  max_arch=3, max_garch=3, ignore_nyblom=TRUE)
  mod[["Dates"]] <- dts
  saveRDS(mod, file=sprintf("data/model_objects/austd_%s.RDS", yrs))
}

for (yrs in names(year_list)) {
  x <- DT[year(DATE) %in% year_list[[yrs]],  candd]
  dts <- DT[year(DATE) %in% year_list[[yrs]], DATE]
  mod <- auto_fit(x)
  mod[["Dates"]] <- dts
  saveRDS(mod, file=sprintf("data/model_objects/candd_%s.RDS", yrs))
}

for (yrs in names(year_list)) {
  x <- DT[year(DATE) %in% year_list[[yrs]],  newz]
  dts <- DT[year(DATE) %in% year_list[[yrs]], DATE]
  mod <- auto_fit(x)
  mod[["Dates"]] <- dts
  saveRDS(mod, file=sprintf("data/model_objects/newz_%s.RDS", yrs))
}

for (yrs in names(year_list)) {
  x <- DT[year(DATE) %in% year_list[[yrs]],  yen]
  dts <- DT[year(DATE) %in% year_list[[yrs]], DATE]
  mod <- auto_fit(x,  max_arch=3, max_garch=3, ignore_nyblom=TRUE)
  mod[["Dates"]] <- dts
  saveRDS(mod, file=sprintf("data/model_objects/yen_%s.RDS", yrs))
}

for (yrs in names(year_list)) {
  x <- DT[year(DATE) %in% year_list[[yrs]],  swfranc]
  dts <- DT[year(DATE) %in% year_list[[yrs]], DATE]
  mod <- auto_fit(x)
  mod[["Dates"]] <- dts
  saveRDS(mod, file=sprintf("data/model_objects/swfranc_%s.RDS", yrs))
}

# -------------------------------------



# Pre-Euro Currencies
# -------------------------------------
for (yrs in names(year_list)[1:2]) {
  x <- DT[year(DATE) %in% year_list[[yrs]],  deutsch]
  dts <- DT[year(DATE) %in% year_list[[yrs]], DATE]
  mod <- auto_fit(x,  max_arch=3, max_garch=3, ignore_nyblom=TRUE)
  mod[["Dates"]] <- dts
  saveRDS(mod, file=sprintf("data/model_objects/deutsch_%s.RDS", yrs))
}

for (yrs in names(year_list)[1:2]) {
  x <- DT[year(DATE) %in% year_list[[yrs]],  ffranc]
  dts <- DT[year(DATE) %in% year_list[[yrs]], DATE]
  mod <- auto_fit(x)
  mod[["Dates"]] <- dts
  saveRDS(mod, file=sprintf("data/model_objects/ffranc_%s.RDS", yrs))
}
# -------------------------------------


# Euro and Sing D
# -------------------------------------

# Euro
for (yrs in names(year_list)[3:4]) {
  x <- DT[year(DATE) %in% year_list[[yrs]],  euro]
  dts <- DT[year(DATE) %in% year_list[[yrs]], DATE]
  mod <- auto_fit(x,  max_arch=3, max_garch=3, ignore_nyblom=TRUE)
  mod[["Dates"]] <- dts
  saveRDS(mod, file=sprintf("data/model_objects/euro_%s.RDS", yrs))
}

# Singapore Dollar - Special case with a first date value of NA 
for (yrs in names(year_list)[3:4]) {
  x <- DT[year(DATE) %in% year_list[[yrs]] & !is.na(singd),  singd]
  dts <- DT[year(DATE) %in% year_list[[yrs]] & !is.na(singd), DATE]
  mod <- auto_fit(x)
  mod[["Dates"]] <- dts
  saveRDS(mod, file=sprintf("data/model_objects/singd_%s.RDS", yrs))
}
# -------------------------------------

find_coef_matrix <- function(obj)
{
  m <- as.data.table(obj@fit$matcoef)
  setnames(m, names(m), trimws(names(m)))
  m[["variable"]] <- rownames(obj@fit$matcoef)
  m[["vcv_type"]] <- "standard"
  rm <- as.data.table(obj@fit$robust.matcoef)
  setnames(rm, names(rm), trimws(names(rm)))
  rm[["variable"]] <- rownames(obj@fit$robust.matcoef)
  rm[["vcv_type"]] <- "robust"
  return(rbindlist(list(m, rm)))
}


files_ <- list.files("./data/model_objects/")
collect_u <- vector("list", length(files_)); names(collect_u) <- files_
collect_i <- collect_e <- collect_u
for (f_ in files_) {
  tst <- readRDS(paste(c("./data/model_objects", f_), collapse = "/"))
  curr <- regmatches(f_, regexpr("[a-zA-Z]+", f_))
  collect_u[[f_]] <- data.table(DATE=tst$Dates, currency=curr, u=as.numeric(rugarch::pit(tst$mod)))
  
  decade <- regmatches(f_, regexpr("[0-9]{2}", f_))
  MT <- verify_marginal_test(marginal_tests(tst$model), alpha = "0.05")
  MT[, `:=`(currency=curr, group=decade)]
  collect_e[[f_]] <- MT
  
  CM <- find_coef_matrix(tst$model)
  CM[, `:=`(currency=curr, group=decade)]
  collect_i[[f_]] <- CM
  # collect_i[[f_]][, `:=`(currency=curr, group=decade))]
}
dtfU <- dcast(rbindlist(collect_u), DATE ~ currency, value.var="u")
dtfEval <- rbindlist(collect_e)
dtfCoef <- rbindlist(collect_i)

mod <- readRDS("data/model_objects/austd_yr10s.RDS")
mod$model
dtfEval[currency=="sterling" & group=="00"]

zz <- DT[year(DATE) %in% 2000:2009, auto_fit(austd, max_arch = 3, max_garch = 3, ignore_nyblom = TRUE)]

vars <- c('mu', 'ar1', 'ar2', 'ma1', 'ma2', 'ma3', 'omega', 'alpha1', 'gamma1',
          'alpha2', 'gamma2', 'beta1', 'beta2', 'beta3', 'eta11', 'eta21', 'shape',
          'skew')
dtfVar <- data.table(variable=factor(vars, levels = vars))
D1 <- dcast(dtfCoef, group + vcv_type + variable + currency ~ .,
            value.var=c("Estimate"))
setnames(D1, ".", "Estimate")
D2 <- dcast(dtfCoef, group + vcv_type + variable + currency ~ .,
            value.var=c("Std. Error"))
setnames(D2, ".", "Std. Error")

D <- merge(D1, D2)
D[, `Std. Error`:=round(`Std. Error`, 3)]
D[, Estimate := round(Estimate, 3)]
D[, variable := factor(variable, levels=vars)]
DD <- merge(dtfVar,
            D[order(variable)][group=="10" & vcv_type=="robust" & currency=="austd"],
      all.x=TRUE)
melt(DD, id.vars = c("variable"), measure.vars=c("Estimate", "Std. Error"),
     variable.name="numeric")[order(variable, numeric)]

# Sterling
# -------------------------------------
subbool <- DT$DATE >= "1980-01-01" & DT$DATE <= "1989-12-31"
mod <- auto_fit(DT[subbool, sterling])
ster_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
  distribution.model = "sstd"
)
sterling_fit_1980_1989 <- ugarchfit(
  spec = ster_mod,
  data = DT$sterling[subbool]
)
marginal_tests(sterling_fit_1980_1989)
DTU <- merge(DTU,
             data.frame(DATE=DT$DATE[subbool], sterling=PIT(sterling_fit_1980_1989)),
             by="DATE", all.x=T)


subbool <- DT$DATE >= "1990-01-01" & DT$DATE <= "1999-12-31"
ster_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
  distribution.model = "std"
)
sterling_fit_1990_1999 <- ugarchfit(
  spec = ster_mod,
  data = DT$sterling[subbool]
)
marginal_tests(sterling_fit_1990_1999)
DTU[subbool, "sterling"] <- PIT(sterling_fit_1990_1999)


subbool <- DT$DATE >= "2000-01-01" & DT$DATE <= "2009-12-31"
ster_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
  distribution.model = "ged"
)
sterling_fit_2000_2009 <- ugarchfit(
  spec = ster_mod,
  data = DT$sterling[subbool]
)
marginal_tests(sterling_fit_2000_2009)
DTU[subbool, "sterling"] <- PIT(sterling_fit_2000_2009)


subbool <- DT$DATE >= "2010-01-01" & DT$DATE <= "2018-12-31"
ster_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
  distribution.model = "std"
)
sterling_fit_2010_2018 <- ugarchfit(
  spec = ster_mod,
  data = DT$sterling[subbool]
)
marginal_tests(sterling_fit_2010_2018)
DTU[subbool, "sterling"] <- PIT(sterling_fit_2010_2018)

# -------------------------------------


# Swiss Franc
# -------------------------------------
# TODO : Hong-Li portmantau fails
subbool <- DT$DATE >= "1980-01-01" & DT$DATE <= "1989-12-31"
swf_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(3,1)),
  mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
  distribution.model = "std"
)
swfranc_1980_1989 <- ugarchfit(
  spec = swf_mod,
  data = DT$swfranc[subbool]
)
marginal_tests(swfranc_1980_1989)
DTU <- merge(DTU,
             data.frame(DATE=DT$DATE[subbool], swfranc=PIT(swfranc_1980_1989)),
             by="DATE", all.x=T)


subbool <- DT$DATE >= "1990-01-01" & DT$DATE <= "1999-12-31"
swf_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,1), include.mean = TRUE),
  distribution.model = "sstd"
)
swfranc_1990_1999 <- ugarchfit(
  spec = swf_mod,
  data = DT$swfranc[subbool]
)
marginal_tests(swfranc_1990_1999)
DTU[subbool, "swfranc"] <- PIT(swfranc_1990_1999)


subbool <- DT$DATE >= "2000-01-01" & DT$DATE <= "2009-12-31"
swf_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
  distribution.model = "sstd"
)
swfranc_2000_2009 <- ugarchfit(
  spec = swf_mod,
  data = DT$swfranc[subbool]
)
marginal_tests(swfranc_2000_2009)
DTU[subbool, "swfranc"] <- PIT(swfranc_2000_2009)


# TODO: Hong-Li portmantaeu test fails
subbool <- DT$DATE >= "2010-01-01" & DT$DATE <= "2018-12-31"
swf_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
swfranc_2010_2018 <- ugarchfit(
  spec = swf_mod,
  data = DT$swfranc[subbool]
)
marginal_tests(swfranc_2010_2018)
DTU[subbool, "swfranc"] <- PIT(swfranc_2010_2018)


# -------------------------------------


# Japanese Yen
# -------------------------------------
# TODO : Hong-Li portmantau fails
subbool <- DT$DATE >= "1980-01-01" & DT$DATE <= "1989-12-31"
yen_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(5,0), include.mean = TRUE),
  distribution.model = "sstd"
)
yen_1980_1989 <- ugarchfit(
  spec = yen_mod,
  data = DT$yen[subbool]
)
marginal_tests(yen_1980_1989)
DTU <- merge(DTU,
             data.frame(DATE=DT$DATE[subbool], yen=PIT(yen_1980_1989)),
             by="DATE", all.x=T)


subbool <- DT$DATE >= "1990-01-01" & DT$DATE <= "1999-12-31"
yen_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
yen_1990_1999 <- ugarchfit(
  spec = yen_mod,
  data = DT$yen[subbool]
)
marginal_tests(yen_1990_1999)
DTU[subbool, "yen"] <- PIT(yen_1990_1999)


# TODO: LM 2nd and 4th moment fails
subbool <- DT$DATE >= "2000-01-01" & DT$DATE <= "2009-12-31"
yen_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(3,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
yen_2000_2009 <- ugarchfit(
  spec = yen_mod,
  data = DT$yen[subbool]
)
marginal_tests(yen_2000_2009)
DTU[subbool, "yen"] <- PIT(yen_2000_2009)



subbool <- DT$DATE >= "2010-01-01" & DT$DATE <= "2018-12-31"
yen_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
yen_2010_2018 <- ugarchfit(
  spec = yen_mod,
  data = DT$yen[subbool]
)
marginal_tests(yen_2010_2018)
DTU[subbool, "yen"] <- PIT(yen_2010_2018)

# -------------------------------------


# Australian Dollar
# -------------------------------------
# TODO: Hong-Li 3rd and 4th Moment fails 
subbool <- DT$DATE >= "1980-01-01" & DT$DATE <= "1989-12-31"
austd_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(3,1), include.mean = TRUE),
  distribution.model = "sstd"
)
austd_1980_1989 <- ugarchfit(
  spec = austd_mod,
  data = DT$austd[subbool]
)
marginal_tests(austd_1980_1989)
DTU <- merge(DTU,
             data.frame(DATE=DT$DATE[subbool], austd=PIT(austd_1980_1989)),
             by="DATE", all.x=T)


subbool <- DT$DATE >= "1990-01-01" & DT$DATE <= "1999-12-31"
austd_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
austd_1990_1999 <- ugarchfit(
  spec = austd_mod,
  data = DT$austd[subbool]
)
marginal_tests(austd_1990_1999)
DTU[subbool, "austd"] <- PIT(austd_1990_1999)


subbool <- DT$DATE >= "2000-01-01" & DT$DATE <= "2009-12-31"
austd_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
austd_2000_2009 <- ugarchfit(
  spec = austd_mod,
  data = DT$austd[subbool]
)
marginal_tests(austd_2000_2009)
DTU[subbool, "austd"] <- PIT(austd_2000_2009)


subbool <- DT$DATE >= "2010-01-01" & DT$DATE <= "2018-12-31"
austd_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
austd_2010_2018 <- ugarchfit(
  spec = austd_mod,
  data = DT$austd[subbool]
)
marginal_tests(austd_2010_2018)
DTU[subbool, "austd"] <- PIT(austd_2010_2018)

# -------------------------------------


# New Zealand Dollar
# -------------------------------------
# TODO: Not Evaluated
subbool <- DT$DATE >= "1980-01-01" & DT$DATE <= "1989-12-31"
newz_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(3,1)),
  mean.model = list(armaOrder = c(5,0), include.mean = TRUE),
  distribution.model = "sstd"
)
newz_1980_1989 <- ugarchfit(
  spec = newz_mod,
  data = DT$newz[subbool]
)
marginal_tests(newz_1980_1989)
DTU <- merge(DTU,
             data.frame(DATE=DT$DATE[subbool], newz=PIT(newz_1980_1989)),
             by="DATE", all.x=T)



# TODO: Hong-Li 3rd and 4th moment fail
subbool <- DT$DATE >= "1990-01-01" & DT$DATE <= "1999-12-31"
newz_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(5,1)),
  mean.model = list(armaOrder = c(5,0), include.mean = TRUE),
  distribution.model = "sstd"
)
newz_1990_1999 <- ugarchfit(
  spec = newz_mod,
  data = DT$newz[subbool]
)
marginal_tests(newz_1990_1999)
DTU[subbool, "newz"] <- PIT(newz_1990_1999)


# TODO: Hong-Li portmantaeu fails & Shapiro-Wilks fails at 10 %
subbool <- DT$DATE >= "2000-01-01" & DT$DATE <= "2009-12-31"
newz_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(5,1)),
  mean.model = list(armaOrder = c(5,0), include.mean = FALSE),
  distribution.model = "sstd"
)
newz_2000_2009 <- ugarchfit(
  spec = newz_mod,
  data = DT$newz[subbool]
)
marginal_tests(newz_2000_2009)
DTU[subbool, "newz"] <- PIT(newz_2000_2009)


subbool <- DT$DATE >= "2010-01-01" & DT$DATE <= "2018-12-31"
newz_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,1), include.mean = TRUE),
  distribution.model = "sstd"
)
newz_2010_2018 <- ugarchfit(
  spec = newz_mod,
  data = DT$newz[subbool]
)
marginal_tests(newz_2010_2018)
DTU[subbool, "newz"] <- PIT(newz_2010_2018)

# -------------------------------------


# French Franc
# -------------------------------------
subbool <- DT$DATE >= "1980-01-01" & DT$DATE <= "1989-12-31"
ff_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(2,1)),
  mean.model = list(armaOrder = c(2,1), include.mean = TRUE),
  distribution.model = "std"
)
ffranc_1980_1989 <- ugarchfit(
  spec = ff_mod,
  data = DT$ffranc[subbool]
)
marginal_tests(ffranc_1980_1989)
DTU <- merge(DTU,
             data.frame(DATE=DT$DATE[subbool], ffranc=PIT(ffranc_1980_1989)),
             by="DATE", all.x=T)



# TODO: LM 4th moment fails
subbool <- DT$DATE >= "1990-01-01" & DT$DATE <= "1999-12-31"
ff_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(5,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
ffranc_1990_1999 <- ugarchfit(
  spec = ff_mod,
  data = DT$ffranc[subbool]
)
marginal_tests(ffranc_1990_1999)
DTU[subbool, "ffranc"] <- PIT(ffranc_1990_1999)


subbool <- DT$DATE >= "2000-01-01" & DT$DATE <= "2001-12-31"
ff_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(3,1)),
  mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
  distribution.model = "sstd"
)
ffranc_2000_2009 <- ugarchfit(
  spec = ff_mod,
  data = DT$ffranc[subbool]
)
marginal_tests(ffranc_2000_2009)
DTU[subbool, "ffranc"] <- PIT(ffranc_2000_2009)

# -------------------------------------


# German Deutschemark
# -------------------------------------
subbool <- DT$DATE >= "1980-01-01" & DT$DATE <= "1989-12-31"
dmk_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(3,1)),
  mean.model = list(armaOrder = c(1,1), include.mean = TRUE),
  distribution.model = "std"
)
deutsch_1980_1989 <- ugarchfit(
  spec = dmk_mod,
  data = DT$deutsch[subbool]
)
marginal_tests(deutsch_1980_1989)
DTU <- merge(DTU,
             data.frame(DATE=DT$DATE[subbool], deutsch=PIT(deutsch_1980_1989)),
             by="DATE", all.x=T)


subbool <- DT$DATE >= "1990-01-01" & DT$DATE <= "1999-12-31"
dmk_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(3,1)),
  mean.model = list(armaOrder = c(1,1), include.mean = TRUE),
  distribution.model = "std"
)
deutsch_1990_1999 <- ugarchfit(
  spec = dmk_mod,
  data = DT$deutsch[subbool]
)
marginal_tests(deutsch_1990_1999)
DTU[subbool, "deutsch"] <- PIT(deutsch_1990_1999)


# TODO: LM 2nd and 4th moment fails
subbool <- DT$DATE >= "2000-01-01" & DT$DATE <= "2001-12-31"
dmk_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(1,5), include.mean = TRUE),
  distribution.model = "sged"
)
deutsch_2000_2009 <- ugarchfit(
  spec = dmk_mod,
  data = DT$deutsch[subbool]
)
marginal_tests(deutsch_2000_2009)
DTU[subbool, "deutsch"] <- PIT(deutsch_2000_2009)

# -------------------------------------


# Euro
# -------------------------------------
# TODO: LM 2nd moment fails
subbool <- DT$DATE >= "1999-01-05" & DT$DATE <= "2009-12-31"
euro_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(5,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
euro_2000_2009 <- ugarchfit(
  spec = euro_mod,
  data = DT$euro[subbool]
)
marginal_tests(euro_2000_2009)
DTU <- merge(DTU,
             data.frame(DATE=DT$DATE[subbool], euro=PIT(euro_2000_2009)),
             by="DATE", all.x=T)


subbool <- DT$DATE >= "2010-01-01" & DT$DATE <= "2018-12-31"
euro_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
euro_2010_2018 <- ugarchfit(
  spec = euro_mod,
  data = DT$euro[subbool]
)
marginal_tests(euro_2010_2018)
DTU[subbool, "euro"] <- PIT(euro_2010_2018)
# -------------------------------------


# Singapore Dollar
# -------------------------------------
subbool <- DT$DATE >= "2000-01-05" & DT$DATE <= "2009-12-31"
sing_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(5,1), include.mean = TRUE),
  distribution.model = "sstd"
)
singd_2000_2009 <- ugarchfit(
  spec = sing_mod,
  data = DT$singd[subbool]
)
marginal_tests(singd_2000_2009)
DTU <- merge(DTU,
             data.frame(DATE=DT$DATE[subbool], singd=PIT(singd_2000_2009)),
             by="DATE", all.x=T)


subbool <- DT$DATE >= "2010-01-01" & DT$DATE <= "2018-12-31"
sing_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
singd_2010_2018 <- ugarchfit(
  spec = sing_mod,
  data = DT$singd[subbool]
)
marginal_tests(singd_2010_2018)
DTU[subbool, "singd"] <- PIT(singd_2010_2018)

# -------------------------------------


# South Korea Won
# -------------------------------------
subbool <- DT$DATE >= "2005-04-04" & DT$DATE <= "2009-12-31"
won_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
won_2000_2009 <- ugarchfit(
  spec = won_mod,
  data = DT$won[subbool]
)
marginal_tests(won_2000_2009)
DTU <- merge(DTU,
             data.frame(DATE=DT$DATE[subbool], won=PIT(won_2000_2009)),
             by="DATE", all.x=T)


subbool <- DT$DATE >= "2010-01-01" & DT$DATE <= "2018-12-31"
won_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
won_2010_2018 <- ugarchfit(
  spec = won_mod,
  data = DT$won[subbool]
)
marginal_tests(won_2010_2018)
DTU[subbool, "won"] <- PIT(won_2010_2018)
# -------------------------------------
