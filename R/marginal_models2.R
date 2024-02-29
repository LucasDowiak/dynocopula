setwd("~/Git/dynocopula")
library("quantmod")
library("xts")
library("zoo")
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

zz <- DT[year(DATE) %in% 2000:2009, auto_fit(austd, max_arch = 2, max_garch = 2, ignore_nyblom = TRUE)]

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
ster_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
sterling_fit_2000_2009 <- DT[year(DATE) %in% 2000:2009, ugarchfit(spec=ster_mod, data=sterling)]
res <- verify_marginal_test(marginal_tests(sterling_fit_2000_2009), alpha="0.05")
saveRDS(sterling_fit_2000_2009, file="data/model_objects/sterling_yr00s.rds")

ster_mod <- ugarchspec(
  variance.model = list(model = "csGARCH", garchOrder = c(2,2)),
  mean.model = list(armaOrder = c(0, 2), include.mean = FALSE),
  distribution.model = "std" # "ged"
)
sterling_fit_2010_2018 <- DT[year(DATE) %in% 2010:2018, ugarchfit(spec=ster_mod, data=sterling)]
res <- verify_marginal_test(marginal_tests(sterling_fit_2010_2018), alpha="0.05")
saveRDS(sterling_fit_2010_2018, file="data/model_objects/sterling_yr10s.rds")

# -------------------------------------


# Japanese Yen
# -------------------------------------
yen_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
yen_fit_2000_2009 <- DT[year(DATE) %in% 2000:2009, ugarchfit(spec=yen_mod, data=yen)]
res <- verify_marginal_test(marginal_tests(yen_fit_2000_2009), alpha="0.05")
saveRDS(yen_fit_2000_2009, file="data/model_objects/yen_yr00s.rds")

yen_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
yen_fit_2010_2018 <- DT[year(DATE) %in% 2010:2018, ugarchfit(spec=yen_mod, data=yen)]
res <- verify_marginal_test(marginal_tests(yen_fit_2010_2018), alpha="0.05")
saveRDS(yen_fit_2010_2018, file="data/model_objects/yen_yr10s.rds")
# -------------------------------------


# Euro
# -------------------------------------

euro_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(2,1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
euro_fit_2000_2009 <- DT[year(DATE) %in% 2000:2009, ugarchfit(spec=euro_mod, data=euro)]
res <- verify_marginal_test(marginal_tests(euro_fit_2000_2009), alpha="0.05")
saveRDS(euro_fit_2000_2009, file="data/model_objects/euro_yr00s.rds")


euro_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)
euro_fit_2010_2018 <- DT[year(DATE) %in% 2010:2018, ugarchfit(spec=euro_mod, data=euro)]
res <- verify_marginal_test(marginal_tests(euro_fit_2010_2018), alpha="0.05")
saveRDS(euro_fit_2010_2018, file="data/model_objects/euro_yr10s.rds")

# -------------------------------------


