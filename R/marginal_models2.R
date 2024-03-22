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

# Lag the euro, sterling, and yen log-returns
DT[, l1_euro := shift(euro)]
DT[, l2_euro := shift(euro, n=2)]
DT[, l1_sterling := shift(sterling)]
DT[, l2_sterling := shift(sterling, n=2)]
DT[, l1_yen := shift(yen)]
DT[, l2_yen := shift(yen, n=2)]

# Lagged Residuals
files_ <- setdiff(list.files("./data/model_objects/"), "tmp1")
collect <- vector("list", length(files_)); names(collect) <- files_
for (f_ in files_) {
  tst <- readRDS(paste(c("./data/model_objects", f_), collapse = "/"))
  curr <- regmatches(f_, regexpr("[a-zA-Z]+", f_))
  dtbl <- data.table(DATE=tst$Dates, currency=curr, resids=as.numeric(rugarch::residuals(tst$model)),
                     res_type="normal_resid")
  dtbl2 <- data.table(DATE=tst$Dates, currency=curr, resids=as.numeric(rugarch::residuals(tst$model, standardize=TRUE)),
                      res_type="standardize_resid")

  collect[[f_]] <- rbindlist(list(dtbl, dtbl2))
}
dtfResids <- rbindlist(collect)
dtfResids <- dcast(rbindlist(collect), DATE ~ currency + res_type, value.var="resids")

dtfResids[, l1_euro_normal_resid := shift(euro_normal_resid)]
dtfResids[, l2_euro_normal_resid := shift(euro_normal_resid, n=2L)]
dtfResids[, l1_sterling_normal_resid := shift(sterling_normal_resid)]
dtfResids[, l2_sterling_normal_resid := shift(sterling_normal_resid, n=2L)]
dtfResids[, l1_yen_normal_resid := shift(yen_normal_resid)]
dtfResids[, l2_yen_normal_resid := shift(yen_normal_resid, n=2L)]
DT <- merge(DT, dtfResids, by="DATE", all.x=TRUE)


# All Currencies
yr80s <- 1980:1989
yr90s <- 1990:1999
yr00s <- 2000:2009
yr10s <- 2010:2018
year_list <- list(yr80s=yr80s, yr90s=yr90s, yr00s=yr00s, yr10s=yr10s)
fxnames <- DT[1:10, names(.SD), .SDcols=which(names(DT) == tolower(names(DT)))]


# Sterling
# -------------------------------------
ster_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
sterling_fit_2000_2009 <- list(model=DT[year(DATE) %in% 2000:2007, ugarchfit(spec=ster_mod, data=sterling)],
                               Dates=DT[year(DATE) %in% 2000:2007, DATE])
res <- verify_marginal_test(marginal_tests(sterling_fit_2000_2009$model), alpha="0.05")
saveRDS(sterling_fit_2000_2009, file="data/model_objects/tmp1/sterling_yr00s.rds")

ster_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(2,2)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std" # "ged"
)
sterling_fit_2010_2018 <- list(model=DT[year(DATE) %in% 2008:2018, ugarchfit(spec=ster_mod, data=sterling)],
                               Dates=DT[year(DATE) %in% 2008:2018, DATE])
res <- verify_marginal_test(marginal_tests(sterling_fit_2010_2018$model), alpha="0.05")
saveRDS(sterling_fit_2010_2018, file="data/model_objects/tmp1/sterling_yr10s.rds")


# -------------------------------------


# Japanese Yen
# -------------------------------------
yen_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
yen_fit_2000_2009 <- list(model=DT[year(DATE) %in% 2000:2007, ugarchfit(spec=yen_mod, data=yen)],
                          Dates=DT[year(DATE) %in% 2000:2007, DATE])
res <- verify_marginal_test(marginal_tests(yen_fit_2000_2009$model), alpha="0.05")
saveRDS(yen_fit_2000_2009, file="data/model_objects/tmp1/yen_yr00s.rds")


yen_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
yen_fit_2010_2018 <- list(model=DT[year(DATE) %in% 2008:2018, ugarchfit(spec=yen_mod, data=yen)],
                          Dates=DT[year(DATE) %in% 2008:2018, DATE])
res <- verify_marginal_test(marginal_tests(yen_fit_2010_2018$model), alpha="0.05")
saveRDS(yen_fit_2010_2018, file="data/model_objects/tmp1/yen_yr10s.rds")


# -------------------------------------


# Euro
# -------------------------------------
euro_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(2,1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
euro_fit_2000_2009 <- list(model=DT[year(DATE) %in% 2000:2007, ugarchfit(spec=euro_mod, data=euro)],
                           Dates=DT[year(DATE) %in% 2000:2007, DATE])
res <- verify_marginal_test(marginal_tests(euro_fit_2000_2009$model), alpha="0.05")
saveRDS(euro_fit_2000_2009, file="data/model_objects/tmp1/euro_yr00s.rds")


euro_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)
euro_fit_2010_2018 <- list(model=DT[year(DATE) %in% 2008:2018, ugarchfit(spec=euro_mod, data=euro)],
                           Dates=DT[year(DATE) %in% 2008:2018, DATE])
res <- verify_marginal_test(marginal_tests(euro_fit_2010_2018$model), alpha="0.05")
saveRDS(euro_fit_2010_2018, file="data/model_objects/tmp1/euro_yr10s.rds")



# -------------------------------------


# Marginal Model Regression Summaries
# -----------------------------------------------------------------------------
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
  cat("----------------------------------------------------------------\n")
  print(sprintf("Object: %s", f_))
  print(tst$model)
  cat("\n----------------------------------------------------------------\n\nn")
}
dtfU <- dcast(rbindlist(collect_u), DATE ~ currency, value.var="u")
dtfEval <- rbindlist(collect_e)
dtfCoef <- rbindlist(collect_i)

mod <- readRDS("data/model_objects/yen_yr00s.RDS")
coef(mod$model)

# 
vars <- c('mu', "ar1",
          'omega', 'alpha1', 'gamma1', 'alpha2', 'gamma2','beta1', 'beta2',
          'shape')
dtfVar <- data.table(variable=factor(vars, levels = vars))
D1 <- dcast(dtfCoef, group + vcv_type + variable + currency ~ .,
            value.var=c("Estimate"))
setnames(D1, ".", "Estimate")
D2 <- dcast(dtfCoef, group + vcv_type + variable + currency ~ .,
            value.var=c("Std. Error"))
setnames(D2, ".", "Std. Error")


# Regression Summaries
D <- merge(D1, D2)
D[, `Std. Error`:=round(`Std. Error`, 3)]
D[, Estimate := round(Estimate, 3)]
D[, variable := factor(variable, levels=vars)]

D2 <- D[group=="00" & vcv_type=="standard" & currency=="euro"]
melt(merge(dtfVar, D2, all.x = TRUE),
     id.vars=c("currency", "variable", "group"), variable.name = "Stat",
           measure.vars = c("Estimate", "Std. Error"))[order(variable, Stat)]

# ----------------------------------------------------


# Validation Summaries
# ----------------------------------------------------
tmp <- dcast(dtfEval[!grepl("NYBLOM", Test)],
                     Test + Spec + group ~ currency,
                     value.var="P-Value")[order(group, Test, Spec)]
tmp[group=="00"] # , round(.SD, 4), .SDcols=4:6]
tmp[group=="10", round(.SD, 3), .SDcols=4:6]

tmp <- dcast(dtfEval[grepl("NYBLOM", Test)],
             Test + Spec + group ~ currency,
             value.var="Stat")[order(group, Test, Spec)]
tmp[group=="00"] # , round(.SD, 4), .SDcols=4:6]
tmp[group=="00", round(.SD, 3), .SDcols=4:6]

# ----------------------------------------------------

# Orthogonality Tests for Conditional Copula
# ------------------------------------------------------------------------------

############ Mean model

# 2000 - 2009
# Not needed since no AR or MA components in the mean equation


# 2010 - 2018
summary(lm(euro_normal_resid ~ l1_sterling, data=DT[year(DATE) %in% 2010:2018]))
summary(lm(yen_normal_resid ~ l1_sterling, data=DT[year(DATE) %in% 2010:2018]))



############ Variance

# 2000 - 2009
summary(lm(I(euro_standardize_resid**2) ~ I(l1_sterling_normal_resid**2),
           data=DT[year(DATE) %in% 2000:2009]))
summary(lm(I(euro_standardize_resid**2) ~ I(l1_yen_normal_resid**2),
           data=DT[year(DATE) %in% 2000:2009]))


summary(lm(I(sterling_standardize_resid**2) ~ I(l1_euro_normal_resid**2) + I(l2_euro_normal_resid**2),
           data=DT[year(DATE) %in% 2000:2009]))
summary(lm(I(sterling_standardize_resid**2) ~ I(l1_yen_normal_resid**2),
           data=DT[year(DATE) %in% 2000:2009]))


summary(lm(I(yen_standardize_resid**2) ~ I(l1_euro_normal_resid**2) + I(l2_euro_normal_resid**2),
           data=DT[year(DATE) %in% 2000:2009]))
summary(lm(I(yen_standardize_resid**2) ~ I(l1_sterling_normal_resid**2),
           data=DT[year(DATE) %in% 2000:2009]))


# 2010 - 2018
summary(lm(I(euro_standardize_resid**2) ~ I(l1_sterling_normal_resid**2) + I(l2_sterling_normal_resid**2),
           data=DT[year(DATE) %in% 2010:2018]))
summary(lm(I(euro_standardize_resid**2) ~ I(l1_yen_normal_resid**2),
           data=DT[year(DATE) %in% 2010:2018]))

summary(lm(I(sterling_standardize_resid**2) ~ I(l1_euro_normal_resid**2),
           data=DT[year(DATE) %in% 2010:2018]))
summary(lm(I(sterling_standardize_resid**2) ~ I(l1_yen_normal_resid**2),
           data=DT[year(DATE) %in% 2010:2018]))


summary(lm(I(yen_standardize_resid**2) ~ I(l1_euro_normal_resid**2),
           data=DT[year(DATE) %in% 2010:2018]))
summary(lm(I(yen_standardize_resid**2) ~ I(l1_sterling_normal_resid**2) + I(l2_sterling_normal_resid**2),
           data=DT[year(DATE) %in% 2010:2018]))


# ------------------------------------------------------------------------------



