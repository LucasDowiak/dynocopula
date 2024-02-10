setwd("~/dynocopula")
source("~/conditional_copula/R/normality_tests.R")
library("rugarch")

# bring in data
dtf <- read.csv("data/foreign_exchange_20180602.csv", skip=9,
                  colClasses=c("character", rep("numeric", 7)),
                  col.names=c("date", "Aust_d", "Won", "Euro", "Pound",
                              "Yen", "Newz_d", "SWfranc"))
dtf$date <- as.Date(dtf$date, format = "%d %b %Y")

# turn spot rates into returns log(X_t) - log(X_t-1)
fxreturns <- lapply(dtf[, -1], function(x) c(NA, 100 * diff(log(x))))
names(fxreturns) <- tolower(names(fxreturns))
dtf <- cbind(dtf, as.data.frame(fxreturns))
rm(fxreturns)

subbool <- c(F, rep(TRUE, nrow(dtf) - 1)) #  dtf$date >= "1999-01-04" & dtf$date <= "2012-12-31"

# Euro
euro_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
euro_fit <- ugarchfit(
  spec = euro_mod,
  data = dtf$euro[subbool]
)

# Pound marginal specification
pound_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model= list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
pound_fit <- ugarchfit(
  spec = pound_mod,
  data = dtf$pound[subbool]
)

# Yen marginal specification
yen_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model= list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
yen_fit <- ugarchfit(
  spec = yen_mod,
  data = dtf$yen[subbool]
)

# Australian dollar marginal specification
"JB test isn't very good"
aust_mod <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
  mean.model= list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
aust_fit <- ugarchfit(
  spec = aust_mod,
  data = dtf$aust_d[subbool]
)

# Yen marginal specification
won_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model= list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "std"
)
won_fit <- ugarchfit(
  spec = won_mod,
  data = dtf$won[subbool]
)

# Yen marginal specification
"JB is not so good here either"
newz_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model= list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
newz_fit <- ugarchfit(
  spec = newz_mod,
  data = dtf$newz_d[subbool]
)

# Yen marginal specification
swf_mod <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(3,3)),
  mean.model= list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "sstd"
)
swf_fit <- ugarchfit(
  spec = swf_mod,
  data = dtf$swfranc[subbool]
)

margins <- list(
  euro = PIT(euro_fit),
  pound = PIT(pound_fit),
  yen = PIT(yen_fit),
  won = PIT(won_fit),
  aust = PIT(aust_fit),
  newz = PIT(newz_fit),
  swf = PIT(swf_fit)
)
