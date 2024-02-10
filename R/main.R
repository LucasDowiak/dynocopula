setwd("~/dynocopula/")
library(VineCopula)
source("R/dynofits.R")
source("R/dynohelpers.R")
source("R/marginal_models.R")


epMS <- MSfit(margins$euro, margins$swf, family=list(c(3, 13))) #,
            # initValues=c(-0.5,0.5, 0.6, 0.6, 0.5))

epST <- STfit(margins$euro, margins$swf, family=2, 2)

epBP <- BPfit(margins$yen, margins$swf, fam1=2, parallel=T)


x <- margins$euro
y <- margins$yen

msfit <- MSfit(x, y, family=list(1,1))
msfit <- MSfit(x, y, family=list(2,2))
msfit <- MSfit(x, y, family=list(c(4,14), c(4, 14)))
msfit <- MSfit(x, y, family=list(c(3,13), c(3,13)))

stfit <- STfit(x, y, family=1, 2)
stfit <- STfit(x, y, family=2, 3)
stfit <- STfit(x, y, family=c(3, 13), 3)
stfit <- STfit(x, y, family=c(4, 14), 2)

bpfit <- BPfit(x, y, fam1=1, parallel=T)
bpfit <- BPfit(x, y, fam1=2, parallel=T)
bpfit <- BPfit(x, y, fam1=3, fam2=13, parallel=T)
bpfit <- BPfit(x, y, fam1=4, fam2=14, parallel=T)