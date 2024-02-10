nms <- names(dtf)[1:8]

dtf2 <- dtf[, .SD, .SDcols = nms]

base <- dtf2[year(date) == 2013, colMeans(.SD), .SDcols=nms[-1]]

dtf3 <- dtf2[, as.data.table(sweep(.SD, 2, base, `/`)), .SDcols=nms[-1]]

dtf2 <- cbind(dtf2[, date], dtf3)

for (n in nms[-1]) {
  if (n == "Aust_d") {
    dtf2[, plot(V1, eval(parse(text=n)), type="l")]
  } else {
    lines(dtf2[, V1], dtf2[,eval(parse(text=n))])
  }
}
grid()
