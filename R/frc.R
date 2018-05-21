library(dplyr)
library(midasr)
data("rvsp500", package = "midasr")
# zz <- read.csv("OxfordManRealizedVolatilityIndices.csv", stringsAsFactors = FALSE, skip = 3)

rr <- rvsp500 %>% select(DateID,SPX2.rv) %>% mutate(rv = 100 * sqrt(252*SPX2.rv))

rv1 <- na.omit(rr$rv)

mr1 <- midas_r(rv1 ~ mls(rv1, 1:20,1, harstep), start = list(rv1=c(1,1,1)))

n <- nrow(mr1$model)
wnd <- round(n/3*2)
frc <- NULL

system.time(
  for (i in wnd:n) {
    mrs <- midas_r_simple(y = mr1$model[(i - wnd):i, 1],
                          X = mr1$model[(i - wnd):i, -(1:2)],
                          z = mr1$model[(i - wnd):i, 2],
                          weight = harstep,
                          startx = mr1$start_opt[2:4]
    )
    frc[[i]] <- sum(mrs$midas_coefficients * mrs$model[wnd, -1])
  }
)





s1 <- midas_r_simple(y = mr1$model[, 1], X = mr1$model[, -(1:2)],
                     z = mr1$model[, 2], weight = harstep,
                     startx = mr1$start_opt[2:4])


mrfc <- midasr::forecast(mr1, newdata = list(rv1 = rep(NA, 1)))

forecast(s1, newdata = list(rv1 = rep(NA, 1)))
# Error in ets(object, lambda = lambda, biasadj = biasadj, allow.multiplicative.trend = allow.multiplicative.trend,  : 
# y should be a univariate time series

