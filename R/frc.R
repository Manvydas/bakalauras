# Accuracy functions --------------------------------------------------
## |-- MSE (Mean Standart Error) 
mse <- function(Yt, Ft) {
  mean((Yt - Ft) ^ 2)
}
## |-- MAPE (Mean Average Percentage Error)
mape <- function(Yt, Ft) {
  (1 / length(Yt)) * (sum(abs((Yt - Ft) / Yt))) * 100
}
## |-- MASE (Mean Absolute Scaled Error)
mase <- function(Yt, Ft) {
  nn <- length(Yt)
  sum(abs(Yt - Ft)) / ((nn / (nn - 1)) * (sum(abs(Yt[2:nn] - Yt[1:nn - 1]))))
}


library(dplyr)
library(midasr)
data("rvsp500", package = "midasr")

rr <- rvsp500 %>% select(DateID,SPX2.rv) %>% mutate(rv = 100 * sqrt(252*SPX2.rv))

rv1 <- na.omit(rr$rv)

mr1 <- midas_r(rv1 ~ mls(rv1, 1:20,1, harstep), start = list(rv1=c(1,1,1)))

n <- nrow(mr1$model)
wnd <- round(n/3*2)
frc <- NULL

nn <- length(rv1)
clmn <- ncol(mr1$model)
yy <- c(mr1$model[, 1], NA)
intr <- c(mr1$model[, 2], mean(mr1$model[, 2]))
dd <- rbind(mr1$model[, 3:clmn], rv1[nn:(nn-19)])

data <- cbind(yy, dd, intr)


system.time(
  for (i in wnd:n) {
    mrs <- midas_r_simple(y = data[(i - wnd):i, 1],
                          X = data[(i - wnd):i, -c(1,clmn)],
                          z = data[(i - wnd):i, clmn],
                          weight = harstep,
                          startx = mr1$start_opt[2:4]
    )
    frc[[i+1]] <- sum(mrs$midas_coefficients * data[i+1, -1])
  }
)


dataaa <- data.frame("real" = yy[wnd:n],
                    "pred" = frc[wnd:n]
)

data.frame("mse" = mse(dataa[, 1], dataa[, 2]),
           "mape" = mape(dataa[, 1], dataa[, 2]),
           "mase" = mase(dataa[, 1], dataa[, 2])
)







#-------------------------------------------------------------------

s1 <- midas_r_simple(y = mr1$model[, 1], X = mr1$model[, -(1:2)],
                     z = mr1$model[, 2], weight = harstep,
                     startx = mr1$start_opt[2:4])


mrfc <- midasr::forecast(mr1, newdata = list(rv1 = rep(NA, 1)))

forecast(s1, newdata = list(rv1 = rep(NA, 1)))
# Error in ets(object, lambda = lambda, biasadj = biasadj, allow.multiplicative.trend = allow.multiplicative.trend,  : 
# y should be a univariate time series

