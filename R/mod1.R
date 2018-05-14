library(data.table)
library(midasr)
setwd("~/Documents/bakalauras/R")
dat <- fread("OxfordManRealizedVolatilityIndices.csv", select = 1:2, header = F)[-c(1:2), ]

setnames(dat, as.character(dat[1,]))
dat <- as.data.frame(dat[-1,])
dat[, 2] <- as.numeric(dat[, 2])
dat[is.na(dat)] <- NA
dat[dat == ""] <- NA
dat <- na.omit(dat)

spx2_rvol2 <- 100 * sqrt(252 * dat[, "SPX2.rv"])

mh2 <- midas_r(rv ~ mls(rv, 1:20, 1, harstep), data = list(rv = spx2_rvol2), start = list(rv = c(1, 1, 1)))
summary(mh2)

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

## |-- MAE (Mean Average Error)
mae <- function(Yt, Ft) {
  (1 / length(Yt)) * sum(abs(Yt - Ft))
}
## |-- RMSE (Root Mean Squared Error)
rmse <- function(Yt, Ft) {
  sqrt(mean((Yt - Ft) ^ 2))
}

# Forecasts --------------------------------------------------
n <- length(spx2_rvol2)
## |-- MIDAS rolling forecast ----------
m_pred <- NULL
system.time(
  for(j in 3000:n){
    m_pred[[j+1]] <- as.numeric(as.data.frame(forecast(midas_r(rv ~ mls(rv, 1:20, 1, harstep),
                                                            data = list(rv = spx2_rvol2[(j-3000):j]),
                                                            start = list(rv = c(1, 1, 1))), newdata = list(rv = rep(NA, 1)))))
  }
)

## |-- Naive rolling forecast ----------
n_pred <- NULL
for(i in 3000:n) {
  n_pred[[i + 1]] <- as.numeric(naive(spx2_rvol2[(i - 3000):i], h = 1)[4][1])
}

# merging observed data and forecasts --------------------------------------------------
dataa <- data.frame("real" = spx2_rvol2[3001:n],
                    "mid_pred" = m_pred[3001:n],
                    "naiv_pred" = n_pred[3001:n]
)
# calculating accuracy measures (MSE, MAPE, MASE) --------------------------------------------------
data.frame("mid_mse" = mse(dataa[, 1], dataa[, 2]),
           "mid_mape" = mape(dataa[, 1], dataa[, 2]),
           "mid_mase" = mase(dataa[, 1], dataa[, 2]),
           "naiv_mse" = mse(dataa[, 1], dataa[, 3]),
           "naiv_mape" = mape(dataa[, 1], dataa[, 3]),
           "naiv_mase" = mase(dataa[, 1], dataa[, 3])
)



system.time(
  forc <- average_forecast(list(midas_r(rv ~ mls(rv, 1:20, 1, harstep),
                                        data = list(rv = spx2_rvol2[(j-3000):j]),
                                        start = list(rv = c(1, 1, 1)))), data = list(rv = spx2_rvol2),
                           insample = 1:3000, outsample = 3001:n, type = "rolling",
                           show_progress = TRUE)
)
forc$accuracy$individual
