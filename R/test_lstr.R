library(midasr)
library(dplyr)

# preparing data
data("rvsp500", package = "midasr")
rr <- rvsp500 %>% select(DateID,SPX2.rv) %>% mutate(rv = 100 * sqrt(252*SPX2.rv))
rv1 <- na.omit(rr$rv)
rv2 <- mls(rv1, 1:20, 1)
dat <- na.omit(cbind("y" = rv1, rv2))

# lstr function
lstr_m <- function(y, x, z=NULL, weight, startx, method=c("Nelder-Mead","BFGS")){
  d <- ncol(x)
  
  # weight function
  all_coef <- function(p) {
    weight(p,d)
  }
  
  GG <- function(p){
   b3 <- p[4]; b4 <- p[5]
    xm <- x%*%(all_coef(p))
    (1 + exp(-b3 * (xm - b4) / sd(x)))^{-1}
  }
  # midasr
  lstr <- function(p){
    b0 <- p[1]; b1 <- p[2]; b2 <- p[3]
    xm <- x%*%(all_coef(p))
    mod <- b0 + b1*xm*(1 + b2 * GG(p))
    sum((y-mod)^2)
  }
  
  cc <- optimx(par=startx, lstr, method=method)
  ccmin <- which.min(cc$value)
  cc <- cc[ccmin, ]
  par <- as.numeric(cc[ ,1:length(startx)])
  fitt <- x%*%all_coef(par)
  ress <- y - x%*%all_coef(par)
  
  out1 <- list("data" = cbind(y,x),
              "coef" = par,
              "mid_coef" = all_coef(par),
              "RSS" = cc$value,
              "fitted" = fitt,
              "residuals" = y-fitt,
              "mse" = mean((y-fitt) ^ 2))
  
  # Galvao
  lstr <- function(p){
    b0 <- p[1]; b1 <- p[2]; b2 <- p[3]
    xm <- x%*%(all_coef(p))
    mod <- b0 + b1*xm + (b2-b1)*GG(p)
    sum((y-mod)^2)
  }
  
  cc <- optimx(par=startx, lstr, method=method)
  ccmin <- which.min(cc$value)
  cc <- cc[ccmin, ]
  par <- as.numeric(cc[ ,1:length(startx)])
  fitt <- x%*%all_coef(par)
  ress <- y - x%*%all_coef(par)
  
  out2 <- list("data" = cbind(y,x),
               "coef" = par,
               "mid_coef" = all_coef(par),
               "RSS" = cc$value,
               "fitted" = fitt,
               "residuals" = y-fitt,
               "mse" = mean((y-fitt) ^ 2))
  
  list("midasr" = out1, "galvao" = out2)
}

mr <- lstr_m(dat[,1], dat[,-1], weight = nbeta, startx = c(1,1,1,1,1))
mrs <- midas_r_simple(dat[,1], dat[,-1], weight = nbeta, startx = c(1,1,1))

# data
summary(mr$midasr$data[ ,1])
# Min.  | 1st Qu. | Median | Mean   | 3rd Qu. | Max. 
# 3.386 | 9.024   | 12.813 | 15.390 | 18.095  | 139.729 

# residuals
summary(mr$midasr$residuals)
summary(mr$galvao$residuals)
#           Min.      | 1st Qu.  | Median   | Mean    | 3rd Qu. | Max. 
# midasr:   -8.175    | 3.299    | 5.503    | 7.176   | 8.898   | 109.331 
# galvao:   -34.65499 | -2.16458 | -0.02902 | 0.50855 | 2.43776 | 81.18002 
# r_simple: -35.9014  | -2.4311  | -0.2589  | 0.2326  | 2.1762  | 80.5053
# mse

mr$midasr$mse # 97.32877
mr$galvao$mse # 31.08381
mean((mrs$model[,1]-mrs$model[,-1]%*%mrs$midas_coefficients)^2) # 30.98337

