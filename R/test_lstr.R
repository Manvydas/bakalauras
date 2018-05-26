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
  
  all_coef <- function(p) {
    weight(p,d)
  }
  
  lstr <- function(p){
    b0 <- p[1]; b1 <- p[2]; b2 <- p[3]; b3 <- p[4]; b4 <- p[5]
    xm <- x%*%(all_coef(p))
    mod <- b0 + b1*xm*(1 + b2 * (1 + exp(-b3 * (xm - b4) / sd(x)))^{-1})
    sum((y-mod)^2)
  }
  
  cc <- optimx(par=startx, lstr, method=method)
  ccmin <- which.min(cc$value)
  cc <- cc[ccmin, ]
  par <- as.numeric(cc[ ,1:length(startx)])
  ress <- y - x%*%all_coef(par)
  
  out <- list("data" = cbind(y,x),
              "coef" = par,
              "mid_coef" = all_coef(par),
              "RSS" = cc$value,
              "residuals" = ress)
  
  out
}

mr <- lstr_m(dat[,1], dat[,-1], weight = nbeta, startx = c(1,1,1,1,1))

summary(mr$data[ ,1])
summary(mr$residuals)
