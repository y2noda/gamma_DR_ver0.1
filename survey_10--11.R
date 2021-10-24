library(nleqslv)
library(tidyverse)
library(ggplot2)
library(bbmle)
library(MASS)
library(Rsolnp)
# fn <- function(x){
#   a <- x[1]
#   b <- x[2]
#   
#  ans <- c() 
#  ans[1] <- a+b-10
#  ans[2] <- a-b-6
#  
#  return(ans)
# }
# nleqslv(c(1,1), fn)

####シュミレーション設定####
rm(list = ls(all.names = TRUE))

# シミュレーションの繰り返し数
kk_T <- 1000
# kk_T <- 1

# サンプルサイズ
n <- 1000

kekka1 <- kekka2 <- matrix(0,ncol=kk_T,nrow=3)
count <- 0


# 共変量の分布の設定、乱数の発生
mu <- c(0,0)

## 相関パラメータ
rho <- 0
Sigma <- matrix(c(1,rho,rho,1),2,2)
XX <- mvrnorm(n,mu=mu,Sigma=Sigma)

# 処置の分布の設定
T_lm <- 0.2+XX[,1]+XX[,2]
p_T <- exp(T_lm)/(1+exp(T_lm))

# Potential outcomesの分布の設定
# アウトカムの回帰パラメータの設定
Y_lm1 <- 0.2+1+XX[,1]+XX[,2]; Y_lm0 <- 0.2+XX[,1]+XX[,2]
p_Y1 <- exp(Y_lm1)/(1+exp(Y_lm1)); p_Y0 <- exp(Y_lm0)/(1+exp(Y_lm0))



# Misclassified outcomesの割合の設定
pp11 <- 0.7; pp10 <- 0.3
# pp11 <- 0.8; pp10 <- 0.2
# pp11 <- 1.0; pp10 <- 0.0


TT <- YY1 <- YY0 <- YY <- Mis_p <- rep(0,n);
for(i in 1:n){
  # 処置、potential outcomesの乱数の発生
  TT[i] <- rbinom(1,size=1,prob=p_T[i])
  YY1[i] <- rbinom(1,size=1,prob=p_Y1[i])
  YY0[i] <- rbinom(1,size=1,prob=p_Y0[i])
  
  # 観測されるoutcome
  YY[i] <- TT[i]*YY1[i]+(1-TT[i])*YY0[i]
  
  # Misclassified outcomesの設定
  if(YY[i]==1){
    Mis_p <- rbinom(1,size=1,prob=(1-pp11))
    if(Mis_p==1) YY[i] <- 0
  }
  
  if(YY[i]==0){
    Mis_p <- rbinom(1,size=1,prob=pp10)
    if(Mis_p==1) YY[i] <- 1 
  }
}

# 真の因果効果 (通常は推定不可能)
tau0 <- mean(YY1)-mean(YY0)

# Error-prone covatiatesの設定
# Errorの分散の設定
# sigma_e <- 1 
# ee <- rnorm(n,sd=sigma_e)

# XX[,1] <- XX[,1]+ee

####work10/10####

pb <- txtProgressBar(min = 1, max = kk_T, style = 3)

results.beta_t <- NULL
gamma_list <- c(0, 0.1, 0.2, 0.3)

for (j in 1:length(gamma_list)) {

beta_t_list <- NULL

for (i in 1:kk_T) {
setTxtProgressBar(pb,i)

#### データ生成 ####
TT <- YY1 <- YY0 <- YY <- Mis_p <- rep(0,n);
for(i in 1:n){
  # 処置、potential outcomesの乱数の発生
  TT[i] <- rbinom(1,size=1,prob=p_T[i])
  YY1[i] <- rbinom(1,size=1,prob=p_Y1[i])
  YY0[i] <- rbinom(1,size=1,prob=p_Y0[i])
  
  # 観測されるoutcome
  YY[i] <- TT[i]*YY1[i]+(1-TT[i])*YY0[i]
  
  # Misclassified outcomesの設定
  if(YY[i]==1){
    Mis_p <- rbinom(1,size=1,prob=(1-pp11))
    if(Mis_p==1) YY[i] <- 0
  }
  
  if(YY[i]==0){
    Mis_p <- rbinom(1,size=1,prob=pp10)
    if(Mis_p==1) YY[i] <- 1 
  }
}

####アウトカム回帰のみのモデル####
fn <- function(beta){
  beta_0 <- beta[1]
  beta_t <- beta[2]
  beta_x <- beta[3:4]
  
  or_ml <- beta_0 + TT*beta_t + XX%*%beta_x
  lp <- exp((gamma+1)*(or_ml))
  lq <- exp(YY*(gamma+1)*(or_ml))
  pi <- lp/(1+lp)
  g <- (lq/(1+lp))^(gamma/(1+gamma))
  return(c(t(g*(YY-pi))%*%rep(1,1000), t(g*(YY-pi))%*%TT, t(g*(YY-pi))%*%XX))
}

# パラメータ推定 c(0,1,1,1)
estimate_pra_fn <- function(gamma){
  return(beta_hat <- nleqslv(c(0,1,1,1), fn)$x)
}

gamma<- gamma_list[j]
beta_hat <- estimate_pra_fn(gamma)

####因果効果の推定####

lm1 <- exp(beta_hat[1] + TT*beta_hat[2] + XX%*%beta_hat[3:4])
lm1 <- lm1[TT==1]
psy1 <- lm1/(lm1+1)
lm2 <- exp(beta_hat[1] + XX%*%beta_hat[3:4])
lm2 <- lm2[TT==0]
psy2 <- lm2/(lm2+1)

tau <- mean(psy1)-mean(psy2)
# tau #0.21

beta_t <- beta_hat[2]
beta_t_list <- c(beta_t_list, beta_t)
}

results.beta_t <- cbind(results.beta_t, beta_t_list)
}
colnames(results.beta_t) <- gamma_list

export_data <- data.frame(results.beta_t)
write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/result.csv")

export_data %>% summary() %>% xtable(label="tb-ref")
####ガンマの推定####

# minus.log.lik <- function(gamma){
#   or_ml <- beta_hat[1] + TT*beta_hat[2] + XX%*%beta_hat[3:4]
#   lp <- exp((gamma+1)*or_ml)
#   # lq <- exp(YY*(gamma+1)*XX%*%beta_hat + phi_hat/ps_ml)
#   pi <- lp/(1+lp)
#   return(-sum(YY*pi + (1-YY)*(1-pi)))
# }
# 
# res <- mle2(minus.log.lik, start = list(gamma=0))
# coef(res)

# gamma_fn <- function(gamma){
#   or_ml <- beta_hat[1] + TT*beta_hat[2] + XX%*%beta_hat[3:4]
#   lp <- exp((gamma+1)*or_ml )
#   lq <- exp(YY*(gamma+1)*or_ml)
#   pi <- lp/(1+lp)
#   return(-sum(pi^(gamma+1)+(1-pi)^(gamma+1))^(1/(gamma+1)))
# }
# gamma_hat <- optimise(gamma_fn,c(0,10))$minimum


# 参考1
# data1.YY <- YY[TT==1]
# data0.YY <- YY[TT==0]
# data1.XX <- XX[TT==1,]
# data0.XX <- XX[TT==0,]
# y1_logit<- glm(data1.YY~data1.XX,family = binomial)$fit
# y0_logit<- glm(data0.YY~data0.XX,family = binomial)$fit

# tau <- mean(y1_logit) - mean(y0_logit)
# tau #0.21
# 
# #　参考2
# TT_XX <- cbind(TT,XX)
# 
# coef <- glm(YY~TT + XX[,1] + XX[,2], family = binomial)$coef
# 
# beta_hat <- coef[2:4]
# lm1 <- exp(TT*beta_hat[1] + XX%*%beta_hat[2:3] + coef[1])
# lm1 <- lm1[TT==1]
# psy1 <- lm1/(lm1+1)
# lm2 <- exp(XX%*%beta_hat[2:3] + coef[1])
# lm2 <- lm2[TT==0]
# psy2 <- lm2/(lm2+1)
# 
# tau <- mean(psy1)-mean(psy2)
# tau #0.21
# 
# # 参考3
# yy_fit <- glm(YY~TT + XX[,1] + XX[,2], family = binomial)$fit
# yy_1 <- yy_fit[TT==1]
# yy_0 <- yy_fit[TT==0]
#     
# tau <- mean(yy_1)- mean(yy_0)
# tau #0.21



####経験尤度法####
#経験尤度法によって，wの分布とガンマハットを求める

# test of solnp
# ObjFunc = function(x) return( - x[1]^(2/5) * x[2]^(3/5))
# 
# ConstFunc = function(x) return( x[1]*4 + x[2] * 6)
# eq.value <- c(100)
# 
# x0 <- c(1,1)
# 
# solution <- solnp(x0, fun = ObjFunc, eqfun = ConstFunc, eqB = eq.value)
# 
# solution$pars


# 実装
# n=1000
# gamma = 0
# ObjFunc = function(x) {
#   omega <- x[5:1004]
#   return(prod(omega)*n)
# }
# 
# ConstFunc = function(x) {
#   beta_0 <- x[1]
#   beta_t <- x[2]
#   beta_x <- x[3:4]
#   omega <- x[5:1004]
#   
#   or_ml <- beta_0 + TT*beta_t + XX%*%beta_x
#   lp <- exp((gamma+1)*(or_ml))
#   lq <- exp(YY*(gamma+1)*(or_ml))
#   pi <- lp/(1+lp)
#   g <- (lq/(1+lp))^(gamma/(1+gamma)) 
#   
#   h <- t(omega * (YY-mean(pi))) %*% TT
#   
#   return(c(t(omega*g*(YY-pi))%*%rep(1,1000), t(omega*g*(YY-pi))%*%TT, t(omega*g*(YY-pi))%*%XX, h))
# }
# 
# eq.value <- rep(0,5)
# 
# x0 <- rep(1, 1004)
# 
# solution <- solnp(x0, fun = ObjFunc, eqfun = ConstFunc, eqB = eq.value)
# 
# pars <- solution$pars
# omega_hat <- pars[5:1004]
# beta_hat <- pars[1:4]
# 
# lm1 <- exp(beta_hat[1] + TT*beta_hat[2] + XX%*%beta_hat[3:4])
# psy1 <- lm1/(lm1+1)
# lm2 <- exp(beta_hat[1] + XX%*%beta_hat[3:4])
# psy2 <- lm2/(lm2+1)
# 
# tau <- mean(psy1)-mean(psy2)
# tau #0.0518
