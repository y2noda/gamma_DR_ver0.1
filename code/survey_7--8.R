library(nleqslv)
library(tidyverse)
library(ggplot2)
library(bbmle)
library(MASS)
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
# kk_T <- 1000
kk_T <- 1

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
pp11 <- 0.9; pp10 <- 0.1
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

####work 8月ぐらい####
####傾向スコア測定のためのロジスティック回帰モデル####
ps_fit <- glm(TT~XX, family=binomial)$fit

y_logit <- glm(YY~XX,family = binomial)$fit

#パラメータ推定
# beta <- rep(1,ncol(XX))
# 
# for (i in 1:n) {
#   fn <- function(beta){
#     lp <- exp(beta%*%XX[i,])
#     pi <- lp/(1+lp)
#     return((YY[i] - pi)%*%XX[i,])
#   }
# }
# 
# i<-1
# fn <- function(beta){
#   lp <- exp(beta%*%XX[i,])
#   pi <- lp/(1+lp)
#   return((YY[i] - pi)%*%XX[i,])
# }
# beta <- nleqslv(c(1,1), fn)$x

####ロジスティック回帰をnleqslvで実行　テスト####
fn <- function(beta){
  lp <- exp(XX%*%beta)
  pi <- lp/(1+lp)
  return(t(YY - pi)%*%XX)
}
beta_hat <- nleqslv(c(1,1), fn)$x
lp <- exp(XX%*%beta_hat)
pi <- lp/(1+lp)

####提案手法####

ps_fit <- glm(TT~XX, family=binomial)$fit

ps_ml <- TT*ps_fit - (1-TT)*(1-ps_fit)

gamma<-0
# gamma <- gamma_hat



# beta_fn <- function(beta){
#   lp <- exp((gamma+1)*XX%*%beta + phi/ps_ml)
#   lq <- exp(YY*(gamma+1)*XX%*%beta + phi/ps_ml)
#   pi <- lp/(1+lp)
#   g <- (lq/(1+lp))^(gamma/(1+gamma))
#   return(t(g*(YY-pi))%*%XX)
# }
# 
# phi_fn <- function(phi){
#   lp <- exp((gamma+1)*XX%*%beta_hat + phi/ps_ml)
#   lq <- exp(YY*(gamma+1)*XX%*%beta_hat + phi/ps_ml)
#   pi <- lp/(1+lp)
#   g <- (lq/(1+lp))^(gamma/(1+gamma))
#   return(t(g*(YY-pi))%*%(ps_ml)^(-1))
# }

# YY1のパラメータ推定
# beta_hat <- nleqslv(c(1,1), beta_fn)$x
# phi_hat <- nleqslv(1, phi_fn)$x

fn <- function(theta){
  beta <- theta[1:2]
  phi <- theta[3]
  
  or_ml <- XX%*%beta + TT
  lp <- exp((gamma+1)*(or_ml + phi/ps_ml))
  lq <- exp(YY*(gamma+1)*(or_ml + phi/ps_ml))
  pi <- lp/(1+lp)
  g <- (lq/(1+lp))^(gamma/(1+gamma))
  return(c(t(g*(YY-pi))%*%XX, t(g*(YY-pi))%*%(ps_ml)^(-1)))
}


# パラメータ推定 c(1,1,0.5)
estimate_pra_fn <- function(gamma){
  return(theta_hat <- nleqslv(c(1,1,0.5), fn)$x)
}
theta_hat <- estimate_pra_fn(gamma)
beta_hat <- theta_hat[1:2]
phi_hat <- theta_hat[3]



# lm1 <- exp(XX%*%beta_hat + TT + phi_hat/ps_fit)
# psy1 <- lm1/(lm1+1)
# lm2 <- exp(XX%*%beta_hat + phi_hat/-(1-ps_fit))
# psy2 <- lm2/(lm2+1)
# 
# tau <- mean(psy1)-mean(psy2)
# tau
# tau0

####gammaの推定####
#ver1 (need bbmle)
minus.log.lik <- function(gamma){
  or_ml <- XX%*%beta_hat + TT
  lp <- exp((gamma+1)*XX%*%beta_hat + phi_hat/ps_ml)
  # lq <- exp(YY*(gamma+1)*XX%*%beta_hat + phi_hat/ps_ml)
  pi <- lp/(1+lp)
  return(-sum(YY*pi + (1-YY)*(1-pi)))
}

res <- mle2(minus.log.lik, start = list(gamma=0))
coef(res)

# gamma ver2 
gamma_fn <- function(gamma){
  or_ml <- XX%*%beta_hat + TT
  lp <- exp((gamma+1)*(or_ml + phi_hat/ps_ml))
  lq <- exp(YY*(gamma+1)*(or_ml + phi_hat/ps_ml))
  pi <- lp/(1+lp)
  return(-sum(pi^(gamma+1)+(1-pi)^(gamma+1))^(1/(gamma+1)))
}
gamma_hat <- optimise(gamma_fn,c(0,10))$minimum

# gamma ver3
# gamma_fn <- function(gamma){
#   lp <- exp((gamma+1)*XX%*%beta_hat + phi_hat/ps_ml)
#   lq <- exp(YY*(gamma+1)*XX%*%beta_hat + phi_hat/ps_ml)
#   pi <- lp/(1+lp)
#   g <- (lq/(1+lp))^(gamma/(1+gamma))
#   return(-sum(g))
# }
# gamma_hat <- optimise(gamma_fn,c(0,10))$minimum

####因果効果の推定####
theta_hat <- estimate_pra_fn(gamma)

beta_hat <- theta_hat[1:2]
phi_hat <- theta_hat[3]

lm1 <- exp(XX%*%beta_hat + TT + phi_hat/ps_fit)
lm1 <- lm1[TT==1]
psy1 <- lm1/(lm1+1)
lm2 <- exp(XX%*%beta_hat + phi_hat/-(1-ps_fit))
lm2 <- lm2[TT==0]
psy2 <- lm2/(lm2+1)

tau <- mean(psy1)-mean(psy2)
tau
tau0

####通常のDR推定値####
ps_fit <- glm(TT~XX, family=binomial)$fit

data1.YY <- YY[TT==1]
data0.YY <- YY[TT==0]
data1.XX <- XX[TT==1,]
data0.XX <- XX[TT==0,]
y1_logit<- glm(data1.YY~data1.XX,family = binomial)$fit
y0_logit<- glm(data0.YY~data0.XX,family = binomial)$fit

dre1 <- (1/1000)*sum(YY+((TT-ps_fit)/ps_fit)*(YY-y1_logit))
dre0 <- (1/1000)*sum(((1-TT)*YY)/(1-ps_fit)+(1-(1-TT)/(1-ps_fit))*y0_logit)

tau_dr <- dre1 - dre0
# 0.18 tau=0.171

































