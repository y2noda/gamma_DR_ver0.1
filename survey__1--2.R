library(tidyverse)
library(ggplot2)
library(MASS)
library(lattice)
library(scatterplot3d)
library(nleqslv)

rm(list = ls(all.names = TRUE))

# サンプルサイズ
n = 1000

# 共変量の分布の設定、乱数の発生
mu <- c(0,0,0,0)

## 相関パラメータ
rho <- 0
Sigma <- matrix(c(1,rho,rho,rho,
                  rho,1,rho,rho,
                  rho,rho,1,rho,
                  rho,rho,rho,1),4,4)

XX <- mvrnorm(n, mu, Sigma)


# 処置の分布の設定,生成
T_lm <- 0.2+XX[,1]+XX[,2]+XX[,3]+XX[,4]
p_T <- exp(T_lm)/(1+exp(T_lm))

TT <- rep(0,n)
for (i in 1:n) {
  TT[i] <- rbinom(1,size=1,prob=p_T[i])
}

betax <- 0.2+TT+XX[,1]+XX[,2]

# 誤判別確率の生成
## アウトカムに誤判別が含まれる割合
mp <- 0.1
mp_size <- n * mp

## 誤判別のアウトカムのインデックス
mp_id <- sort(betax, decreasing = F, index=T)$ix %>% head(mp_size/2) 
mp_id <- mp_id %>% append(sort(betax, decreasing = T, index=T)$ix %>% head(mp_size/2))


# Potential outcomesの分布の設定
# アウトカムの回帰パラメータの設定
Y_lm1 <- 0.2+1+XX[,1]+XX[,2]; Y_lm0 <- 0.2+XX[,1]+XX[,2]
p_Y1 <- exp(Y_lm1)/(1+exp(Y_lm1)); p_Y0 <- exp(Y_lm0)/(1+exp(Y_lm0))


#　誤判別を含むアウトカムデータの生成
## Misclassified outcomesの割合の設定

YY1 <- YY0 <- YY <- rep(0,n);
for(i in 1:n){
  # potential outcomesの乱数の発生
  YY1[i] <- rbinom(1,size=1,prob=p_Y1[i])
  YY0[i] <- rbinom(1,size=1,prob=p_Y0[i])
 
  # 観測されるoutcome
  YY[i] <- TT[i]*YY1[i]+(1-TT[i])*YY0[i]
  
  # Misclassified outcomesの設定
  if(i %in% mp_id){
    if(YY[i]==1) {
      YY[i] <-  0
    }else if(YY[i]==0){
      YY[i] <- 1
    }
  }
}

# sum((TT*YY1 + (1-TT)*YY0) == YY) 

# 真の因果効果 (通常は推定不可能)
tau0 <- mean(YY1)-mean(YY0)
tau_naive <- mean(YY[TT==1])-mean(YY[TT==0])

# 傾向スコアの算出
ps_fit <- glm(TT~XX, family=binomial)$fit

# ZZ
ZZ <- XX[,1:2]
ZZ_tilde <- XX[,3:4]

####アウトカム回帰のみのモデル####
fn <- function(beta){
  beta_0 <- beta[1]
  beta_t <- beta[2]
  beta_z <- beta[3:4]
  
  or_ml <- beta_0 + TT*beta_t + ZZ%*%beta_z
  lp <- exp((gamma+1)*(or_ml))
  lq <- exp(YY*(gamma+1)*(or_ml))
  pi <- lp/(1+lp)
  g <- (lq/(1+lp))^(gamma/(1+gamma))
  w <- ifelse(TT==1, 1/ps_fit, 1/(1-ps_fit))
  
  return(c(t(w*g*(YY-pi))%*%rep(1,1000), t(w*g*(YY-pi))%*%TT, t(w*g*(YY-pi))%*%ZZ))
}

# パラメータ推定 c(0,1,1,1)
estimate_pra_fn <- function(gamma){
  return(beta_hat <- nleqslv(c(0,0,0,0), fn, method = "Newton")$x)
}

gamma <- 0.5
beta_hat <- estimate_pra_fn(gamma)

lm1 <- exp(beta_hat[1] + TT*beta_hat[2] + ZZ%*%beta_hat[3:4])
lm1 <- lm1[TT==1]
psy1 <- lm1/(lm1+1)
lm2 <- exp(beta_hat[1] + ZZ%*%beta_hat[3:4])
lm2 <- lm2[TT==0]
psy2 <- lm2/(lm2+1)

tau <- mean(psy1)-mean(psy2)
beta_t <- beta_hat[2]

#  glm(YY~TT+XX[,1:2], family = binomial)
#0.3463       0.6481       0.9504       0.9646 

w <- ifelse(TT==1, 1/ps_fit, 1/(1-ps_fit))
wY <- w*YY
glm(wY~TT+XX[,1:2], family = binomial)

