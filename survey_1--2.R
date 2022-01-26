library(tidyverse)
library(ggplot2)
library(MASS)
library(lattice)
library(scatterplot3d)
library(nleqslv)

rm(list = ls(all.names = TRUE))

# サンプルサイズ
n = 200

# シミュレーションの繰り返し数
kk_T <- 1000

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

# Potential outcomesの分布の設定
# アウトカムの回帰パラメータの設定
Y_lm1 <- 0.2+1+XX[,1]+XX[,2]; Y_lm0 <- 0.2+XX[,1]+XX[,2]
p_Y1 <- exp(Y_lm1)/(1+exp(Y_lm1)); p_Y0 <- exp(Y_lm0)/(1+exp(Y_lm0))

#　誤判別を含むアウトカムデータの生成
## Misclassified outcomesの割合の設定
pp01_t <- 0.1; pp10_t <- 0.1
pp01_c <- 0.2; pp10_c <- 0.2

YY1 <- YY1_t <- YY0 <- YY0_t <- YY <- YY_t <- rep(0,n);
for(i in 1:n){
  # potential outcomesの乱数の発生
  YY1[i] <- YY1_t[i] <-rbinom(1, size=1, prob=p_Y1[i])
  YY0[i] <- YY0_t[i] <- rbinom(1, size=1, prob=p_Y0[i])
  
  # 観測されるoutcome
  YY[i] <- TT[i]*YY1[i]+(1-TT[i])*YY0[i] #汚染あり
  YY_t[i] <- TT[i]*YY1_t[i]+(1-TT[i])*YY0_t[i] #なし
  
  # Misclassified outcomesの設定
  if(TT[i]==1){
    if(YY_t[i]==1){
      Mis_p <- rbinom(1,size=1,prob=pp01_t)
      if(Mis_p==1) YY[i] <- 0
    }
    
    if(YY_t[i]==0){
      Mis_p <- rbinom(1,size=1,prob=pp10_t)
      if(Mis_p==1) YY[i] <- 1
    }
  }
  if(TT[i]==0){
    if(YY_t[i]==1){
      Mis_p <- rbinom(1,size=1,prob=pp01_c)
      if(Mis_p==1) YY[i] <- 0
    }
    
    if(YY_t[i]==0){
      Mis_p <- rbinom(1,size=1,prob=pp10_c)
      if(Mis_p==1) YY[i] <- 1
    }
  }
  
}

sum(YY == YY_t)

# 真の因果効果 (通常は推定不可能)
tau0 <- mean(YY1_t)-mean(YY0_t)
tau_naive <- mean(YY[TT==1])-mean(YY[TT==0])

# 傾向スコアの算出
ps_fit <- glm(TT~XX, family=binomial)$fit

# ZZ
ZZ <- XX[,1:2]
ZZ_tilde <- XX[,3:4]

####アウトカム回帰のみのモデル####

# 推定方程式
estimate_eq <- function(beta){
  beta_0 <- beta[1]
  beta_t <- beta[2]
  beta_z <- beta[3:4]
  
  or_ml <- beta_0 + TT*beta_t + ZZ%*%beta_z
  lp <- exp((gamma+1)*(or_ml))
  lq <- exp(YY*(gamma+1)*(or_ml))
  pi <- lp/(1+lp)
  g <- (lq/(1+lp))^(gamma/(1+gamma))
  w <- ifelse(TT==1, 1/ps_fit, 1/(1-ps_fit))
  
  return(c(t(w*g*(YY-pi))%*%rep(1,n), t(w*g*(YY-pi))%*%TT, t(w*g*(YY-pi))%*%ZZ))
}

# パラメータ推定 初期値：c(0,0,0,0)
estimate_pra_fn <- function(gamma){
  gamma <<- gamma
  beta_hat <- nleqslv(c(0,0,0,0), estimate_eq, method = "Newton")$x
  
  return(list(beta_hat=beta_hat, gamma=gamma))
}


## 重みのプロット
plot_fn <- function(beta,gamma,type="density"){
  beta_0 <- beta[1]
  beta_t <- beta[2]
  beta_z <- beta[3:4]
  
  or_ml <- beta_0 + TT*beta_t + ZZ%*%beta_z
  lp <- exp((gamma+1)*(or_ml))
  lq <- exp(YY*(gamma+1)*(or_ml))
  pi <- lp/(1+lp)
  g <- (lq/(1+lp))^(gamma/(1+gamma))
  
  if(type == "density"){
    plot(density(g))
  }else if(type == "hist"){
    hist(g)
  }else{
    print("error")
  }
}

## だいたい全体の何割のデータを使っているか？
propotion_fn <- function(beta,gamma){
  beta_0 <- beta[1]
  beta_t <- beta[2]
  beta_z <- beta[3:4]
  
  or_ml <- beta_0 + TT*beta_t + ZZ%*%beta_z
  lp <- exp((gamma+1)*(or_ml))
  lq <- exp(YY*(gamma+1)*(or_ml))
  pi <- lp/(1+lp)
  g <- (lq/(1+lp))^(gamma/(1+gamma))
  
  return(mean(g))
}

## 結果

gamma <- 0.3
beta_hat <- estimate_pra_fn(gamma)$beta_hat
print(beta_hat)
plot_fn(beta_hat, gamma, type = "hist")

propotion_fn(beta_hat, gamma)

# gamma <- 1
# beta_hat <- estimate_pra_fn(gamma)
# 
# beta_t <- beta_hat[2]

# 
# ガンマの選択について
# gamma_0 <- 0.1
# gamma <- 0.1
# 
# gamma_fn1 <- function(gamma){
#   beta_hat <- estimate_pra_fn(gamma)$beta_hat
# 
#   or_ml <- beta_hat[1] + TT*beta_hat[2] + ZZ%*%beta_hat[3:4]
#   lp <- exp(or_ml)
#   pi <- lp/(1+lp)
#   return(-sum(pi^(gamma_0+1)+(1-pi)^(gamma_0+1))^(1/(gamma_0+1)))
# }
# 
# gamma_hat <- optimise(gamma_fn1,c(0,1))$minimum
# 
# ## plot
# curve(Vectorize(gamma_fn1)(x), 0, 10)
# # plot(gamma_fn1, 0, 10)


# Sugasawa2021
DD <- function(expr, name, order = 1) {
  if(order < 1) stop("'order' must be >= 1")
  if(order == 1) D(expr, name)
  else DD(D(expr, name), name, order - 1)
}

f <- expression(pi^YY * (1-pi)^(1-YY))
DD(f, "YY", 2)   

gamma_fn2 <- function(gamma){
  
  beta_hat <- estimate_pra_fn(gamma)$beta_hat

  or_ml <- beta_hat[1] + TT*beta_hat[2] + ZZ%*%beta_hat[3:4]
  lp <- exp(or_ml)
  pi <- lp/(1+lp)
  
  C_gamma <- (pi^(1+gamma) + (1-pi)^(1+gamma))^(gamma/(gamma+1))
  
  f0 <- (pi^2/(1-pi))^gamma - (1-pi)^gamma
  f1 <- pi^gamma - ((1-pi)^2/pi)^gamma
  
  H_gamma <- ifelse(YY==0, f0*pi*log(pi/(1-pi))/(C_gamma*2*gamma), 
                    (1)*f1*(1-pi)*log(pi/(1-pi))/(C_gamma*2*gamma))

  return(sum(H_gamma))
}

## plot(gamma_fn2, 0, 10)
curve(Vectorize(gamma_fn2)(x), 0, 0.7)



gamma_hat <- optimise(gamma_fn2,c(0,2))$minimum


# 新たなガンマによる推定結果

beta_hat_gamma_hat <- estimate_pra_fn(gamma_hat)

beta_hat_gamma_hat

# まとめた推定

## 最適化のための関数にパラメータを入れられない

gamma_0 <- 0.1
estimate_fn <- function(gamma_0){
  
  gamma_fn1 <- function(gamma){
    beta_hat <- estimate_pra_fn(gamma)$beta_hat
    
    or_ml <- beta_hat[1] + TT*beta_hat[2] + ZZ%*%beta_hat[3:4]
    lp <- exp(or_ml)
    pi <- lp/(1+lp)
    return(-sum(pi^(gamma_0+1)+(1-pi)^(gamma_0+1))^(1/(gamma_0+1)))
  }
  
  gamma_hat <- optimise(gamma_fn1,c(0,2.5))$minimum
  beta_hat_gamma <- estimate_pra_fn(gamma_hat)$beta_hat
  
  return(list(beta_hat_gamma=beta_hat_gamma, gamma_hat=gamma_hat))
}

