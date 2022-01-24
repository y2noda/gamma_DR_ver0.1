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
pp01 <- 0.35; pp10 <- 0.35

YY1 <- YY1_t <- YY0 <- YY0_t <- YY <- YY_t <- rep(0,n);
for(i in 1:n){
  # potential outcomesの乱数の発生
  YY1[i] <- YY1_t[i]<- rbinom(1, size=1, prob=p_Y1[i])
  Mis_p1 <- rbinom(1, size=1, prob=pp10)
  
  YY0[i] <- YY0_t[i] <- rbinom(1, size=1, prob=p_Y0[i])
  Mis_p0 <- rbinom(1, size=1, prob=pp01)
  
  if(Mis_p1==1){
    YY1_t[i] <- ifelse(YY1[i]==1, 0, 1)
  }
  
  if(Mis_p0==1){
    YY0_t[i] <- ifelse(YY0[i]==1, 0, 1)
  }

  # 観測されるoutcome
  YY[i] <- TT[i]*YY1[i]+(1-TT[i])*YY0[i]
  YY_t[i] <- TT[i]*YY1_t[i]+(1-TT[i])*YY0_t[i]
}

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

## 結果

gamma <- 0.1
beta_hat <- estimate_pra_fn(gamma)$beta_hat
print(beta_hat)
plot_fn(beta_hat, gamma, type = "hist")



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


## Sugasawa2021
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
  
  H_gamma <- ifelse(YY==0, f0*pi*log(pi/(1-pi))/(C_gamma*2*gamma), f1*(1-pi)*log(pi/(1-pi))/(C_gamma*2*gamma))

  return(sum(H_gamma))
}

# plot(gamma_fn2, 0, 10)
curve(Vectorize(gamma_fn2)(x), 0, 2)


gamma_hat <- optimise(gamma_fn2,c(0,2))$minimum


# 新たなガンマによる推定結果

beta_hat_gamma_hat <- estimate_pra_fn(gamma_hat)


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


# シュミレーション開始
pb <- txtProgressBar(min = 1, max = kk_T, style = 3)

results.beta_t <- NULL
beta_t_list <- NULL

for (i in 1:kk_T) {
  setTxtProgressBar(pb,i)
  
  #　誤判別を含むアウトカムデータの生成
  ## Misclassified outcomesの割合の設定
  pp01 <- 0.3; pp10 <- 0.3
  
  YY1 <- YY1_t <- YY0 <- YY0_t <- YY <- YY_t <- rep(0,n);
  for(i in 1:n){
    # potential outcomesの乱数の発生
    YY1[i] <- YY1_t[i]<- rbinom(1, size=1, prob=p_Y1[i])
    Mis_p1 <- rbinom(1, size=1, prob=pp10)
    
    YY0[i] <- YY0_t[i] <- rbinom(1, size=1, prob=p_Y0[i])
    Mis_p0 <- rbinom(1, size=1, prob=pp01)
    
    if(Mis_p1==1){
      YY1_t[i] <- ifelse(YY1[i]==1, 0, 1)
    }
    
    if(Mis_p0==1){
      YY0_t[i] <- ifelse(YY0[i]==1, 0, 1)
    }
    
    # 観測されるoutcome
    YY[i] <- TT[i]*YY1[i]+(1-TT[i])*YY0[i]
    YY_t[i] <- TT[i]*YY1_t[i]+(1-TT[i])*YY0_t[i]
  }  
  
  gamma <- 0.1
  
  beta_hat <- estimate_pra_fn(gamma)$beta_hat
  # print(beta_hat)
  # plot_fn(beta_hat, gamma, type = "density")
  
  beta_t <- beta_hat[2]
  
  results.beta_t <- append(results.beta_t, beta_t)
}

results.beta_t %>% summary()
boxplot(results.beta_t)

## アウトプット
# export_data <- data.frame(results.beta_t)
# write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/gamma0-5.csv")
# write_csv2(export_data %>% summary() %>% data.frame(), file = "~/Projects/gamma_DR_ver0.1/results/gamma0-5summary.csv")

gamma <- 0
beta_hat <- estimate_pra_fn(gamma)$beta_hat

lm1 <- exp(beta_hat[1] + TT*beta_hat[2] + ZZ%*%beta_hat[3:4])
psy1 <- lm1/(lm1+1)
lm2 <- exp(beta_hat[1] + ZZ%*%beta_hat[3:4])
psy2 <- lm2/(lm2+1)


tau <- mean(psy1)-mean(psy2)

#  glmでgamma=0のときと一致するかの確認
##0.3678035 0.5080564 0.2977938 0.3353576

w <- ifelse(TT==1, 1/ps_fit, 1/(1-ps_fit))
result <- glm(YY~TT+XX[,1:2], family = binomial, weights = w)
