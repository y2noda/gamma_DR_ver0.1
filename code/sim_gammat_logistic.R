library(tidyverse)
library(ggplot2)
library(MASS)

# シュミレーション開始
rm(list = ls(all.names = TRUE))
source("~/Projects/gamma_DR_ver0.1/code/gammat_logistic_nt.R")

# サンプルサイズ
n = 500

# シミュレーションの繰り返し数
kk_T <- 500

pb <- txtProgressBar(min = 1, max = kk_T, style = 3)

# 推定量の数
nval <- 5

results.beta_hat <- data.frame(matrix(vector(), kk_T, nval))
colnames(results.beta_hat) <- c("t","b0","b1","b2","b3")
# results.gamma <- NULL

for (i in 1:kk_T) {
  setTxtProgressBar(pb,i)
  
  # 共変量の分布の設定、乱数の発生
  mu<- c(0,0,0,0,0,0,0,0,0)
  
  ## 相関パラメータ
  rho <- 0
  Sigma <- matrix(c(1,rho,rho,rho,rho,rho,rho,rho,rho,
                    rho,1,rho,rho,rho,rho,rho,rho,rho,
                    rho,rho,1,rho,rho,rho,rho,rho,rho,
                    rho,rho,rho,1,rho,rho,rho,rho,rho,
                    rho,rho,rho,rho,1,rho,rho,rho,rho,
                    rho,rho,rho,rho,rho,1,rho,rho,rho,
                    rho,rho,rho,rho,rho,rho,1,rho,rho,
                    rho,rho,rho,rho,rho,rho,rho,1,rho,
                    rho,rho,rho,rho,rho,rho,rho,rho,1),9,9)
  
  XX <- mvrnorm(n, mu, Sigma)
  
  # 処置の分布の設定,生成
  T_lm <- 0.2+XX[,1]+XX[,2]+XX[,3]+XX[,4]+XX[,5]+XX[,6]+XX[,7]+XX[,8]+XX[,9]
  p_T <- exp(T_lm)/(1+exp(T_lm))

  TT <- rep(0,n)
  for (j in 1:n) {
    TT[j] <- rbinom(1,size=1,prob=p_T[j])
  }

  # Potential outcomesの分布の設定
  # アウトカムの回帰パラメータの設定
  Y_lm <- 0.2 + TT + XX[,1] - XX[,2] + XX[,3]
  p_Y <- exp(Y_lm)/(1+exp(Y_lm))
  
  # ZZ
  ZZ <- XX[,1:3]
  ZZ_tilde <- XX[,4:9]
  
  #　誤判別を含むアウトカムデータの生成
  ## Misclassified outcomesの割合の設定
  pp10 <- 0.1; pp01 <- 0.05
  
  YY <- rep(0,n)
  for(j in 1:n){
    p_Yc <- pp10*(1-p_Y[j])+(1-pp01)*p_Y[j]
    YY[j] <-rbinom(1, size=1, prob=p_Yc)
  }
  
  # 切片追加
  ZZ <- rep(1,n) %>% cbind(ZZ)
  
  # 初期値とガンマの設定
  # gamma <- 0.00001
  gamma <- 2
  b1 <- rep(1,nval)
  
  # 推定
  res_df <- data.frame(t(gammat_logistic_nt(YY, ZZ, TT, gamma, b1)$b1))
  colnames(res_df) <- c("t","b0","b1","b2","b3")
  
  results.beta_hat[i,] <- res_df
  
  
}


results.beta_hat %>% summary()



## アウトプット
export_data <- results.beta_hat

# write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/beta_g0_m02.csv")
# 
# write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/beta_g4_m01.csv")

# export_data <- data.frame(results.gamma)
# write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/gamma_hat_m00.csv")


df <- read_csv2("~/Projects/gamma_DR_ver0.1/results/beta_g0_m01.csv") 
df %>% summary()
df %>% summary %>% xtable::xtable()


