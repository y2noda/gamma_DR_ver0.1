library(tidyverse)
library(ggplot2)
library(MASS)

# シュミレーション開始
rm(list = ls(all.names = TRUE))

# サンプルサイズ
n = 500

# シミュレーションの繰り返し数
kk_T <- 100

pb <- txtProgressBar(min = 1, max = kk_T, style = 3)

# 推定量の数
nval <- 1

results.beta_hat <- data.frame(matrix(vector(), kk_T, nval))
colnames(results.beta_hat) <- c("b0")
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
    # TT[j] <- rbinom(1,size=1,prob=p_T[j])
    TT[j] <- rbinom(1,size=1,prob=0.5)
    if(TT[j]==0){
      TT[j]=-1
    }
  }
  
  
  # Potential outcomesの分布の設定
  # アウトカムの回帰パラメータの設定

  Y_lm <- 5*(XX[,1]+XX[,2]+XX[,3]+XX[,4]+XX[,5]+XX[,6]+XX[,7]+XX[,8]+XX[,9])*TT/2
  p_Y <- exp(Y_lm)/(1+exp(Y_lm))
  
  # ZZ
  # ZZ <- XX[,1]
  # ZZ_tilde <- XX[,2:9]
  
  #　誤判別を含むアウトカムデータの生成
  ## Misclassified outcomesの割合の設定
  pp10 <- 0.1; pp01 <- 0.1
  
  YY <- rep(0,n)
  for(j in 1:n){
    p_Yc <- pp10*(1-p_Y[j])+(1-pp01)*p_Y[j]
    YY[j] <-rbinom(1, size=1, prob=p_Yc)
  }
  
  
  # # 初期値とガンマの設定
  # # gamma <- 0.00001
  # gamma <- 2
  # b1 <- rep(2,nval)
  # 
  # 
  # # 推定
  # res_df <- data.frame(t(gammat_logistic_nt(YY, ZZ, TT, gamma, b1)$b1))
  # colnames(res_df) <- c("t","b0","b1")
  # 
  # results.beta_hat[i,] <- res_df
  
  
  # optim関数を使った場合
  # gamma <- 2
  gamma <- 0.0001
  b1 <- rep(1,nval)
  
  # 傾向スコアの算出
  # ps_fit <- glm(TT~XX, family=binomial)$fit
  # MSMのためのウェイト
  # ps_w <- ifelse(TT==1, 1/ps_fit, 1/(1-ps_fit))
  

  
  W <- (XX[,1]+XX[,2]+XX[,3]+XX[,4]+XX[,5]+XX[,6]+XX[,7]+XX[,8]+XX[,9])*TT/2
  
  
  loss <- function(beta){
    ei_gamma <- exp((gamma+1)*W*beta)
    val <- sum( (ei_gamma^YY)/(1+ei_gamma))^(gamma/(gamma+1) )
    return(val)
  }
  
  gr <- function(beta){
    n <- dim(XX)[1]
    p <- dim(XX)[2]
    ei_gamma <- exp((gamma+1)*(W*beta))
    wi <- ( (ei_gamma^YY)/(1+ei_gamma) )^(gamma/(gamma+1))
    pi <- ei_gamma/(1+ei_gamma)
    S <- W%*%diag(as.vector(wi/n)) %*% (YY-pi)
    return(S) 
  }
  
  # 準ニュートン法
  res <- optim(b1, loss, method="BFGS", gr=gr, control = list(fnscale = -1))
  res_df <- data.frame(t(res$par))
  colnames(res_df) <- c("b0")
  
  
  results.beta_hat[i,] <- res_df
  
}


results.beta_hat %>% summary()