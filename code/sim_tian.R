library(tidyverse)
library(ggplot2)
library(MASS)
library(glmnet)

# シュミレーション開始
rm(list = ls(all.names = TRUE))　

# set.seed(12345)
# サンプルサイズ
n = 100

# シミュレーションの繰り返し数
kk_T <- 100
# kk_T <- 500

pb <- txtProgressBar(min = 1, max = kk_T, style = 3)

# 共変量行列の次元
## 100以上でglmはエラーになる
nval <- 50

# results.beta_hat <- data.frame(matrix(vector(), kk_T, nval))
results.corr <- data.frame(matrix(NaN,nrow = kk_T, ncol=2))
colnames(results.corr) <- c("Full Regression", "new")
# results.gamma <- NULL

for (i in 1:kk_T) {
  setTxtProgressBar(pb,i)
  
  # 共変量の分布の設定、乱数の発生
  mu<- rep(0,nval)
  
  ## 相関パラメータ
  rho <- 0
  
  Sigma <- diag(nval)
  
  XX <- mvrnorm(n, mu, Sigma)
  
  # 処置の分布の設定,生成
  # T_lm <- 0.2+XX[,1]+XX[,2]+XX[,3]+XX[,4]+XX[,5]+XX[,6]+XX[,7]+XX[,8]+XX[,9]+XX[,10]
  # p_T <- exp(T_lm)/(1+exp(T_lm))
  
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
  
  beta_0 <- replace(rep(0,nval+1), c(1,2,3,4,5), c(0.4, 0.8, -0.8, 0.8, -0.8))
  
  alpha_0 <- replace(rep(0,nval+1), 
                     c(1,4,5,6,7,8,9,10,11),
                     c(sqrt(6)^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1)
  )
  
  # 連続の場合
  # Y_lm <- (rep(alpha_0[1], n) + (XX %*% alpha_0[-1]))^2 + (rep(beta_0[1], n) + (XX %*% beta_0[-1]) + (XX[,1]*XX[,2] *0.8))*TT + sqrt(2)*rnorm(n)
  # # p_Y <- exp(Y_lm)/(1+exp(Y_lm))
  # 
  # score_true <- 1.6 * (0.5 + XX[,1] - XX[,2] + XX[,3] - XX[,4] + XX[,1]*XX[,2])
  
  # ２値の場合
  Y_lm <- (rep(alpha_0[1], n) + (XX %*% alpha_0[-1]))^2 + (rep(beta_0[1], n) + (XX %*% beta_0[-1]) + (XX[,1]*XX[,2] *0.8))*TT + sqrt(2)*rnorm(n)
  # Y <- ifelse(Y_lm >=0, 1,0)
  
  lm <- 1.6 * (0.5 + XX[,1] - XX[,2] + XX[,3] - XX[,4] + XX[,1]*XX[,2])
  score_true <- (exp(lm/2)-1) / (exp(lm/2)+1)
  
  
  
  #　誤判別を含むアウトカムデータの生成
  ## Misclassified outcomesの割合の設定
  p_Y <- exp(Y_lm)/(1+exp(Y_lm))
  pp10 <- 0; pp01 <- 0

  Y <- rep(0,n)
  for(j in 1:n){
    p_Yc <- pp10*(1-p_Y[j])+(1-pp01)*p_Y[j]
    Y[j] <-rbinom(1, size=1, prob=p_Yc)
  }
  
  
  # 傾向スコアの算出
  # ps_fit <- glm(TT~XX, family=binomial)$fit
  # MSMのためのウェイト
  # ps_w <- ifelse(TT==1, 1/ps_fit, 1/(1-ps_fit))
  
  
  ## 連続の場合
  # modified covariance method
  # W <- cbind(rep(1, n), XX)
  # W_star <- W * TT/2
  # 
  # W_star.scaled <- scale(W_star)
  # # beta_hat <- glm(Y_lm ~ W_star +0 , family = gaussian)$coef
  # lasso.model.cv <- cv.glmnet(x = W_star, y = Y_lm, family = "gaussian", alpha = 1, nfolds = 10)
  # 
  # lasso.model <- glmnet(x = W_star, y = Y_lm, family = "gaussian", alpha = 1)
  # 
  # mc_pred <- predict(lasso.model.cv,
  #                      newx = W,
  #                      alpha=1,
  #                      s = lasso.model.cv$lambda.min)
  # 
  # corr.mc <- cor(score_true, mc_pred)
  # 
  # results.corr[i,2] <- corr.mc
  # # 
  # 
  # # Full regression
  # XX_1 <- cbind(rep(1, n), XX)
  # 
  # XX_full <- XX_1
  # for (j in seq(ncol(XX_1))) {
  #   XX_full <- cbind(XX_full, XX_1[,j]*TT)
  # }
  # 
  # XX_full.scaled <- scale(XX_full)
  # 
  # lasso.model.cv <- cv.glmnet(x = XX_full, y = Y_lm, family = "gaussian", alpha = 1, nfolds = 10) 
  # 
  # lasso.model <- glmnet(x = XX_full, y = Y_lm, family = "gaussian", alpha = 1)
  # 
  # XX_full[,1:51] <- 0
  # XX_full <- XX_full * TT
  # 
  # full_pred <- predict(lasso.model.cv,
  #                 newx = XX_full,
  #                 alpha=1,
  #                 s = lasso.model.cv$lambda.min)
  # 
  # # beta <-coef(lasso.model, s = lasso.model.cv$lambda.min)
  # 
  # corr.full <- cor(score_true, full_pred)
  # 
  # results.corr[i,1] <- corr.full
  
  
  
  ## ２値の場合
  # modified covariance method
  W <- cbind(rep(1, n), XX)
  W_star <- W * TT/2
  
  W_star.scaled <- scale(W_star)
  # beta_hat <- glm(Y_lm ~ W_star +0 , family = gaussian)$coef
  lasso.model.cv <- cv.glmnet(x = W_star, y = Y, family = "binomial", alpha = 1, nfolds = 10) %>% plot
  
  lasso.model <- glmnet(x = W_star, y = Y, family = "binomial", alpha = 1)
  
  mc_pred <- predict(lasso.model.cv,
                     newx = W,
                     alpha=1,
                     s = lasso.model.cv$lambda.min)
  
  score_mc <- (exp(mc_pred/2)-1) / (exp(mc_pred/2)+1)
  
  corr.mc <- cor(score_true, score_mc)
  
  results.corr[i,2] <- corr.mc
  # 
  
  # Full regression
  XX_1 <- cbind(rep(1, n), XX)
  
  XX_full <- XX_1
  for (j in seq(ncol(XX_1))) {
    XX_full <- cbind(XX_full, XX_1[,j]*TT)
  }
  
  XX_full.scaled <- scale(XX_full)
  
  lasso.model.cv <- cv.glmnet(x = XX_full, y = Y, family = "binomial", alpha = 1, nfolds = 10) 
  
  lasso.model <- glmnet(x = XX_full, y = Y, family = "binomial", alpha = 1)
  
  XX_full[,1:51] <- 0
  XX_full <- XX_full * TT
  
  full_pred <- predict(lasso.model.cv,
                       newx = XX_full,
                       alpha=1,
                       s = lasso.model.cv$lambda.min)
  
  score_full <- (exp(full_pred/2)-1) / (exp(full_pred/2)+1)
  # beta <-coef(lasso.model, s = lasso.model.cv$lambda.min)
  
  corr.full <- cor(score_true, score_full)
  
  results.corr[i,1] <- corr.full

  
}


results.corr %>% boxplot()


## アウトプット
export_data <- results.corr

# write_csv2(export_data, file = "~/gamma_DR_ver0.1/results/tian/misbin_nval50.csv")

# 
# df <- read_csv2("~/Projects/gamma_DR_ver0.1/results/0423/g0_m01.csv")
# df <- read_csv2("~/gamma_DR_ver0.1/results/tian/misbin_nval50.csv")
# df %>% summary()
# df %>% boxplot()
