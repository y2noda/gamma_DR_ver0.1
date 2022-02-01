library(tidyverse)
library(ggplot2)
library(MASS)
library(lattice)
library(scatterplot3d)
library(nleqslv)

# シュミレーション開始
rm(list = ls(all.names = TRUE))

# サンプルサイズ
n = 200

# シミュレーションの繰り返し数
kk_T <- 1000

pb <- txtProgressBar(min = 1, max = kk_T, style = 3)

results.beta_t <- NULL
results.gamma <- NULL

for (i in 1:kk_T) {
  setTxtProgressBar(pb,i)
  
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
  Y_lm1 <- 0.2+1+XX[,1]+XX[,2]+XX[,3]+XX[,4]; Y_lm0 <- 0.2+XX[,1]+XX[,2]+XX[,3]+XX[,4]
  p_Y1 <- exp(Y_lm1)/(1+exp(Y_lm1)); p_Y0 <- exp(Y_lm0)/(1+exp(Y_lm0))
  
  # 傾向スコアの算出
  ps_fit <- glm(TT~XX, family=binomial)$fit
  
  # ZZ
  ZZ <- XX[,1:2]
  ZZ_tilde <- XX[,3:4]
  
  
  
  
  
  #　誤判別を含むアウトカムデータの生成
  ## Misclassified outcomesの割合の設定
  pp01_t <- 0.1; pp10_t <- 0.1
  pp01_c <- 0.1; pp10_c <- 0.1
  
  YY1 <- YY1_t <- YY0 <- YY0_t <- YY <- YY_t <- rep(0,n);
  for(i in 1:n){
    # potential outcomesの乱数の発生
    YY1[i] <- YY1_t[i] <-rbinom(1, size=1, prob=p_Y1[i])
    YY0[i] <- YY0_t[i] <- rbinom(1, size=1, prob=p_Y0[i])
    
  #   Mis_p <- rbinom(1,size=1,prob=pp01_t)
  #   if(Mis_p==1) {
  #     YY1[i] <- ifelse(YY1[i]==1, 0, 1)
  # }
    
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
  
  
  # 推定方程式
  estimate_eq <- function(beta){
    beta_0 <- beta[1]
    beta_t <- beta[2]
    beta_x <- beta[3:6]
    
    or_ml <- beta_0 + TT*beta_t + XX%*%beta_x
    lp <- exp((gamma+1)*(or_ml))
    lq <- exp(YY*(gamma+1)*(or_ml))
    pi <- lp/(1+lp)
    g <- (lq/(1+lp))^(gamma/(1+gamma))
    # w <- ifelse(TT==1, 1/ps_fit, 1/(1-ps_fit))
    
    return(c(t(g*(YY-pi))%*%rep(1,n), t(g*(YY-pi))%*%TT, t(g*(YY-pi))%*%XX))
  }
  
  
  gamma <- 0
  # パラメータ推定 初期値：c(0,0,0,0)
  estimate_pra_fn <- function(gamma){
    gamma <<- gamma
    beta_hat <- nleqslv(c(0.2,1,1,1,1,1), estimate_eq, method = "Newton")$x
    
    return(list(beta_hat=beta_hat, gamma=gamma))
  }
  
  # gamma_fn2 <- function(gamma){
  #   
  #   beta_hat <- estimate_pra_fn(gamma)$beta_hat
  #   
  #   or_ml <- beta_hat[1] + TT*beta_hat[2] + XX%*%beta_hat[3:6]
  #   lp <- exp(or_ml)
  #   pi <- lp/(1+lp)
  #   
  #   C_gamma <- (pi^(1+gamma) + (1-pi)^(1+gamma))^(gamma/(gamma+1))
  #   
  #   f0 <- (pi^2/(1-pi))^gamma - (1-pi)^gamma
  #   f1 <- pi^gamma - ((1-pi)^2/pi)^gamma
  #   
  #   H_gamma <- ifelse(YY==0, f0*pi*log(pi/(1-pi))/(C_gamma*2*gamma), 
  #                     (1)*f1*(1-pi)*log(pi/(1-pi))/(C_gamma*2*gamma))
  #   
  #   return(sum(H_gamma))
  # }
  # 
  # gamma_hat <- 0
  
  # gamma_hat <- optimise(gamma_fn2,c(0,2))$minimum
  beta_hat <- estimate_pra_fn(2)$beta_hat
  
  beta_t <- beta_hat[2]
  
  # results.gamma <- append(results.gamma, gamma_hat)
  results.beta_t <- append(results.beta_t, beta_t)
}

results.beta_t %>% summary()
boxplot(results.beta_t)

results.gamma %>% summary()

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4548  0.6047  0.6601  0.6642  0.7230  0.8835 

## アウトプット
export_data <- data.frame(results.beta_t)
# write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/beta_hat_m00.csv")
# 
export_data <- data.frame(results.gamma)
# write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/gamma_hat_m00.csv")

## input
df <- read_csv2("~/Projects/gamma_DR_ver0.1/results/beta_hat_m03.csv") 
df %>% summary()


## data加工
make_df_fn <- function(df, est_type, mp_num, file_path){
  
  beta <- read_csv2(file_path) 
  
  n <- nrow(beta)
  
  type <- rep(est_type, n)
  mp <- rep(mp_num, n)
  
  new_df <- data.frame(type, mp, beta)
  
  names(new_df) <- c("Type","MP","Beta")
  
  results <- rbind(df, new_df) 
  
  return(results)
}

est_types <- c("naive", "hat")
mp_nums <- c("00","005","01","015","02","03")

df <- data.frame()
for (est_type in est_type) {
  for(mp_num in mp_nums){
    file_path <- paste("~/Projects/gamma_DR_ver0.1/results/", "beta_",est_type ,"_m" ,mp_num ,".csv" ,sep="")
    df <- make_df_fn(df, est_type,mp_num,file_path)
  }
}
# write_csv2(df, file = "~/Projects/gamma_DR_ver0.1/results/all_beta.csv")

df %>% dim


## plot

library(gt)

# path <- "~/gamma_DR_ver0.1/results/all_beta.csv"
path <- "~/Projects/gamma_DR_ver0.1/results/all_beta.csv"
dat <- read_csv2(path)


p2 <- ggplot(dat, aes(x = MP, y = Beta, fill = factor(Type))) +
  geom_boxplot() +
  scale_y_continuous(name = "Beta") +
  scale_x_discrete(labels = abbreviate, name = "Miss probability")

p2

dat <- dat %>% 
  group_by(Type, MP) %>% 
  summarise(mean = mean(Beta),
            sd = sd(Beta)) %>% 
  ungroup()

beta_table <- dat %>%
  gt(rowname_col = "Type", groupname_col = "MP")

beta_table
