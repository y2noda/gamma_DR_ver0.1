library(tidyverse)
library(ggplot2)
library(MASS)
library(lattice)
library(scatterplot3d)
library(nleqslv)

# シュミレーション開始
rm(list = ls(all.names = TRUE))

# サンプルサイズ
n = 500

# シミュレーションの繰り返し数
kk_T <- 500

pb <- txtProgressBar(min = 1, max = kk_T, style = 3)

results.beta_hat <- data.frame(matrix(vector(), kk_T, 9))
colnames(results.beta_hat) <- c("b1","b2","b3","b4","b5","b6","b7","b8","b9")
# results.gamma <- NULL

# library(mlbench)
# df <- PimaIndiansDiabetes2

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
  # T_lm <- 0.2+XX[,1]+XX[,2]+XX[,3]+XX[,4]
  # p_T <- exp(T_lm)/(1+exp(T_lm))
  # 
  # TT <- rep(0,n)
  # for (i in 1:n) {
  #   TT[i] <- rbinom(1,size=1,prob=p_T[i])
  # }
  
  # Potential outcomesの分布の設定
  # アウトカムの回帰パラメータの設定
  # Y_lm1 <- 0.2+1+XX[,1]+XX[,2]+XX[,3]+XX[,4]; Y_lm0 <- 0.2+XX[,1]+XX[,2]+XX[,3]+XX[,4]
  # p_Y1 <- exp(Y_lm1)/(1+exp(Y_lm1)); p_Y0 <- exp(Y_lm0)/(1+exp(Y_lm0))
  
  Y_lm <- XX[,1] - XX[,2] + XX[,3]
  p_Y <- exp(Y_lm)/(1+exp(Y_lm))
  
  # 傾向スコアの算出
  # ps_fit <- glm(TT~XX, family=binomial)$fit
  # 
  # # ZZ
  # ZZ <- XX[,1:2]
  # ZZ_tilde <- XX[,3:4]
  # 
  
  
  
  
  #　誤判別を含むアウトカムデータの生成
  ## Misclassified outcomesの割合の設定
  # pp01_t <- 0.1; pp10_t <- 0.1
  # pp01_c <- 0.1; pp10_c <- 0.1
  
  pp10 <- 0.05; pp01 <- 0.1
  
  YY <- rep(0,n)
  for(j in 1:n){
    # # potential outcomesの乱数の発生
    # YY1[i] <- YY1_t[i] <-rbinom(1, size=1, prob=p_Y1[i])
    # YY0[i] <- YY0_t[i] <- rbinom(1, size=1, prob=p_Y0[i])
    
  #   Mis_p <- rbinom(1,size=1,prob=pp01_t)
  #   if(Mis_p==1) {
  #     YY1[i] <- ifelse(YY1[i]==1, 0, 1)
  # }
    
    # # 観測されるoutcome
    # YY[i] <- TT[i]*YY1[i]+(1-TT[i])*YY0[i] #汚染あり
    # YY_t[i] <- TT[i]*YY1_t[i]+(1-TT[i])*YY0_t[i] #なし
    # 
    # # Misclassified outcomesの設定
    # if(TT[i]==1){
    #   if(YY_t[i]==1){
    #     Mis_p <- rbinom(1,size=1,prob=pp01_t)
    #     if(Mis_p==1) YY[i] <- 0
    #   }
    #   
    #   if(YY_t[i]==0){
    #     Mis_p <- rbinom(1,size=1,prob=pp10_t)
    #     if(Mis_p==1) YY[i] <- 1
    #   }
    # }
    # if(TT[i]==0){
    #   if(YY_t[i]==1){
    #     Mis_p <- rbinom(1,size=1,prob=pp01_c)
    #     if(Mis_p==1) YY[i] <- 0
    #   }
    #   
    #   if(YY_t[i]==0){
    #     Mis_p <- rbinom(1,size=1,prob=pp10_c)
    #     if(Mis_p==1) YY[i] <- 1
    #   }
    # }
    
    p_Yc <- pp10*(1-p_Y[j])+(1-pp01)*p_Y[j]
    
    YY[j] <-rbinom(1, size=1, prob=p_Yc)
  }
  
  
#matlabのコードを再現
  
  #input
  # y: ２値アウトカム(n x 1)
  # X: 共変量行列(n x p)
  # gamma: パラメータ(1 x 1)
  # b1: 初期値(p x 1)
  
  #output
  # b1: 推定量(p x 1)
  # wi: 重み(n x 1)
  # cov: 漸近共分散行列(p x p)
  
  gamma_logistic_nt <- function(y, X, gamma, b1){
    
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    lam = log(p)/n^1.5
    
    p <- p - 1
    DI <- lam*diag(p+1)
    # DI[1,1]=0
    
    val1 <- obj(y, X, gamma, lam, b1)
    d <- 1; sg <- 1; itr <-100
    
    while (d > 10^(-3) && sg <= itr) {
      b0 <- b1
      val0 <- val1
      ei_gamma <- exp((gamma+1)*(X%*%b0))
      wi <- ( (ei_gamma^y)/(1+ei_gamma) )^(gamma/(gamma+1))
      pi <- ei_gamma/(1+ei_gamma)
      vi <- pi*(1-pi)
      fi <- ((1+ei_gamma)/(1+exp(X%*%b1))^(gamma+1))^(1/(gamma+1))
      
      A <- (t(X)%*%diag(as.vector(fi*vi/n))%*%X + DI)
      b <- (t(X)%*%diag(as.vector(wi/n)) %*% (y-pi)-DI%*%b0)
      HG <- solve(A, b)
      
      val1 <- amj(y, X, gamma, lam, val0, b0, HG)$val1
      b1 <- amj(y,X, gamma, lam, val0, b0, HG)$b1
      
      d <- norm(matrix(b1-b0))/norm(matrix(b0))
      sg <- sg+1
    }
    
    if(sg > itr){
      print(paste('Not converge at gamma (nt) = ', toString(gamma)))
    }
    
    # asymptotic covariance
    ei_gamma <- exp((gamma+1)*(X%*%b1))
    wi <- ( (ei_gamma^y)/(1+ei_gamma) )^(gamma/(gamma+1))
    pi <- ei_gamma/(1+ei_gamma)
    vi <- pi*(1-pi)
    fi <- ((1+ei_gamma)/(1+exp(X%*%b1))^(gamma+1))^(1/(gamma+1))
    delta <- gamma*(t(X)%*%diag(as.vector(wi*(vi-(y-pi)^2) / n))%*%X)
    
    U1 <- (t(X)%*%diag(as.vector((wi*(y-pi))^2 / n)) %*%X)
    H2 <- (t(X)%*%diag(as.vector(fi*vi / n)) %*%X) + delta + DI
    cov = solve(H2, U1/H2 /n)
    
    return(list(b1=b1, wi=wi, cov=cov))
  }
  
  
obj <- function(y, X, gamma, lam, beta){
  ei_gamma <- exp((gamma+1)*X%*%beta)
  val <- mean( ((ei_gamma^y)/(1+ei_gamma))^(gamma/(gamma+1)) )/gamma - 0.5*lam*norm(matrix(beta))^2
  return(val)
}

amj <- function(y, X, gamma, lam, val0, b0, HG){
  rate <- 1
  amj <- 1
  amj_num <- 30
  
  b1 <- b0 + rate * HG
  val1 <- obj(y, X, gamma, lam, b1)
  while(val1 < val0 - 10^-8 && amj <= amj_num){
    b1 <- b0 + (0.5^amj * rate) * HG
    val1 <- obj(y, X, gamma, lam, b1)
    amj <- amj + 1
  }
  if(amj > amj_num){
    print('gamma-logi (nt): armijo limit')
  }
  return(list(val1=val1, b1=b1))
}
    
  
  # 推定方程式
  estimate_eq <- function(beta){
    # beta_0 <- beta[1]
    # beta_t <- beta[2]
    beta_x <- beta[1:9]
    
    # or_ml <- beta_0 + TT*beta_t + XX%*%beta_x
    or_ml <- XX%*%beta_x
    lp <- exp((gamma+1)*(or_ml))
    lq <- exp(YY*(gamma+1)*(or_ml))
    pi <- lp/(1+lp)
    g <- (lq/(1+lp))^(gamma/(1+gamma))
    # w <- ifelse(TT==1, 1/ps_fit, 1/(1-ps_fit))
    
    # return(c(t(g*(YY-pi))%*%rep(1,n), t(g*(YY-pi))%*%TT, t(g*(YY-pi))%*%XX))
    return(t(g*(YY-pi))%*%XX)
  }
  
  gamma <- 2
  b1 <- rep(1,9)
  res_df <- data.frame(t(gamma_logistic_nt(YY, XX, gamma, b1)$b1))
  colnames(res_df) <- c("b1","b2","b3","b4","b5","b6","b7","b8","b9")
  
  results.beta_hat[i,] <- res_df
  
  # gamma <- 2
  # 
  # YY <- rep(0,n)
  # 
  # g <- function(gamma){
  #   beta <- c(0.2,1,1,1,1,1)
  #   beta_0 <- beta[1]
  #   beta_t <- beta[2]
  #   beta_x <- beta[3:6]
  #   
  #   # or_ml <- beta_0 + TT*beta_t + XX%*%beta_x
  #   
  #   or_ml <- rep(1,n)
  #   lp <- exp((gamma+1)*(or_ml))
  #   lq <- exp(YY*(gamma+1)*(or_ml))
  #   pi <- lp/(1+lp)
  #   g <- (lq/(1+lp))^(gamma/(1+gamma))
  #   return(mean(g))
  # }
  # 
  # # h <- Vectorize(g); curve(h,0,100)
  # curve(Vectorize(g)(x), 0, 10)
  # # plot(g_function,0,30)
  
  # gamma <- 0
  # # パラメータ推定 初期値：c(0,0,0,0)
  # estimate_pra_fn <- function(gamma){
  #   gamma <<- gamma
  #   beta_hat <- nleqslv(c(1,1,1,1,1,1,1,1,1), estimate_eq, method = "Newton")$x
  #   
  #   return(list(beta_hat=beta_hat, gamma=gamma))
  # }
  
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
  # beta_hat <- estimate_pra_fn(2)$beta_hat
  
  # beta_t <- beta_hat
  
  # results.gamma <- append(results.gamma, gamma_hat)
  # results.beta_hat <- append(results.beta_hat, beta_hat)
}

results.beta_hat %>% summary()
# boxplot(results.beta_t)

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
