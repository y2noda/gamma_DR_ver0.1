library(tidyverse)
library(MASS)

#input
# y: ２値アウトカム(n x 1)
# X: 共変量行列(n x p) 一列目が切片　
# TT: 治療変数(n x 1)
# gamma: パラメータ(1 x 1)
# b1: 初期値(p x 1)

#output
# b1: 推定量(p x 1)
# wi: 重み(n x 1)
# cov: 漸近共分散行列(p x p)

gammat_logistic_nt = function(y, X, TT, gamma, b1){
  
  # 傾向スコアの算出
  ps_fit <- glm(TT~X, family=binomial)$fit
  
  # 共変量行列の一列目に治療変数を追加
  X <- cbind(TT,X)
  
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
    
    # MSMのためのウェイト
    ps_w <- ifelse(TT==1, 1/ps_fit, 1/(1-ps_fit))
    # ps_w <- rep(1,n)
    
    #||f||_(gamma+1) (n x 1)
    fi <- ((1+ei_gamma)/(1+exp(X%*%b1))^(gamma+1))^(1/(gamma+1))
    
    #ヘッセ行列 (p x p)
    H <- (t(X)%*%diag(as.vector(ps_w*fi*vi/n))%*%X + DI)
    
    #スコア関数 (p x 1)
    S <- (t(X)%*%diag(as.vector(ps_w*wi/n)) %*% (y-pi)-DI%*%b0)
    
    #変化量H^(-1)*S
    HG <- solve(H, S)
    
    val1 <- amj(y, X, gamma, lam, val0, b0, HG)$val1
    b1 <- amj(y, X, gamma, lam, val0, b0, HG)$b1
    
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


obj = function(y, X, gamma, lam, beta){
  ei_gamma <- exp((gamma+1)*X%*%beta)
  val <- mean( ((ei_gamma^y)/(1+ei_gamma))^(gamma/(gamma+1)) )/gamma - 0.5*lam*norm(matrix(beta))^2
  return(val)
}

amj = function(y, X, gamma, lam, val0, b0, HG){
  rate <- 0.5
  amj <- 1
  amj_num <- 30
  
  b1 <- b0 + rate * HG
  val1 <- obj(y, X, gamma, lam, b1)
  while(val1 < val0 - 10^(-8) && amj <= amj_num){
    b1 <- b0 + (0.5^amj * rate) * HG
    val1 <- obj(y, X, gamma, lam, b1)
    amj <- amj + 1
  }
  if(amj > amj_num){
    print('gamma-logi (nt): armijo limit')
  }
  return(list(val1=val1, b1=b1))
}

