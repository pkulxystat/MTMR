# Hypothesis-testing,Power Calculation and Sample Size decision
# for multi-reader, multi-test design

# This app includes these following functions for practical use: 
# (1) mtmr_correlation: function to estimate 11 correlatoins using rating data 
# (2) mtmr_power: function to compute the power 
# (3) mtmr_test: function for hypothesis-testing using rating data
# (4) mtmr_samplesize: function to decide the samplesize given requirements of
# siginificant level and power


##Here is a method for practical use
Phi <- function(x, y){
  tmp = (x > y) + (x == y)/2
  return(tmp)
}

  
U_statistic <- function(x_rate, y_rate, detail=TRUE){
  colnum <- ncol(x_rate)
  m <- nrow(x_rate)
  n <- nrow(y_rate)
  theta = rep(0, colnum)
  S_10 <- matrix(rep(0, colnum**2), ncol=colnum)
  S_01 <- matrix(rep(0, colnum**2), ncol=colnum)
  S_11 <- matrix(rep(0, colnum**2), ncol=colnum)
  for (k in 1:colnum){
    for (i in 1:m){
      theta[k] <- theta[k] + sum(Phi(x_rate[i,k], y_rate[,k]))
    }
  }
  theta <- theta/m/n
  if (detail == FALSE){return(theta)}
  for (i in 1:colnum){
    for (j in 1:colnum){
      V_10_i <- rep(0, m)
      V_01_i <- rep(0, n)
      V_10_j <- rep(0, m)
      V_01_j <- rep(0, n)
      for (l in 1:m){
        V_10_i[l] <- sum(Phi(x_rate[l,i], y_rate[,i]))
        V_10_j[l] <- sum(Phi(x_rate[l,j], y_rate[,j]))
      }
      for (v in 1:n){
        V_01_i[v] <- sum(Phi(x_rate[,i], y_rate[v,i]))
        V_01_j[v] <- sum(Phi(x_rate[,j], y_rate[v,j]))
      }
      V_10_i <- V_10_i/n
      V_01_i <- V_01_i/m
      V_10_j <- V_10_j/n
      V_01_j <- V_01_j/m
      S_10[i,j] <- sum((V_10_i-theta[i])*(V_10_j-theta[j]))
      S_01[i,j] <- sum((V_01_i-theta[i])*(V_01_j-theta[j]))
      for (v in 1:n){
        S_11[i,j] <- S_11[i,j] + (Phi(x_rate[,i], y_rate[v,i])-theta[i])%*%(Phi(x_rate[,j], y_rate[v,j])-theta[j])
      }
    }
  }
  S_11 <- S_11/(m-1)/(n-1)
  return(list(theta, S_10/(m-1), S_01/(n-1), S_11))
}


#correlation estimate
get_rho <- function(x_rate, y_rate, h, theta, S_10, S_01, S_11){
  m <- nrow(x_rate)
  n <- nrow(y_rate)
  colnum <- ncol(x_rate)
  V_k <- matrix(0, colnum)
  p <- 0
  for (k in 1:colnum)
  {
    sum <- 0
    for (i in 1:m){
      for(j in 1:n){
        if (x_rate[i,k]==y_rate[j,k]){sum <- sum+1}
      }
    }
    prob_eq <- sum/m/n                    
    V_k[k] <- theta[k] -0.25*prob_eq- theta[k]^2
    p <- p + prob_eq
  }  
  p <- p/colnum
  each_column <- colnum/h
  tmp11 <- 0; tmp12 <- 0; tmp13 <- 0; tmp14 <- 0
  tmp21 <- 0; tmp22 <- 0; tmp23 <- 0; tmp24 <- 0
  tmp31 <- 1; tmp32 <- 0; tmp33 <- 0; tmp34 <- 0
  for (k1 in 1:colnum){
    for (k2 in 1:colnum){
      division <- sqrt(V_k[k1]*V_k[k2])
      if (division==0){division <- 1}
      if (k1==k2){
        tmp11 <- tmp11 + S_10[k1, k2]/division/colnum
        tmp21 <- tmp21 + S_01[k1, k2]/division/colnum
      }
      else if((k1-1)%/%each_column == (k2-1)%/%each_column){
        tmp12 <- tmp12 + S_10[k1, k2]/division/colnum/(each_column-1)
        tmp22 <- tmp22 + S_01[k1, k2]/division/colnum/(each_column-1)
        tmp32 <- tmp32 + S_11[k1, k2]/division/colnum/(each_column-1)
      }
      else if((k1-k2)%%each_column==0){
        tmp13 <- tmp13 + S_10[k1, k2]/division/colnum/(h-1)
        tmp23 <- tmp23 + S_01[k1, k2]/division/colnum/(h-1)
        tmp33 <- tmp33 + S_11[k1, k2]/division/colnum/(h-1)
      }
      else {
        tmp14 <- tmp14 + S_10[k1, k2]/division/colnum/(h-1)/(each_column-1)
        tmp24 <- tmp24 + S_01[k1, k2]/division/colnum/(h-1)/(each_column-1)
        tmp34 <- tmp34 + S_11[k1, k2]/division/colnum/(h-1)/(each_column-1)
      }
    }
  }
  return(c(tmp11,tmp12,tmp13,tmp14,tmp21,tmp22,tmp23,tmp24,tmp31,tmp32,tmp33,tmp34,p))
}


var_estimate <- function(rho_exp, m, n, r, V, p=0){
  tmp <- (n-1) * rho_exp[1] + (m-1) * rho_exp[5] + rho_exp[9]
  tmp <- tmp + (r-1) * ((n-1) * rho_exp[2] + (m-1) * rho_exp[6] + rho_exp[10])
  tmp <- tmp - (n-1) * rho_exp[3] - (m-1) * rho_exp[7] - rho_exp[11]
  tmp <- tmp - (r-1) * ((n-1) * rho_exp[4] + (m-1) * rho_exp[8] + rho_exp[12])
  return(tmp*2*(V-V**2-p/4)/(m*n*r))
}


make_decision <- function(x_rate, y_rate){
  w_1 <- rep(1/m, m)
  w_2 <- matrix(rep(1/r, r), ncol=1)
  a11 <- var(c(w_1 %*% x_rate[,1:r]))
  a12 <- var(c(x_rate[,1:r] %*% w_2))
  b11 <- var(c(w_1 %*% x_rate[,(r+1):(2*r)]))
  b12 <- var(c(x_rate[,(r+1):(2*r)] %*% w_2))
  w_1 <- rep(1/n, n)
  w_2 <- matrix(rep(1/r, r), ncol=1)
  a21 <- var(c(w_1 %*% y_rate[,1:r]))
  a22 <- var(c(y_rate[,1:r] %*% w_2))
  b21 <- var(c(w_1 %*% y_rate[,(r+1):(2*r)]))
  b22 <- var(c(y_rate[,(r+1):(2*r)] %*% w_2))
  if ((a11+a21)/(a12+a22) > 0.1  | (b11+b21)/(b12+b22) > 0.1){return('bootstrap')}
  else{return('simple')}
}


bootstrap <- function(x_rate, y_rate, times, r){
  auc_1 <- c()
  auc_2 <- c()
  m <- nrow(x_rate)
  n <- nrow(y_rate)
  diff_auc <- c()
  #select readers randomly with return
  for (simulation in 1:times){
    x_part1 <- c(); x_part2 <- c()
    y_part1 <- c(); y_part2 <- c()
    x_index <- c()
    y_index <- c()
    for (j in 1:m){
      x_index <- append(x_index, sample(c(1:m),1))
    }
    for (j in 1:n){
      y_index <- append(y_index, sample(c(1:n),1))
    }
    for (j in 1:r){
      index <- sample(c(1:r),1) 
      x_part1 <- cbind(x_part1, x_rate[x_index,index])
      x_part2 <- cbind(x_part2, x_rate[x_index,index+r])
      y_part1 <- cbind(y_part1, y_rate[y_index,index])
      y_part2 <- cbind(y_part2, y_rate[y_index,index+r])
    }
    x_data <- cbind(x_part1, x_part2)
    y_data <- cbind(y_part1, y_part2)
    theta_hat <- U_statistic(x_data, y_data, detail=FALSE)  #get estimation of AUCs
    theta_hat_h <- rep(0, h)
    for (i in 1:h){
      theta_hat_h[i] <- sum(theta_hat[((i-1)*r+1):(i*r)])/r
    }
    auc_1 <- append(auc_1, theta_hat_h[1])
    auc_2 <- append(auc_2, theta_hat_h[2])
    diff_auc <- append(diff_auc, theta_hat_h[1] - theta_hat_h[2])
  }
  sd_auc_1 <- sqrt(var(auc_1))
  sd_auc_2 <- sqrt(var(auc_2))
  sd_diff_auc <- sqrt(var(diff_auc))
  return(c(sd_auc_1, sd_auc_2, sd_diff_auc))
}


##estimate the correlation from data
mtmr_correlation <- function(x_rate, y_rate, h=2){
  U_result <- U_statistic(x_rate, y_rate)
  theta_hat <- U_result[[1]]    #get wilcoxon estimation of AUCs
  #print(theta_hat[1])
  S_10 <- U_result[[2]]
  S_01 <- U_result[[3]]
  S_11 <- U_result[[4]]
  rho_exp <- get_rho(x_rate, y_rate, h, theta_hat,S_10,S_01,S_11)
  return(rho_exp)
}


##power calculation
mtmr_power <- function(r, m, n, AUC, delta, level=0.95, rho_exp, equal_p=0){
  z_norm <- qnorm((1-level)/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  var_theta_h <- var_estimate(rho_exp, m, n, r, AUC, equal_p)
  power <- pnorm(z_norm + delta/sqrt(var_theta_h), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+
    pnorm(z_norm - delta/sqrt(var_theta_h), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  return(power)
}


##hypothesis-testing
mtmr_test <- function(x_rate, y_rate, level=0.95, times = 20000, h=2, fixed_reader=TRUE){
  m <- nrow(x_rate)
  n <- nrow(y_rate)
  r <- ncol(x_rate)/h
  if (fixed_reader){method='simple'}else {method <- make_decision(x_rate, y_rate)}
  z_norm <- qnorm((1-level)/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  if (method == 'simple'){
    U_result <- U_statistic(x_rate, y_rate)
    theta_hat <- U_result[[1]]    #get wilcoxon estimation of AUCs
    S_10 <- U_result[[2]]
    S_01 <- U_result[[3]]
    S_11 <- U_result[[4]]
    S <- S_10/m + S_01/n 
    theta_hat_h <- rep(0, h)
    for (i in 1:h){
      theta_hat_h[i] <- sum(theta_hat[((i-1)*r+1):(i*r)])/r
    } 
    
    C <- c();C_1 <- c();C_2 <- c()
    C[1:r] <- 1/r; C_1[1:r] <- 1/r; C_2[1:r] <- 0
    C[(r+1):(2*r)] <- -1/r; C_1[(r+1):(2*r)] <- 0; C_2[(r+1):(2*r)] <- 1/r;
    sd_1 <- sqrt(C_1 %*% S %*% C_1)
    sd_2 <- sqrt(C_2 %*% S %*% C_2)
    sd_diff <- sqrt(C %*% S %*% C)
    interval_1 <- c(theta_hat_h[1]+z_norm*sd_1, theta_hat_h[1]-z_norm*sd_1)
    interval_2 <- c(theta_hat_h[2]+z_norm*sd_2, theta_hat_h[2]-z_norm*sd_2)
    interval_diff <- c(theta_hat_h[1]-theta_hat_h[2]+z_norm*sd_diff, theta_hat[1]-theta_hat_h[2]-z_norm*sd_diff)
    Z <- (C %*% theta_hat)/sd_diff  #test for the equality of AUCs
    if (Z < z_norm){result <- 'method 2 is better'}
    else if (Z > -z_norm){result <- 'method 1 is better'}
    else{result <- 'no difference'}
  }
  else{
    theta_hat <- U_statistic(x_rate, y_rate, detail=FALSE)
    theta_hat_h <- rep(0, h)
    for (i in 1:h){
      theta_hat_h[i] <- sum(theta_hat[((i-1)*r+1):(i*r)])/r
    } 
    bootstrap_result <- bootstrap(x_rate, y_rate, times, r)
    sd_1 <- bootstrap_result[1]
    sd_2 <- bootstrap_result[2]
    sd_diff <- bootstrap_result[3]
    interval_1 <- c(theta_hat_h[1]+z_norm*sd_1, theta_hat_h[1]-z_norm*sd_1)
    interval_2 <- c(theta_hat_h[2]+z_norm*sd_2, theta_hat_h[2]-z_norm*sd_2)
    interval_diff <- c(theta_hat_h[1]-theta_hat_h[2]+z_norm*sd_diff, theta_hat_h[1]-theta_hat_h[2]-z_norm*sd_diff)
    Z <- (theta_hat[1]-theta_hat_h[2])/sd_diff  #test for the equality of AUCs
    if (Z < z_norm){result <- 'method 2 is better'}
    else if (Z > -z_norm){result <- 'method 1 is better'}
    else{result <- 'no difference'}
  }
  return(list(result,interval_diff, interval_1,interval_2))
}


##sample-size decision, ratio is the proportion of diseased patients
mtmr_samplesize <- function(r, AUC, delta, power, rho_exp, equal_p=0, ratio='best', level=0.95){
  z_norm <- qnorm((1-level)/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  V <- AUC - AUC^2 - equal_p/4
  var_tolerance <- (delta/(qnorm(power, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) - z_norm))^2
  C_0 <- 1 - rho_exp[1] - rho_exp[5] + (r-1) * (rho_exp[10]-rho_exp[2]-rho_exp[6]) - rho_exp[11] +
    rho_exp[3] + rho_exp[7] - (r-1) * (rho_exp[12]-rho_exp[4]-rho_exp[8])
  if (ratio == 'best'){
    C_1 <- rho_exp[1] + (r-1) * rho_exp[2] - rho_exp[3] - (r-1) * rho_exp[4]
    C_2 <- rho_exp[5] + (r-1) * rho_exp[6] - rho_exp[7] - (r-1) * rho_exp[8]
    tmp1 <- r * var_tolerance/8/V
    tmp2 <- -(C_1 + C_2)/2
    tmp3 <- V*(C_1-C_2)^2/2/r/var_tolerance - C_0
    tmp <- sqrt(tmp2^2 - 4*tmp1*tmp3)
    N <- (-tmp2+tmp)/2/tmp1
    ratio <- 1/2 + V*(C_1-C_2)/N/r/var_tolerance
  }
  C_1 <- (1-ratio) * rho_exp[1] + ratio * rho_exp[5] + (r-1) * ((1-ratio) * rho_exp[2] + ratio * rho_exp[6])-
      (1-ratio) * rho_exp[3] - ratio * rho_exp[7] - (r-1) * ((1-ratio) * rho_exp[4] + ratio * rho_exp[8])
  C_2 <- r * ratio * (1-ratio) * var_tolerance/2/V
  tmp <- sqrt(C_1^2 + 4*C_0*C_2)
  N <- ceiling((C_1 + tmp)/2/C_2)
  m <- round(N*ratio)
  n <- N-m
  return(c(N, m, n))
}

#An example for practical use:
#AUC=0.825;r=10;delta=0.05;power=0.9;rho_exp=c(0.35709644,0.08953128,0.27397317,
#0.07368219,0.25279248,0.06770141,0.19881069,0.05571047,1.00000000,0.16964908,0.60065536,
#0.13833997); equal_p=0; m=50; n=50
#mtmr_samplesize(r, AUC, delta, power, rho_exp, ratio='best', level=0.95)
#mtmr_power(r, m, n, AUC, delta, level=0.95, rho_exp, equal_p)
