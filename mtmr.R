library(MASS)

#the dimension of the vector v is 4
get_cov <- function(v, r, h){
    result <- matrix(rep(0, (h * r) ** 2), nrow=h*r)
    for (i in 1:(h * h)) 
    {
        row <- floor((i-1)/h) + 1
        col <- i - h * (row-1)
        if (row == col)
        {
          tmp <- matrix(rep(v[2], r*r), nrow = r)+ diag(r)* (v[1] - v[2])
        }
        else
        {
          tmp <- matrix(rep(v[4], r*r), nrow = r)+ diag(r)* (v[3] - v[4])
        }
        result[((row-1)*r+1):(row*r),((col-1)*r+1):(col*r)] <- tmp
    }
    return(result)
}

#m:diseased, n:non-diseased, h:modalities, r:readers, mu:h-dimension vector, 
generate_data <- function(m, n, h, r, sigma_r, mu, cov_coef){
    #r_x <- rnorm(r*h, sd=sqrt(sigma_r[1]))  ##random readers
    #r_y <- rnorm(r*h, sd=sqrt(sigma_r[2]))
    r_x <- fix_r_x  ##fixed readers
    r_y <- fix_r_y
    cov_x <- get_cov(cov_coef[1,], r, h)
    cov_y <- get_cov(cov_coef[2,], r, h)
    mu_x <- c()
    for (i in 1:h){
      mu_x <- append(mu_x, rep(mu[i],r))
    }
    epsilon_x <- mvrnorm(m, mu_x, cov_x)
    epsilon_y <- mvrnorm(n, rep(0,r*h), cov_y)
    x_rate <- t(r_x + t(epsilon_x))
    y_rate <- t(r_y + t(epsilon_y))
    return (list(x_rate, y_rate))
}


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
    p <- p + prob_eq
    V_k[k] <- theta[k] -0.25*prob_eq- theta[k]^2
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

power_simulation <- function(m, n, h, r, sigma_r, cov_coef, mu, all, level=0.95){
  different <- 0
  true_theta_hat_h <- rep(0, h)
  rho_exp <- rep(0, 13)
  for (simulation in 1:all){
    rate <- generate_data(m, n, h, r, sigma_r, mu, cov_coef)
    x_rate <- rate[[1]]   #diseased
    y_rate <- rate[[2]]   #non-diseased
    U_result <- U_statistic(x_rate, y_rate)
    theta_hat <- U_result[[1]]    #get wilcoxon estimation of AUCs
    #print(theta_hat[1])
    S_10 <- U_result[[2]]
    S_01 <- U_result[[3]]
    S_11 <- U_result[[4]]
    S <- S_10/m + S_01/n 
    theta_hat_h <- rep(0, h)
    for (i in 1:h){
      theta_hat_h[i] <- sum(theta_hat[((i-1)*r+1):(i*r)])/r
      true_theta_hat_h[i] <- ((simulation - 1) * true_theta_hat_h[i] + theta_hat_h[i])/simulation
    } 
    C <- c()
    C[1:r] <- 1/r
    C[(r+1):(2*r)] <- -1/r
    Z <- (C %*% theta_hat)/sqrt(C %*% S %*% C)  #test for the equality of AUCs
    if (abs(Z) > qnorm(1.95/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)){
      different <- different + 1
    }
    rho_exp <- (get_rho(x_rate,y_rate,2,theta_hat,S_10,S_01,S_11) + rho_exp * (simulation - 1))/simulation
  }
  #rho_exp <- c(0.31, 0.08, 0.24, 0.06, 0.22, 0.06, 0.17, 0.05, 1, 0.15, 0.55, 0.12, 0)
  #rho_exp <- c(0.528, 0.243, 0.267, 0.213, 0.24, 0.102, 0.118, 0.096, 1, 0.382, 0.431, 0.341, 0)
  V <- mean(true_theta_hat_h)
  diff_auc <- true_theta_hat_h[1] - true_theta_hat_h[2]
  var_theta_h <- var_estimate(rho_exp[1:12], m, n, r, V, rho_exp[13])
  z_norm <- qnorm((1-level)/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  delta <- true_theta_hat_h[1] - true_theta_hat_h[2]
  power <- pnorm(z_norm + delta/sqrt(var_theta_h), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+
    pnorm(z_norm - delta/sqrt(var_theta_h), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  return(list(true_theta_hat_h, power, different/all))  #auc_estimation, theoritical power, empirical power
}

##here is a sample for simulation(1000 times) 
m <- 67
n <- 133
h <- 2
r <- 15
sigma_r <- c(0.02^2, 0.03^2)  #small variance for the reader effect
cov_coef <- matrix(c(0.98, 0.3, 0.8, 0.25, 0.72, 0.225, 0.6, 0.1875), byrow= T,nrow = 2) #ncol = 4
####Power calculation
mu <- c(1.37, 1.12)
all <- 1000
level <- 0.95 #significant level
delta = 0.05   #the tolerance of the difference of AUCs
set.seed(12345)
fix_r_x <- rnorm(r*h, sd=sqrt(sigma_r[1]))
set.seed(52431)
fix_r_y <- rnorm(r*h, sd=sqrt(sigma_r[2]))

power_result <- power_simulation(m, n, h, r, sigma_r, cov_coef, mu, all, level=0.95)
true_theta_hat_h <- power_result[[1]]
diff_auc <- true_theta_hat_h[1] - true_theta_hat_h[2]
print(power_result)
print(diff_auc)

#check the estimation of auc
valid <- 0
valid_1 <- 0
valid_2 <- 0
diff_record <- c()
sd_record <- c()
z_norm <- qnorm((1-level)/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
for (simulation in 1:all){
  rate <- generate_data(m, n, h, r, sigma_r, mu, cov_coef)
  x_rate <- rate[[1]]   #diseased
  y_rate <- rate[[2]]   #non-diseased
  U_result <- U_statistic(x_rate, y_rate)
  theta_hat <- U_result[[1]]    #get wilcoxon estimation of AUCs
  S_10 <- U_result[[2]]
  S_01 <- U_result[[3]]
  S_11 <- U_result[[4]]
  S <- S_10/m + S_01/n 
  C <- c();C_1 <- c();C_2 <- c()
  C[1:r] <- 1/r; C_1[1:r] <- 1/r; C_2[1:r] <- 0
  C[(r+1):(2*r)] <- -1/r; C_1[(r+1):(2*r)] <- 0; C_2[(r+1):(2*r)] <- 1/r;
  sd_1 <- sqrt(C_1 %*% S %*% C_1)
  sd_2 <- sqrt(C_2 %*% S %*% C_2)
  sd_diff <- sqrt(C %*% S %*% C)
  theta_hat_h <- rep(0, h)
  for (i in 1:h){
    theta_hat_h[i] <- sum(theta_hat[((i-1)*r+1):(i*r)])/r
  } 
  if (theta_hat_h[1]-theta_hat_h[2] + z_norm * sd_diff < diff_auc & theta_hat_h[1]-theta_hat_h[2] - z_norm * sd_diff > diff_auc){
    valid <- valid + 1
  }
  sd_record <- append(sd_record, sd_diff)
  diff_record <- append(diff_record, theta_hat_h[1]-theta_hat_h[2])
  if (theta_hat_h[1] + z_norm * sd_1 < true_theta_hat_h[1] & theta_hat_h[1] - z_norm * sd_1 > true_theta_hat_h[1]){
    valid_1 <- valid_1 + 1
  }
  if (theta_hat_h[2] + z_norm * sd_2 < true_theta_hat_h[2] & theta_hat_h[2] - z_norm * sd_2 > true_theta_hat_h[2]){
    valid_2 <- valid_2 + 1
  }
}
print(valid/all)
print(valid_1/all)
print(valid_2/all)
print(sqrt(var(diff_record)))
print(mean(sd_record))
qqnorm(diff_record)

##test for significant level
mu <- c(1.12, 1.12)
different <- 0
all <- 1000
true_theta_hat_h <- rep(0, h)
rho_exp <- rep(0, 13)
for (simulation in 1:all){
  rate <- generate_data(m, n, h, r, sigma_r, mu, cov_coef)
  x_rate <- rate[[1]]   #diseased
  y_rate <- rate[[2]]   #non-diseased
  U_result <- U_statistic(x_rate, y_rate)
  theta_hat <- U_result[[1]]    #get estimation of AUCs
  S_10 <- U_result[[2]]
  S_01 <- U_result[[3]]
  S_11 <- U_result[[4]]
  S <- S_10/m + S_01/n + S_11/m/n
  theta_hat_h <- rep(0, h)
  for (i in 1:h){
    theta_hat_h[i] <- sum(theta_hat[((i-1)*r+1):(i*r)])/r
    true_theta_hat_h[i] <- ((simulation - 1) * true_theta_hat_h[i] + theta_hat_h[i])/simulation
  } 
  C <- c()
  C[1:r] <- 1/r
  C[(r+1):(2*r)] <- -1/r
  Z <- (C %*% theta_hat)/sqrt(C %*% S %*% C)  #test for the equality of AUCs
  if (abs(Z) > qnorm(1.95/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)){
    different <- different + 1
  }
  rho_exp <- (get_rho(x_rate,y_rate,2,theta_hat,S_10,S_01,S_11) + rho_exp * (simulation - 1))/simulation
}

print(1 - different/all)  #significant_level

#When sigma_r is large, the method mentioned above is not correct when the readers are not fixed
#Here is the bootstrap method for variance estimation and getting confidence interval

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

#now let the variance caused by the readers become larger
sigma_r <- c(0.01, 0.02) 
mu <- c(1.37, 1.12)
all <- 2500
times <- 30000
true_theta_hat_h <- rep(0, h)
rate <- generate_data(m, n, h, r, sigma_r, mu, cov_coef)
x_rate <- rate[[1]]   #diseased
y_rate <- rate[[2]]   #non-diseased


bootstrap_result <- bootstrap(x_rate, y_rate, times, r)
sd_auc_1 <- bootstrap_result[1]
sd_auc_2 <- bootstrap_result[2]
sd_diff_auc <- bootstrap_result[3]

#confidence interval for the aucs
theta_hat <- U_statistic(x_rate, y_rate, detail=FALSE) #get wilcoxon estimation of AUCs
theta_hat_h <- rep(0, h)
for (i in 1:h){
  theta_hat_h[i] <- sum(theta_hat[((i-1)*r+1):(i*r)])/r
} 
interval_diff <- c(theta_hat_h[1]-theta_hat_h[2]+z_norm*sd_diff_auc,
                   theta_hat_h[1]-theta_hat_h[2]-z_norm*sd_diff_auc)
interval_1 <- c(theta_hat_h[1]+z_norm*sd_auc_1,
                theta_hat_h[1]-z_norm*sd_auc_1)
interval_2 <- c(theta_hat_h[2]+z_norm*sd_auc_1,
                theta_hat_h[2]-z_norm*sd_auc_2)

diff_record <- c()
auc_1_record <- c()
auc_2_record <- c()
for (simulation in 1:all){
  rate <- generate_data(m, n, h, r, sigma_r, mu, cov_coef)
  x_rate <- rate[[1]]   #diseased
  y_rate <- rate[[2]]   #non-diseased
  theta_hat <- U_statistic(x_rate, y_rate, detail=FALSE) #get wilcoxon estimation of AUCs
  theta_hat_h <- rep(0, h)
  for (i in 1:h){
    theta_hat_h[i] <- sum(theta_hat[((i-1)*r+1):(i*r)])/r
  } 
  auc_1_record <- append(auc_1_record, theta_hat_h[1])
  auc_2_record <- append(auc_2_record, theta_hat_h[2])
  diff_record <- append(diff_record, theta_hat_h[1]-theta_hat_h[2])
}

##bootstrap method is useful for this case
print(sqrt(var(auc_1_record)))
print(sqrt(var(auc_2_record)))
print(sqrt(var(diff_record)))
print(bootstrap_result)
qqnorm(diff_record) ##test normality


#judge whether the influence of readers is under control 
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

