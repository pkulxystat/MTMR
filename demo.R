source('D:/tmp/mtmr_app.R', encoding = 'UTF-8')   #请将路径替换为mtmr_app.R保存在本地的路径

get_score <- function(scores, low, high){
  tmp <- exp(sum(log((high - scores)/(high-low))))
  score <- high * (1-tmp) + low * tmp
  return(score)
}

raw_data <- read.csv('D:/tmp/cancer_data.csv', header=TRUE)  #请将路径替换为cancer_data.csv保存在本地的路径
raw_data <- as.matrix(raw_data)
r <- 6   #six readers
partition <- 10  #10个检测区域
low <- 0; high <- 100  #打分范围
h <- 2   #影像方法数量
N <- nrow(raw_data)/r
m <- 0


for (i in 1:N){
  if (sum(raw_data[i,1:partition]==1)>0){
    m <- m+1
  }
}
n <- N-m

x_rate <- matrix(rep(0, m*r*h), nrow=m)
y_rate <- matrix(rep(0, n*r*h), nrow=n)
i <- 0  #which benign case
j <- 0  #which malign case

for (index in 1:nrow(raw_data)){
  k <- (index-1) %/% N + 1  #which reader
  if (sum(raw_data[index,1:partition]==1)>0){   ##malign case
    i <- i %% m +1
    for (l in 1:h){
      tmp <- raw_data[index, (partition * l+1):(partition * (l+1))]
      score <- get_score(tmp, low, high)
      x_rate[i, (l-1)*r + k] <- score
    }
  }
  else {
    j <- j %% n + 1
    for (l in 1:h){
      tmp <- raw_data[index, (partition * l+1):(partition * (l+1))]
      score <- get_score(tmp, low, high)
      y_rate[j, (l-1)*r + k] <- score
    }
  }
}

#here is the part of analysis tool
fix_result <- mtmr_test(x_rate, y_rate, level=0.95, times = 20000, h=2, fixed_reader=TRUE)  #不考虑随机选择医生带来的随机性
not_fix_result <- mtmr_test(x_rate, y_rate, level=0.95, times = 20000, h=2, fixed_reader=FALSE)   #不考虑随机选择医生带来的随机性

coefs <- mtmr_correlation(x_rate, y_rate)   #估计数据对应的相关系数和正负样本相等的概率
AUC <- (mean(not_fix_result[[3]]) + mean(not_fix_result[[4]]))/2   #平均AUC估计
rho_exp <- coefs[1:12]
equal_p <- coefs[13]

delta <- 0.05   #区分AUC时的差值
power <- 0.9    #检验功效的要求
level <- 0.95   #置信度
true_fix_power <- mtmr_power(r, m, n, AUC, delta, level=0.95, rho_exp, equal_p)   #预计功效

best_plan <- mtmr_samplesize(r, AUC, delta, power, rho_exp, equal_p, ratio='best', level=0.95)   #fix reader effect对应的正负样本分配方案

