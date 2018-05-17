library(glmnet)
library(selectiveInference)

setwd("/Users/Jelena/Dropbox/kevin/jelena/crohns")
data=readRDS("real_data.rds")
X=data$X
y=data$y

X=scale(X)
y=y-mean(y)

set.seed(1)
p=ncol(X)
loss="ls"

#CV = cv.glmnet(X, y, standardize=FALSE, intercept=FALSE, family="gaussian")
#sigma_est=selectiveInference:::estimate_sigma(X,y,coef(CV, s="lambda.min")[-1]) # sigma via Reid et al.
#print(c("sigma est via Reid et al.", sigma_est))   # "sigma est via Reid et al." "0.620021797352246" 

#print(c("std error of y", sqrt(var(y))))  # "std error of y: 0.488729008112197"

estimate_sigma = function(X,y){
  n=nrow(X)
  p=ncol(X)
  loss="ls"
  m=floor(n/2)
  subsample = sample(1:n,m, replace=FALSE)
  leftover = setdiff(1:n, subsample)
  CV = cv.glmnet(X[subsample,], y[subsample], standardize=FALSE, intercept=FALSE, family="gaussian")
  beta_hat = coef(CV, s="lambda.min")[-1]
  #lambda = 0.8*selectiveInference:::theoretical.lambda(X[subsample,], loss, sqrt(var(y)))
  #beta_hat = selectiveInference:::solve_problem_glmnet(X[subsample,], y[subsample], lambda, penalty_factor=rep(1,p), loss=loss)
  
  selected=which(beta_hat!=0)
  nselected=length(selected)
  print(c("nselected",nselected))
  if (nselected==0){
    return(0)
  }
  LM = lm(y[leftover]~X[leftover,][,selected])
  sigma_est = sigma(LM)
  print(sigma_est)
  return(sigma_est) 
}

repeat_data_splitting = function(X,y){
  nreps=10
  sigma_est=0
  nest = 0
  for (b in 1:nreps){
    curr_est = estimate_sigma(X,y)
    sigma_est = sigma_est+curr_est
    if (curr_est>0){
      nest=nest+1
    }
  }

  sigma_est = sigma_est/nest
  cat("sigma est via data splitting:", sigma_est, "\n")  # 0.4753636 
  return(sigma_est)
}

sigma_est = repeat_data_splitting(X,y)

lambda = 0.8*selectiveInference:::theoretical.lambda(X, loss, sigma_est) # 0.02701314
cat("lambda:", lambda, "\n")


beta_hat = selectiveInference:::solve_problem_glmnet(X, y, lambda, penalty_factor=rep(1,p), loss=loss)
selected = which(beta_hat!=0)
nselected=length(selected)
print(c("nselected",nselected))  # 99


