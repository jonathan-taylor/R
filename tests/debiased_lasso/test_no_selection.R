library(selectiveInference)
source("tests/debiased_lasso/test_debiased_coverage.R")


# comparing the length of debiased lasso and OLS without selection

test_methods = function(n, p, s, rho, nrep=50){
  
  snr = sqrt(2*log(p)/n)
  
  coverages_X=NULL
  coverages_popX = NULL
  coverages_dl = NULL
  length_X=NULL
  length_popX=NULL
  length_dl = NULL
  
  for (i in 1:nrep){
    times = c(1:p)
    Sigma <- rho^abs(outer(times, times, "-"))
    data = selectiveInference:::gaussian_instance(n=n, p=p, s=s, rho=rho, sigma=1, snr=snr)
    popX = n*Sigma # E[X^TX]
    
    X=data$X
    y=data$y
    beta=data$beta
    
    lmX=lm(y~X-1)
    int_X = confint(lmX, level=0.9)  # standard OLS intervals
    print(dim(int_X))
    sigma_est = sigma(lmX)
    print(c("sigma est:", sigma_est))
    
    
    popX_inverse = solve(popX)
    est_popX = popX_inverse %*% t(X) %*% y
    mat = diag(rep(n*sum(beta^2),p))
    print(dim(mat))
    var_popX = popX_inverse*sigma_est^2+popX_inverse %*% mat %*% popX_inverse
    int_popX = matrix(nrow=p, ncol=2)
    for (i in 1:p){
      int_popX[i,] = selectiveInference:::naive_CI(est_popX[i], var_popX[i,i], alpha=0.1)
    }
   
    
    lambda_frac=0.8
    loss="ls"
    lambda = lambda_frac*selectiveInference:::theoretical.lambda(X, loss, sigma_est)  # theoretical lambda
    print(c("lambda", lambda))
    soln = selectiveInference:::solve_problem_glmnet(X, y, lambda, penalty_factor=rep(1,p), loss=loss)
    print(c("nactive", length(which(soln!=0))))
    cat("selected", which(soln!=0), "\n")
    results_dl = debiased_lasso_inference(X, y, soln=soln, loss="ls", sigma_est=sigma_est, lambda=lambda, debias_mat="JM")
    int_dl = t(results_dl$naive_intervals)
    
    
    coverages_X=c(coverages_X, selectiveInference:::compute_coverage(t(int_X), beta))
    length_X = c(length_X, as.vector(int_X[,2]-int_X[,1]))
    coverages_popX=c(coverages_popX, selectiveInference:::compute_coverage(t(int_popX), beta))
    length_popX = c(length_popX, as.vector(int_popX[,2]-int_popX[,1]))
    coverages_dl=c(coverages_dl, selectiveInference:::compute_coverage(t(int_dl), beta))
    length_dl = c(length_dl, as.vector(int_dl[,2]-int_dl[,1]))
    
    print(c("X coverage:",mean(coverages_X)))
    print(c("X length: ",mean(length_X)))
    print(c("pop X coverage:",mean(coverages_popX)))
    print(c("pop X length: ", mean(length_popX)))
    print(c("ratio:", mean(length_X)/mean(length_popX)))
    print(c("predicted ratio:", sqrt(n/(n-p-1))))
    
    print(c("DL coverage:",mean(coverages_dl)))
    
    print(c("DL length: ",mean(length_dl)))
    
  }
}

test_methods(100, 80, 30, 0.)


