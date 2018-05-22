library(glmnet)
library(selectiveInference)
library(MASS)
library(knockoff)

cluster=FALSE

if (cluster==TRUE){
  setwd(getwd())
  args = commandArgs(trailingOnly=TRUE)
  method = toString(args[1])
  outdir = "/scratch/users/jelenam/full/"
  label = paste(method, "_real_data_results", sep="")
  outfile = file.path(outdir, paste(sep="",label, ".rds"))
} else{
  setwd("/Users/Jelena/Dropbox/kevin/jelena/crohns")
  method="liu"
  outfile="liu_real_data_results.rds"
}

set.seed(1)
loss="ls"
sigma_est = 0.4753636
lambda = 0.02701314


liu_full = function(outfile){
  if (cluster==TRUE){
    data=readRDS("real_data1.rds")
  } else{
    data=readRDS("real_data.rds")
  }
  X=data$X
  y=data$y
  n=nrow(X)
  X=scale(X)
  y=y-mean(y)
  
  penalty_factor = rep(1, ncol(X))
  
  soln = selectiveInference:::solve_problem_glmnet(X, y, lambda, penalty_factor=penalty_factor, loss=loss)
  PVS = selectiveInference:::inference_debiased_full(X, y, soln, lambda=lambda, penalty_factor=penalty_factor, 
                                                 sigma_est, loss=loss, algo="Q", construct_ci = TRUE,
                                                 verbose=TRUE)
  
  cat("nactive:", length(PVS$active_vars), "\n")
  cat("active vars:", PVS$active_vars, "\n")
  
  saveRDS(list(active_vars=PVS$active_vars, 
               sel_intervals=PVS$sel_intervals, naive_intervals=PVS$naive_intervals,
               pvalues=PVS$pvalues, naive_pvalues=PVS$naive_pvalues), file=outfile)
  
  return(NULL)
}



lee = function(outfile, type){
  if (type=="full"){
    data=readRDS("real_data2.rds")
  } else if (type=="partial"){
    data=readRDS("real_data3.rds")
  }
  X=data$X
  y=data$y
  n=nrow(X)
  X=scale(X)
  y=y-mean(y)
  
  lasso = glmnet(X, y, family=selectiveInference:::family_label(loss), alpha=1, standardize=FALSE, intercept=FALSE, thresh=1e-12)
  soln = as.numeric(coef(lasso,x=X,y=y, family=selectiveInference:::family_label(loss), s=lambda, exact=TRUE))[-1]
  PVS = selectiveInference:::fixedLassoInf(X,y,soln, intercept=FALSE, lambda*n, family=selectiveInference:::family_label(loss),
                                           type=type,sigma=sigma_est)
  
  abs_soln = abs(soln)
  beta_threshold = abs_soln[order(abs_soln,decreasing=TRUE)][length(PVS$pv)]
  active_vars = which(abs_soln>=beta_threshold)
  cat("nactive:", length(active_vars), "\n")
  cat("active vars:", active_vars, "\n")
  
  pvalues = PVS$pv
  sel_intervals = t(PVS$ci)
  
  saveRDS(list(active_vars=active_vars, sel_intervals=sel_intervals, pvalues=pvalues), file=outfile)
  
}


knockoff = function(method, outfile){
  if (method=="knockoff"){
    data=readRDS("real_data6.rds")
  } else if (method=="knockoff+"){
    data=readRDS("real_data7.rds")
  }
  
  X=data$X
  y=data$y
  
  offset=0
  if (method=="knockoff+"){
    offset=1
  }
  filter = knockoff.filter(X, y, fdr=q, offset=offset)
  saveRDS(list(active_vars = filter$rejected), file=outfile)
}


randomized = function(type, outfile){
  
  if (type=="full"){
    data=readRDS("real_data4.rds")
  } else if (type=="partial"){
    data=readRDS("real_data5.rds")
  }
  
  X=data$X
  y=data$y
  n=nrow(X)
  
  rand_lasso_soln = selectiveInference:::randomizedLasso(X, 
                                                         y, 
                                                         lambda*n, 
                                                         family=selectiveInference:::family_label(loss),
                                                         condition_subgrad=TRUE)
  
  full_targets=selectiveInference:::compute_target(rand_lasso_soln, type=type, sigma_est=sigma_est)
  
  PVS = selectiveInference:::randomizedLassoInf(rand_lasso_soln,
                                                targets=full_targets,
                                                sampler = "norejection", #"adaptMCMC", #
                                                level=0.9, 
                                                burnin=1000, 
                                                nsample=10000)
  active_vars=rand_lasso_soln$active_set
  cat("active_vars:",active_vars,"\n")
  pvalues = PVS$pvalues
  sel_intervals = t(PVS$ci)  # matrix with two rows
  saveRDS(list(active_vars=active_vars, sel_intervals=sel_intervals, pvalues=pvalues), file=outfile)
  
  return(NULL)
}

if (method=="liu"){
  liu_full(outfile)
} else if (method=="lee_full"){
  lee(outfile, type="full")
} else if (method=="lee_partial"){
  lee(outfile, type="partial")
} else if (method=="knockoff"){
  knockoff(method, outfile)
} else if (method=="knockoff+"){
  knockoff(method, outfile)
} else if (method=="randomized_full"){
  randomized("full", outfile)
} else if (method=="randomized_partial"){
  randomized("partial", outfile)
}





