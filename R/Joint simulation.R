library(dplyr)
library(R2jags)
library(splines2)
## write out true values from true model
beta_true <- 81.083
sig_e_true <- 1.8
psi_e_true <- 0.75
sig_b_true <- 5.0
psi_b_true <- 1
alpha_true <- c(-0.6, 1, -1.5)
rho_true <- 0.998

## write up a function to create time knots per nth quantile
time_knot_sim_mat <- function(occasion, nq){
  time_knots <- quantile(occasion,
                         probs = seq(nq/100, 1-nq/100, by = nq/100))
  return(time_knots)
}

## write up a function to simulate HR data
simulation <- function(beta, sig_e, psi_e, sig_b, psi_b,
                       n_obs, nq, n_subject, alpha, rho){
  simulated_data <- list()
  occasion <- 1:n_obs
  knots <- time_knot_sim_mat(occasion, nq = nq)
  bspline_mat <- bSpline(x = occasion, knots = knots,
                         degree = 3, intercept = TRUE)
  ## generate variance of random effects
  var_random <- rlnorm(n = n_subject, meanlog = sig_b,
                       sdlog = psi_b)
  sd_random <- sqrt(var_random)
  
  ## generate variance of random errors
  var_error <- rlnorm(n = n_subject, meanlog = sig_e,
                      sdlog = psi_e)
  sd_error <- sqrt(var_error)
  
  for(i in 1:n_subject){
    ## generate random effects
    b <- rnorm(n = ncol(bspline_mat), mean = 0, sd = sd_random[i])
    
    ## simulate expected HR
    HR_sim <- beta+bspline_mat%*%b
    
    ## generate random errors
    w <- rnorm(n = n_obs, mean = 0, sd = sd_error[i])
    e <- stats::filter(w, filter = c(rho), method = "recursive", init = 0)
    
    ## add in errors to expected HR
    HR_sim <- HR_sim+e
    
    data <- data.frame("ID" <- i, "HR" <- HR_sim)
    simulated_data <- append(simulated_data, list(data))
  }
  simulated_data <- bind_rows(simulated_data)
  colnames(simulated_data) <- c("ID", "HR")
  
  ## simulate outcome
  p <- pnorm(alpha[1]+alpha[2]*log(sd_error)+alpha[3]*log(sd_random/10))
  group <- rbinom(n = n_subject, size = 1, prob = p)
  
  return(list(simulated_data, sd_error, sd_random, rho, p, group))
}

## write up a function to create time knots per nth quantile
time_knot_pernquantile <- function(data, nq){
  ## only keep occasions where HR were not missing
  HR_end <- ifelse(sum(is.na(data$HR)) > 0,
                   min(data$occasion[which(is.na(data$HR))])-1,
                   length(data$HR))
  time_knots <- quantile(1:HR_end,
                         probs = seq(nq/100, 1-nq/100, by = nq/100))
  return(time_knots)
}

## B-spline of degree 3: joint model
## prepare data for jointly modeling HR of all subjects
joint_data <- function(combined_data, n_subject, nq, group){
  subject_ID <- unique(combined_data$ID)
  sample_ID <- sample(subject_ID, size = n_subject, replace = FALSE)
  print(sample_ID)
  N <- n_subject  ## number of subjects
  n <- rep(0, N)  ## a vector of number of observations
  for(i in 1:N){
    data <- combined_data[combined_data$ID %in% sample_ID[i],]
    n[i] <- length(data$HR[!is.na(data$HR)])
  }
  n_max <- max(n)  ## max number of observations
  HR <- matrix(Inf,nrow = N, ncol = n_max)  ## a matrix of HR
  ## each row represents one subject's HR data
  n_knots <- rep(0, N)  ## a vector of number of time knots
  for(i in 1:N){
    data <- combined_data[combined_data$ID %in% sample_ID[i],]
    data$occasion <- 1:nrow(data)
    HR[i,1:n[i]] <- data$HR[1:n[i]]
    time_knots <- time_knot_pernquantile(data, nq)
    n_knots[i] <- length(time_knots)+1+3
  }
  knot_max <- max(n_knots)  ## max number of time knots
  Bspline <- array(Inf,dim = c(n_max, knot_max, N))
  for(i in 1:N){
    data <- combined_data[combined_data$ID %in% sample_ID[i],]
    data$occasion <- 1:nrow(data)
    time_knots <- time_knot_pernquantile(data, nq)
    Bspline[1:n[i],1:n_knots[i],i] <- bSpline(x = data$occasion,
                                              knots = time_knots, degree = 3,
                                              intercept = TRUE)
  }
  group <- group[sample_ID]
  return(list(N=N, n=n, n_max=n_max, HR=HR, knot_max=knot_max,
              n_knots=n_knots, Bspline=Bspline, group=group))
}

## write up the joint model
joint.model <- function(){
  for(i in 1:N){
    for(j in 1:n[i]){
      HR[i,j] ~ dnorm(mu[i,j], 1/sigma_e[i]^2)
      fit[i,j] <- beta+Bspline[j,1:n_knots[i],i]%*%b[i,1:n_knots[i]]
    }
    e0[i] ~ dnorm(0, (1-rho^2)/sigma_e[i]^2)
    mu[i,1] <- fit[i,1]+rho*e0[i]
    e[i,1] <- HR[i,1]-fit[i,1]
    for(k in 2:n[i]){
      e[i,k] <- HR[i,k]-fit[i,k]
      mu[i,k] <- fit[i,k]+rho*e[i,k-1]
    }
  }
  beta ~ dnorm(0,1e-6)
  for(i in 1:N){
    for(j in 1:n_knots[i]){
      b[i,j] ~ dnorm(0,1/sigma_b[i]^2)
    }
    for(k in (n_knots[i]+1):knot_max){
      b[i,k] <- 0
    }
  }
  for(i in 1:N){
    var_e[i] ~ dlnorm(sig_e,1/psi_e^2)
    var_b[i] ~ dlnorm(sig_b,1/psi_b^2)
    sigma_e[i] <- sqrt(var_e[i])
    sigma_b[i] <- sqrt(var_b[i])
  }
  rho ~ dunif(0,1)
  sig_e ~ dnorm(0, 1e-6)
  psi_e ~ dt(0, 1/2.5^2,1);T(0,)
  sig_b ~ dnorm(0, 1e-6)
  psi_b ~ dt(0, 1/2.5^2,1);T(0,)
  
  ## outcome model
  for(i in 1:N){
    p[i] <- pnorm(alpha[1]+alpha[2]*log(sigma_e[i])+alpha[3]*log(sigma_b[i]/10), 0, 1)
    group[i] ~ dbern(p[i])
  }
  for(i in 1:3){
    alpha[i] ~ dnorm(0, 1e-2)
  }
}


## start simulation and save 95% CI of posteriors
beta_simulate <- list()
sig_e_simulate <- list()
psi_e_simulate <- list()
sig_b_simulate <- list()
psi_b_simulate <- list()
alpha_simulate <- list()
rho_simulate <- list()

beta_estimate <- c()
sig_e_estimate <- c()
psi_e_estimate <- c()
sig_b_estimate <- c()
psi_b_estimate <- c()
rho_estimate <- c()
alpha_estimate <- list()

k <- Sys.getenv("SLURM_ARRAY_TASK_ID")

## simulate a few subjects
sim_bspline_data <- simulation(beta = beta_true, sig_e = sig_e_true,
                               psi_e = psi_e_true, sig_b = sig_b_true,
                               psi_b = psi_b_true, n_obs = 600,
                               nq = 2, n_subject = 150, alpha = alpha_true,
                               rho = rho_true)
table(sim_bspline_data[[6]])
hist(sim_bspline_data[[5]], xlab = "Probability inside the Bernoulli distribution", main = "Probability Distribution")

sim_bspline_data.model <- joint_data(combined_data = sim_bspline_data[[1]],
                                     n_subject = 150, nq = 2,
                                     group = sim_bspline_data[[6]])
model.params <- c("sigma_e", "sigma_b", "sig_e", "psi_e", "sig_b", "psi_b",
                  "beta", "b", "rho", "alpha")
model.out <- jags.parallel(data = sim_bspline_data.model, 
                           parameters.to.save = model.params,
                           model.file = joint.model, n.burnin = 4000,
                           n.iter = 8000, n.thin = 1,
                           n.chains = 3)

print(model.out)
traceplot(model.out, varname=c("beta","sigma_b","sigma_e","sig_b","psi_b",
                               "sig_e","psi_e","rho","alpha"))

attach.jags(model.out)
quantile(beta, probs = c(0.025, 0.975))
quantile(sig_b, probs = c(0.025, 0.975))
quantile(psi_b, probs = c(0.025, 0.975))
quantile(sig_e, probs = c(0.025, 0.975))
quantile(psi_e, probs = c(0.025, 0.975))
apply(sigma_b, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
apply(sigma_e, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
apply(alpha, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
mean(beta)
mean(sig_b)
mean(psi_b)
mean(sig_e)
mean(psi_e)
quantile(rho, probs = c(0.025, 0.975))
mean(rho)


## save the results
beta_simulate <- append(beta_simulate, 
                        list(quantile(beta, probs = c(0.025,0.975))))
sig_e_simulate <- append(sig_e_simulate, 
                         list(quantile(sig_e, probs = c(0.025,0.975))))
psi_e_simulate <- append(psi_e_simulate, 
                         list(quantile(psi_e, probs = c(0.025,0.975))))
sig_b_simulate <- append(sig_b_simulate, 
                         list(quantile(sig_b, probs = c(0.025,0.975))))
psi_b_simulate <- append(psi_b_simulate, 
                         list(quantile(psi_b, probs = c(0.025,0.975))))
rho_simulate <- append(rho_simulate, 
                       list(quantile(rho, probs = c(0.025,0.975))))                      
alpha_simulate <- append(alpha_simulate,
                         list(apply(alpha, MARGIN = 2, FUN = quantile, probs = c(0.025,0.975))))
  
beta_estimate <- c(beta_estimate, mean(beta))
sig_e_estimate <- c(sig_e_estimate, mean(sig_e))
psi_e_estimate <- c(psi_e_estimate, mean(psi_e))
sig_b_estimate <- c(sig_b_estimate, mean(sig_b))
psi_b_estimate <- c(psi_b_estimate, mean(psi_b))
rho_estimate <- c(rho_estimate, mean(rho))
alpha_estimate <- append(alpha_estimate, list(colMeans(alpha)))


save(beta_simulate, file = paste0(k, "beta_simulate.Rdata"))
save(sig_e_simulate, file = paste0(k, "sig_e_simulate.Rdata"))
save(psi_e_simulate, file = paste0(k, "psi_e_simulate.Rdata"))
save(sig_b_simulate, file = paste0(k, "sig_b_simulate.Rdata"))
save(psi_b_simulate, file = paste0(k, "psi_b_simulate.Rdata"))
save(rho_simulate, file = paste0(k, "rho_simulate.Rdata"))
save(alpha_simulate, file = paste0(k, "alpha_simulate.Rdata"))

save(beta_estimate, file = paste0(k, "beta_estimate.Rdata"))
save(rho_estimate, file = paste0(k, "rho_estimate.Rdata"))
save(sig_e_estimate, file = paste0(k, "sig_e_estimate.Rdata"))
save(psi_e_estimate, file = paste0(k, "psi_e_estimate.Rdata"))
save(sig_b_estimate, file = paste0(k, "sig_b_estimate.Rdata"))
save(psi_b_estimate, file = paste0(k, "psi_b_estimate.Rdata"))
save(alpha_estimate, file = paste0(k, "alpha_estimate.Rdata"))

