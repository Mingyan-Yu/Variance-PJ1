library(R2jags)
library(dplyr)
library(ggplot2)
library(haven)
library(splines2)

## read in EMA data
load("combined_all.Rdata")

## read in survey data
survey <- read_dta("working covariates only for empatica linkage 061122.dta")

## extract subject ID in group 2 and group 3
comparison_group_ID <- unique(survey$idno[survey$group %in% c(2,3)])
comparison_group_ID <- comparison_group_ID[-which(comparison_group_ID==8707438592
)]

## subject IDs
subject_ID <- unique(combined_all$ID)
subject_ID <- intersect(subject_ID, comparison_group_ID)
combined_study <- combined_all[combined_all$ID %in% subject_ID,]

## cut-off time points: time tag1 to time tag6
## total observed HR of each subject
obs_time <- c()
for(i in 1:length(subject_ID)){
  data <- combined_study[combined_study$ID %in% subject_ID[i],]
  obs_time <- c(obs_time, length(data$HR[!is.na(data$HR)]))
}
## time when subjects reached time tag 6
tag6_time <- c()
for(i in 1:length(subject_ID)){
  data <- combined_study[combined_study$ID %in% subject_ID[i],]
  tag6_time <- c(tag6_time, which(data$tags=="tag 6"))
}
## which subject's HR ends before time tag 6
which(obs_time<tag6_time)
## four subjects whose observations ended before time tag 6
subject_ID <- subject_ID[-which(obs_time<tag6_time)]
combined_study <- combined_all[combined_all$ID %in% subject_ID,]

## cut off all observations after time tag 6
combined_study_cut <- list()
for(i in 1:length(subject_ID)){
  data <- combined_study[combined_study$ID %in% subject_ID[i],]
  data <- data[1:(which(data$tags=="tag 6")-1),]
  combined_study_cut <- append(combined_study_cut,
                               list(data))
}
combined_study_cut <- bind_rows(combined_study_cut)

## create indicator for each period in data
## five period: prepare, measure1, task, measure2, recovery
combined_study_cut_period <- list()
for(i in 1:length(subject_ID)){
  data <- combined_study_cut[combined_study_cut$ID %in% subject_ID[i],]
  rownames(data) <- 1:nrow(data)
  tags <- which(!is.na(data$tags))
  data$I_prepare <- rep(0, nrow(data))
  data$I_prepare[tags[1]:(tags[2]-1)] <- 1
  data$I_measure1 <- rep(0, nrow(data))
  data$I_measure1[tags[2]:(tags[3]-1)] <- 1
  data$I_task <- rep(0, nrow(data))
  data$I_task[tags[3]:(tags[5]-1)] <- 1
  data$I_measure2 <- rep(0, nrow(data))
  data$I_measure2[tags[5]:(tags[5]+300)] <- 1
  data$I_recovery <- 1-data$I_prepare-data$I_measure1-data$I_task-data$I_measure2
  combined_study_cut_period <- append(combined_study_cut_period,
                                      list(data))
}
combined_study_cut <- bind_rows(combined_study_cut_period)

## write up a function to create time knots per nth percentile
time_knot_pernquantile <- function(data, nq){
  ## only keep occasions where HR were not missing
  HR_end <- ifelse(sum(is.na(data$HR)) > 0,
                   min(data$occasion[which(is.na(data$HR))])-1,
                   length(data$HR))
  time_knots <- quantile(1:HR_end,
                         probs = seq(nq/100, 1-nq/100, by = nq/100))
  return(time_knots)
}

## write up a function to create B-spline basis matrix
## and period indicators for random effects
bSpline_mat_period <- function(data, nq){
  time_knots <- time_knot_pernquantile(data, nq)
  rownames(data) <- 1:nrow(data)
  tags2 <- data$occasion[data$tags %in% "tag 2"]
  tags3 <- data$occasion[data$tags %in% "tag 3"]
  tags5 <- data$occasion[data$tags %in% "tag 5"]
  tags5m <- tags5+301
  bSpline_mat <- bSpline(x = data$occasion, knots = time_knots,
                         degree = 3, intercept = TRUE)
  ## separate B-spline matrix before tag2 and after tag2
  tag2_sep <- c()
  for(i in 1:ncol(bSpline_mat)){
    if(bSpline_mat[tags2,i]==0){
      tag2_sep <- tag2_sep
    } else {
      tag2_sep <- c(tag2_sep, i)
    }
  }
  bSpline_mat_sep2 <- cbind(bSpline_mat[,1:tag2_sep[length(tag2_sep)]],
                            bSpline_mat[,tag2_sep[1]:tag2_sep[length(tag2_sep)]],
                            bSpline_mat[,(tag2_sep[length(tag2_sep)]+1):ncol(bSpline_mat)])
  colnames(bSpline_mat_sep2) <- 1:ncol(bSpline_mat_sep2)
  bSpline_mat_sep2[tags2:nrow(bSpline_mat_sep2),tag2_sep[1]:tag2_sep[length(tag2_sep)]] <- 0
  bSpline_mat_sep2[1:(tags2-1), (tag2_sep[length(tag2_sep)]+1):(tag2_sep[length(tag2_sep)]+length(tag2_sep))] <- 0
  
  ## separate B-spline matrix before tag3 and after tag3
  tag3_sep <- c()
  for(i in 1:ncol(bSpline_mat_sep2)){
    if(bSpline_mat_sep2[tags3,i]==0){
      tag3_sep <- tag3_sep
    } else {
      tag3_sep <- c(tag3_sep, i)
    }
  }
  bSpline_mat_sep3 <- cbind(bSpline_mat_sep2[,1:tag3_sep[length(tag3_sep)]],
                            bSpline_mat_sep2[,tag3_sep[1]:tag3_sep[length(tag3_sep)]],
                            bSpline_mat_sep2[,(tag3_sep[length(tag3_sep)]+1):ncol(bSpline_mat_sep2)])
  colnames(bSpline_mat_sep3) <- 1:ncol(bSpline_mat_sep3)
  bSpline_mat_sep3[tags3:nrow(bSpline_mat_sep3),tag3_sep[1]:tag3_sep[length(tag3_sep)]] <- 0
  bSpline_mat_sep3[1:(tags3-1), (tag3_sep[length(tag3_sep)]+1):(tag3_sep[length(tag3_sep)]+length(tag3_sep))] <- 0
  
  ## separate B-spline matrix before tag5 and after tag5
  tag5_sep <- c()
  for(i in 1:ncol(bSpline_mat_sep3)){
    if(bSpline_mat_sep3[tags5,i]==0){
      tag5_sep <- tag5_sep
    } else {
      tag5_sep <- c(tag5_sep, i)
    }
  }
  bSpline_mat_sep5 <- cbind(bSpline_mat_sep3[,1:tag5_sep[length(tag5_sep)]],
                            bSpline_mat_sep3[,tag5_sep[1]:tag5_sep[length(tag5_sep)]],
                            bSpline_mat_sep3[,(tag5_sep[length(tag5_sep)]+1):ncol(bSpline_mat_sep3)])
  colnames(bSpline_mat_sep5) <- 1:ncol(bSpline_mat_sep5)
  bSpline_mat_sep5[tags5:nrow(bSpline_mat_sep5),tag5_sep[1]:tag5_sep[length(tag5_sep)]] <- 0
  bSpline_mat_sep5[1:(tags5-1), (tag5_sep[length(tag5_sep)]+1):(tag5_sep[length(tag5_sep)]+length(tag5_sep))] <- 0
  
  ## separate B-spline matrix before tag5m and after tag5m
  tag5m_sep <- c()
  for(i in 1:ncol(bSpline_mat_sep5)){
    if(bSpline_mat_sep5[tags5m,i]==0){
      tag5m_sep <- tag5m_sep
    } else {
      tag5m_sep <- c(tag5m_sep, i)
    }
  }
  bSpline_mat_sep5m <- cbind(bSpline_mat_sep5[,1:tag5m_sep[length(tag5m_sep)]],
                             bSpline_mat_sep5[,tag5m_sep[1]:tag5m_sep[length(tag5m_sep)]],
                             bSpline_mat_sep5[,(tag5m_sep[length(tag5m_sep)]+1):ncol(bSpline_mat_sep5)])
  colnames(bSpline_mat_sep5m) <- 1:ncol(bSpline_mat_sep5m)
  bSpline_mat_sep5m[tags5m:nrow(bSpline_mat_sep5m),tag5m_sep[1]:tag5m_sep[length(tag5m_sep)]] <- 0
  bSpline_mat_sep5m[1:(tags5m-1), (tag5m_sep[length(tag5m_sep)]+1):(tag5m_sep[length(tag5m_sep)]+length(tag5m_sep))] <- 0
  
  ## create period random effects indicator
  knot_prepare <- rep(0, ncol(bSpline_mat_sep5m))
  knot_prepare[1:tag2_sep[length(tag2_sep)]] <- 1
  knot_measure1 <- rep(0, ncol(bSpline_mat_sep5m))
  knot_measure1[(tag2_sep[length(tag2_sep)]+1):tag3_sep[length(tag3_sep)]] <- 1
  knot_task <- rep(0, ncol(bSpline_mat_sep5m))
  knot_task[(tag3_sep[length(tag3_sep)]+1):tag5_sep[length(tag5_sep)]] <- 1
  knot_measure2 <- rep(0, ncol(bSpline_mat_sep5m))
  knot_measure2[(tag5_sep[length(tag5_sep)]+1):tag5m_sep[length(tag5m_sep)]] <- 1
  knot_recovery <- 1-knot_prepare-knot_measure1-knot_task-knot_measure2
  
  return(list(bSpline_mat_sep5m, knot_prepare, knot_measure1,
              knot_task, knot_measure2, knot_recovery))
}

## B-spline of degree 3: joint model
## prepare data for jointly modeling HR of all subjects
joint_data_period_outcome <- function(combined_data, n_subject, nq, survey){
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
  I_prepare <- matrix(Inf, nrow = N, ncol = n_max)  
  ## indicator for preparation period
  I_measure1 <- matrix(Inf, nrow = N, ncol = n_max)
  ## indicator for first measurement period
  I_task <- matrix(Inf, nrow = N, ncol = n_max)
  ## indicator for task period
  I_measure2 <- matrix(Inf, nrow = N, ncol = n_max)
  ## indicator for second measurement period
  I_recovery <- matrix(Inf, nrow = N, ncol = n_max)
  ## indicator for recovery period
  n_knots <- rep(0, N)  ## a vector of number of time knots
  for(i in 1:N){
    data <- combined_data[combined_data$ID %in% sample_ID[i],]
    data$occasion <- 1:nrow(data)
    HR[i,1:n[i]] <- data$HR[1:n[i]]
    I_prepare[i,1:n[i]] <- data$I_prepare[1:n[i]]
    I_measure1[i,1:n[i]] <- data$I_measure1[1:n[i]]
    I_task[i,1:n[i]] <- data$I_task[1:n[i]]
    I_measure2[i,1:n[i]] <- data$I_measure2[1:n[i]]
    I_recovery[i,1:n[i]] <- data$I_recovery[1:n[i]]
    n_knots[i] <- ncol(bSpline_mat_period(data, nq)[[1]])
  }
  knot_max <- max(n_knots)  ## max number of time knots
  Bspline <- array(Inf,dim = c(n_max, knot_max, N))
  knot_prepare <- matrix(Inf, nrow = N, ncol = knot_max)
  knot_measure1 <- matrix(Inf, nrow = N, ncol = knot_max)
  knot_task <- matrix(Inf, nrow = N, ncol = knot_max)
  knot_measure2 <- matrix(Inf, nrow = N, ncol = knot_max)
  knot_recovery <- matrix(Inf, nrow = N, ncol = knot_max)
  for(i in 1:N){
    data <- combined_data[combined_data$ID %in% sample_ID[i],]
    data$occasion <- 1:nrow(data)
    Bspline[1:n[i],1:n_knots[i],i] <- bSpline_mat_period(data, nq)[[1]][1:n[i],1:n_knots[i]]
    knot_prepare[i, 1:n_knots[i]] <- bSpline_mat_period(data, nq)[[2]]
    knot_measure1[i, 1:n_knots[i]] <- bSpline_mat_period(data, nq)[[3]]
    knot_task[i, 1:n_knots[i]] <- bSpline_mat_period(data, nq)[[4]]
    knot_measure2[i, 1:n_knots[i]] <- bSpline_mat_period(data, nq)[[5]]
    knot_recovery[i, 1:n_knots[i]] <- bSpline_mat_period(data, nq)[[6]]
  }
  
  ## insert group
  group <- c()
  for(i in 1:N){
    group <- c(group, unique(survey$group[survey$idno %in% sample_ID[i]]))
  }
  group <- group-2
  return(list(N=N, n=n, n_max=n_max, HR=HR, knot_max=knot_max,
              n_knots=n_knots, Bspline=Bspline,
              I_prepare=I_prepare, I_measure1=I_measure1,
              I_task=I_task, I_measure2=I_measure2, I_recovery=I_recovery,
              knot_prepare=knot_prepare, knot_measure1=knot_measure1,
              knot_task=knot_task, knot_measure2=knot_measure2,
              knot_recovery=knot_recovery,
              group=group))
}

## write up the joint model with periods
## different hyperpriors for different periods for now
joint.period.model.outcome <- function(){
  for(i in 1:N){
    for(j in 1:n[i]){
      HR[i,j] ~ dnorm(mu[i,j], 1/(sigma_wp[i]^2*I_prepare[i,j]+sigma_wm1[i]^2*I_measure1[i,j]+sigma_wt[i]^2*I_task[i,j]+sigma_wm2[i]^2*I_measure2[i,j]+
                                    sigma_wr[i]^2*I_recovery[i,j]))
      fit[i,j] <- beta+Bspline[j,1:n_knots[i],i]%*%b[i,1:n_knots[i]]
    }
    e0[i] ~ dnorm(0, (1-rho[i]^2)/sigma_wp[i]^2)
    mu[i,1] <- fit[i,1]+rho[i]*e0[i]
    e[i,1] <- HR[i,1]-fit[i,1]
    for(k in 2:n[i]){
      e[i,k] <- HR[i,k]-fit[i,k]
      mu[i,k] <- fit[i,k]+rho[i]*e[i,k-1]
    }
  }
  
  beta ~ dnorm(0,1e-6)
  for(i in 1:N){
    for(j in 1:n_knots[i]){
      b[i,j] ~ dnorm(0,1/(sigma_bp[i]^2*knot_prepare[i,j]+sigma_bm1[i]^2*knot_measure1[i,j]+sigma_bt[i]^2*knot_task[i,j]+sigma_bm2[i]^2*knot_measure2[i,j]+
                            sigma_br[i]^2*knot_recovery[i,j]))
    }
    for(k in (n_knots[i]+1):knot_max){
      b[i,k] <- 0
    }
  }
  for(i in 1:N){
    var_wp[i] ~ dlnorm(sig_wp,1/psi_wp^2)
    var_wm1[i] ~ dlnorm(sig_wm1,1/psi_wm1^2)
    var_wt[i] ~ dlnorm(sig_wt,1/psi_wt^2)
    var_wm2[i] ~ dlnorm(sig_wm2,1/psi_wm2^2)
    var_wr[i] ~ dlnorm(sig_wr,1/psi_wr^2)
    sigma_wp[i] <- sqrt(var_wp[i])
    sigma_wm1[i] <- sqrt(var_wm1[i])
    sigma_wt[i] <- sqrt(var_wt[i])
    sigma_wm2[i] <- sqrt(var_wm2[i])
    sigma_wr[i] <- sqrt(var_wr[i])
    var_bp[i] ~ dlnorm(sig_bp,1/psi_bp^2)
    var_bm1[i] ~ dlnorm(sig_bm1,1/psi_bm1^2)
    var_bt[i] ~ dlnorm(sig_bt,1/psi_bt^2)
    var_bm2[i] ~ dlnorm(sig_bm2,1/psi_bm2^2)
    var_br[i] ~ dlnorm(sig_br,1/psi_br^2)
    sigma_bp[i] <- sqrt(var_bp[i])
    sigma_bm1[i] <- sqrt(var_bm1[i])
    sigma_bt[i] <- sqrt(var_bt[i])
    sigma_bm2[i] <- sqrt(var_bm2[i])
    sigma_br[i] <- sqrt(var_br[i])
    rho[i] ~ dunif(0,1)
  }
  sig_wp ~ dnorm(0, 1e-6)
  psi_wp ~ dt(0, 1/2.5^2,1);T(0,)
  sig_bp ~ dnorm(0, 1e-6)
  psi_bp ~ dt(0, 1/2.5^2,1);T(0,)
  sig_wm1 ~ dnorm(0, 1e-6)
  psi_wm1 ~ dt(0, 1/2.5^2,1);T(0,)
  sig_bm1 ~ dnorm(0, 1e-6)
  psi_bm1 ~ dt(0, 1/2.5^2,1);T(0,)
  sig_wt ~ dnorm(0, 1e-6)
  psi_wt ~ dt(0, 1/2.5^2,1);T(0,)
  sig_bt ~ dnorm(0, 1e-6)
  psi_bt ~ dt(0, 1/2.5^2,1);T(0,)
  sig_wm2 ~ dnorm(0, 1e-6)
  psi_wm2 ~ dt(0, 1/2.5^2,1);T(0,)
  sig_bm2 ~ dnorm(0, 1e-6)
  psi_bm2 ~ dt(0, 1/2.5^2,1);T(0,)
  sig_wr ~ dnorm(0, 1e-6)
  psi_wr ~ dt(0, 1/2.5^2,1);T(0,)
  sig_br ~ dnorm(0, 1e-6)
  psi_br ~ dt(0, 1/2.5^2,1);T(0,)
  
  ## outcome model
  for(i in 1:N){
    p[i] <- phi(alpha[1]+alpha[2]*log(sigma_wp[i]*10)+alpha[3]*log(sigma_wm1[i]*10)+
                  alpha[4]*log(sigma_wt[i]*10)+alpha[5]*log(sigma_wm2[i]*10)+
                  alpha[6]*log(sigma_wr[i]*10)+
                  alpha[7]*log(sigma_bp[i]/10)+alpha[8]*log(sigma_bm1[i]/10)+
                  alpha[9]*log(sigma_bt[i]/10)+alpha[10]*log(sigma_bm2[i]/10)+
                  alpha[11]*log(sigma_br[i]/10))
    group[i] ~ dbern(p[i])
  }
  for(i in 1:11){
    alpha[i] ~ dnorm(0, 1e-2)
  }
}

## run the model on a few subjects
set.seed(500)
HR.data <- joint_data_period_outcome(combined_data = combined_study_cut,
                                     n_subject = 104, nq = 2, survey = survey)
table(HR.data$group)
HR.params <- c("sigma_wp", "sigma_wm1", "sigma_wt", "sigma_wm2", "sigma_wr",
               "sigma_bp", "sigma_bm1", "sigma_bt", "sigma_bm2", "sigma_br",
               "sig_wp", "psi_wp", "sig_bp", "psi_bp",
               "sig_wm1", "psi_wm1", "sig_bm1", "psi_bm1",
               "sig_wt", "psi_wt", "sig_bt", "psi_bt",
               "sig_wm2", "psi_wm2", "sig_bm2", "psi_bm2",
               "sig_wr", "psi_wr", "sig_br", "psi_br",
               "beta", "b", "rho", "alpha")
HR.out <- jags(data = HR.data, 
               parameters.to.save = HR.params,
               model.file = joint.period.model.outcome, 
               n.burnin = 4000, n.iter = 8000, n.thin = 1,
               n.chains = 3)

print(HR.out)
traceplot(HR.out, varname=c("sigma_wp", "sigma_wm1", "sigma_wt", "sigma_wm2", "sigma_wr",
                            "sigma_bp", "sigma_bm1", "sigma_bt", "sigma_bm2", "sigma_br",
                            "sig_wp", "psi_wp", "sig_bp", "psi_bp",
                            "sig_wm1", "psi_wm1", "sig_bm1", "psi_bm1",
                            "sig_wt", "psi_wt", "sig_bt", "psi_bt",
                            "sig_wm2", "psi_wm2", "sig_bm2", "psi_bm2",
                            "sig_wr", "psi_wr", "sig_br", "psi_br",
                            "beta", "rho", "alpha"))
attach.jags(HR.out)
save(alpha, file = "Joint_outcome_model_Bspline_period_AR1_alpha.Rdata")

save(sigma_wp, file = "Joint_outcome_model_Bspline_period_AR1_sigma_wp.Rdata")
save(sigma_wm1, file = "Joint_outcome_model_Bspline_period_AR1_sigma_wm1.Rdata")
save(sigma_wt, file = "Joint_outcome_model_Bspline_period_AR1_sigma_wt.Rdata")
save(sigma_wm2, file = "Joint_outcome_model_Bspline_period_AR1_sigma_wm2.Rdata")
save(sigma_wr, file = "Joint_outcome_model_Bspline_period_AR1_sigma_wr.Rdata")

save(sigma_bp, file = "Joint_outcome_model_Bspline_period_AR1_sigma_bp.Rdata")
save(sigma_bm1, file = "Joint_outcome_model_Bspline_period_AR1_sigma_bm1.Rdata")
save(sigma_bt, file = "Joint_outcome_model_Bspline_period_AR1_sigma_bt.Rdata")
save(sigma_bm2, file = "Joint_outcome_model_Bspline_period_AR1_sigma_bm2.Rdata")
save(sigma_br, file = "Joint_outcome_model_Bspline_period_AR1_sigma_br.Rdata")

save(sig_wp, file = "Joint_outcome_model_Bspline_period_AR1_sig_wp.Rdata")
save(sig_wm1, file = "Joint_outcome_model_Bspline_period_AR1_sig_wm1.Rdata")
save(sig_wt, file = "Joint_outcome_model_Bspline_period_AR1_sig_wt.Rdata")
save(sig_wm2, file = "Joint_outcome_model_Bspline_period_AR1_sig_wm2.Rdata")
save(sig_wr, file = "Joint_outcome_model_Bspline_period_AR1_sig_wr.Rdata")

save(sig_bp, file = "Joint_outcome_model_Bspline_period_AR1_sig_bp.Rdata")
save(sig_bm1, file = "Joint_outcome_model_Bspline_period_AR1_sig_bm1.Rdata")
save(sig_bt, file = "Joint_outcome_model_Bspline_period_AR1_sig_bt.Rdata")
save(sig_bm2, file = "Joint_outcome_model_Bspline_period_AR1_sig_bm2.Rdata")
save(sig_br, file = "Joint_outcome_model_Bspline_period_AR1_sig_br.Rdata")

save(psi_wp, file = "Joint_outcome_model_Bspline_period_AR1_psi_wp.Rdata")
save(psi_wm1, file = "Joint_outcome_model_Bspline_period_AR1_psi_wm1.Rdata")
save(psi_wt, file = "Joint_outcome_model_Bspline_period_AR1_psi_wt.Rdata")
save(psi_wm2, file = "Joint_outcome_model_Bspline_period_AR1_psi_wm2.Rdata")
save(psi_wr, file = "Joint_outcome_model_Bspline_period_AR1_psi_wr.Rdata")

save(psi_bp, file = "Joint_outcome_model_Bspline_period_AR1_psi_bp.Rdata")
save(psi_bm1, file = "Joint_outcome_model_Bspline_period_AR1_psi_bm1.Rdata")
save(psi_bt, file = "Joint_outcome_model_Bspline_period_AR1_psi_bt.Rdata")
save(psi_bm2, file = "Joint_outcome_model_Bspline_period_AR1_psi_bm2.Rdata")
save(psi_br, file = "Joint_outcome_model_Bspline_period_AR1_psi_br.Rdata")

save(beta, file = "fixed_effect.Rdata")
save(b, file = "random_effects.Rdata")
save(rho, file = "autocorrelation.Rdata")













