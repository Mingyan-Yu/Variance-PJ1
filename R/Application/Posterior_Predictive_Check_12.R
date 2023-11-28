library(R2jags)
library(dplyr)
library(ggplot2)
library(haven)
library(splines2)

## read in EMA data
load("combined_all.Rdata")

## read in survey data
survey <- read_dta("working covariates only for empatica linkage 061122.dta")

## extract subject ID in group 1 and group 2
comparison_group_ID <- unique(survey$idno[survey$group %in% c(1,2)])

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
## two subjects whose observations ended before time tag 6
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
  group <- group-1
  return(list(N=N, n=n, n_max=n_max, HR=HR, knot_max=knot_max,
              n_knots=n_knots, Bspline=Bspline,
              I_prepare=I_prepare, I_measure1=I_measure1,
              I_task=I_task, I_measure2=I_measure2, I_recovery=I_recovery,
              knot_prepare=knot_prepare, knot_measure1=knot_measure1,
              knot_task=knot_task, knot_measure2=knot_measure2,
              knot_recovery=knot_recovery,
              group=group))
}

## run the model on a few subjects
set.seed(500)
HR.data <- joint_data_period_outcome(combined_data = combined_study_cut,
                                     n_subject = 95, nq = 2, survey = survey)
table(HR.data$group)

## Posterior predictive check for the predictors
## load application results
load("autocorrelation.Rdata")
load("fixed_effect.Rdata")
load("random_effects.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_sigma_wm1.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_sigma_wp.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_sigma_wm2.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_sigma_wt.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_sigma_wr.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_alpha.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_sigma_bm1.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_sigma_bp.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_sigma_bm2.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_sigma_bt.Rdata")
load("Joint_outcome_model_Bspline_period_AR1_sigma_br.Rdata")

## calculate fitted values
set.seed(500)
random_effects <- b[sample(1:dim(b)[1], size = 500),,]
fixed_effects <- beta[sample(1:dim(beta)[1], size = 500),]
AR1 <- rho[sample(1:dim(rho)[1], size = 500),]
sigma_wp <- sigma_wp[sample(1:dim(sigma_wp)[1], size = 500), ]
alpha <- alpha[sample(1:dim(alpha)[1], size = 500), ]

Sys.time()

fitted_values_12 <- list()
f_values_12 <- list()
for(i in 1:HR.data$N){
  fit_i <- matrix(data = NA, nrow = 500, ncol = HR.data$n_max)
  mu_i <- matrix(data = NA, nrow = 500, ncol = HR.data$n_max)
  for(j in 1:500){
    for(k in 1:HR.data$n[i]){
      fit_i[j,k] <- fixed_effects[j] + HR.data$Bspline[k,1:HR.data$n_knots[i],i]%*%random_effects[j,i,1:HR.data$n_knots[i]]
    }
    e_i0 <- rnorm(n = 1, mean = 0, 
                  sd = sqrt(sigma_wp[j,i]^2/(1-AR1[j,i]^2)))
    w <- rnorm(n = HR.data$n[i], mean = 0,
               sd = sigma_wp[j,i]*HR.data$I_prepare[i,1:HR.data$n[i]]+
                 sigma_wm1[j,i]*HR.data$I_measure1[i,1:HR.data$n[i]]+
                 sigma_wt[j,i]*HR.data$I_task[i,1:HR.data$n[i]]+
                 sigma_wm2[j,i]*HR.data$I_measure2[i,1:HR.data$n[i]]+
                 sigma_wr[j,i]*HR.data$I_recovery[i,1:HR.data$n[i]])
    e_i <- rep(0, HR.data$n[i])
    e_i[1] <- AR1[j,i]*e_i0 + w[1]
    mu_i[j,1] <- fit_i[j,1] + AR1[j,i]*e_i0
    for(s in 2:HR.data$n[i]){
      e_i[s] <- AR1[j,i]*e_i[s-1] + w[s]
      mu_i[j,s] <- fit_i[j,s] + AR1[j,i]*e_i[s-1]
    }
  }
  fitted_values_12 <- append(fitted_values_12, list(mu_i))
  f_values_12 <- append(f_values_12, list(fit_i))
}

Sys.time()

## calculate replicated values
sigma_wm1 <- sigma_wm1[sample(1:dim(sigma_wm1)[1], size = 500), ]
sigma_wt <- sigma_wt[sample(1:dim(sigma_wt)[1], size = 500), ]
sigma_wm2 <- sigma_wm2[sample(1:dim(sigma_wm2)[1], size = 500), ]
sigma_wr <- sigma_wr[sample(1:dim(sigma_wr)[1], size = 500), ]

sigma_bp <- sigma_bp[sample(1:dim(sigma_bp)[1], size = 500), ]
sigma_bm1 <- sigma_bm1[sample(1:dim(sigma_bm1)[1], size = 500), ]
sigma_bt <- sigma_bt[sample(1:dim(sigma_bt)[1], size = 500), ]
sigma_bm2 <- sigma_bm2[sample(1:dim(sigma_bm2)[1], size = 500), ]
sigma_br <- sigma_br[sample(1:dim(sigma_br)[1], size = 500), ]

Sys.time()

rep_values_12 <- list()
for(i in 1:HR.data$N){
  rep_i <- matrix(data = NA, nrow = 500, ncol = HR.data$n_max)
  for(j in 1:500){
    e_i0 <- rnorm(n = 1, mean = 0, 
                  sd = sqrt(sigma_wp[j,i]^2/(1-AR1[j,i]^2)))
    w <- rnorm(n = HR.data$n[i], mean = 0,
               sd = sigma_wp[j,i]*HR.data$I_prepare[i,1:HR.data$n[i]]+
                 sigma_wm1[j,i]*HR.data$I_measure1[i,1:HR.data$n[i]]+
                 sigma_wt[j,i]*HR.data$I_task[i,1:HR.data$n[i]]+
                 sigma_wm2[j,i]*HR.data$I_measure2[i,1:HR.data$n[i]]+
                 sigma_wr[j,i]*HR.data$I_recovery[i,1:HR.data$n[i]])
    e_i <- rep(0, HR.data$n[i])
    e_i[1] <- AR1[j,i]*e_i0 + w[1]
    rep_i[j,1] <- f_values_12[[i]][j,1] + e_i[1]
    for(s in 2:HR.data$n[i]){
      e_i[s] <- AR1[j,i]*e_i[s-1] + w[s]
      rep_i[j,s] <- f_values_12[[i]][j,s] + e_i[s]
    }
  }
  rep_values_12 <- append(rep_values_12, list(rep_i))
}

Sys.time()

## Compute p-values
T_obs_12 <- list()
for(i in 1:HR.data$N){
  T_obs_i <- rep(0, 500)
  for(j in 1:500){
    T_obs_i[j] <- sum((HR.data$HR[i,1:HR.data$n[i]]-fitted_values_12[[i]][j,1:HR.data$n[i]])^2/
                        (sigma_wp[j,i]^2*HR.data$I_prepare[i,1:HR.data$n[i]]+
                           sigma_wm1[j,i]^2*HR.data$I_measure1[i,1:HR.data$n[i]]+
                           sigma_wt[j,i]^2*HR.data$I_task[i,1:HR.data$n[i]]+
                           sigma_wm2[j,i]^2*HR.data$I_measure2[i,1:HR.data$n[i]]+
                           sigma_wr[j,i]^2*HR.data$I_recovery[i,1:HR.data$n[i]]))
  }
  T_obs_12 <- append(T_obs_12, list(T_obs_i))
}

T_rep_12 <- list()
for(i in 1:HR.data$N){
  T_rep_i <- rep(0, 500)
  for(j in 1:500){
    T_rep_i[j] <- sum((rep_values_12[[i]][j,1:HR.data$n[i]]-fitted_values_12[[i]][j,1:HR.data$n[i]])^2/
                        (sigma_wp[j,i]^2*HR.data$I_prepare[i,1:HR.data$n[i]]+
                           sigma_wm1[j,i]^2*HR.data$I_measure1[i,1:HR.data$n[i]]+
                           sigma_wt[j,i]^2*HR.data$I_task[i,1:HR.data$n[i]]+
                           sigma_wm2[j,i]^2*HR.data$I_measure2[i,1:HR.data$n[i]]+
                           sigma_wr[j,i]^2*HR.data$I_recovery[i,1:HR.data$n[i]]))
  }
  T_rep_12 <- append(T_rep_12, list(T_rep_i))
}

PPD12 <- rep(0, HR.data$N)
for(i in 1:HR.data$N){
  PPD12[i] <- mean(T_obs_12[[i]]<T_rep_12[[i]])
}
hist(PPD12)
summary(PPD12)

save(PPD12, file = "PPD12.Rdata")
save(fitted_values_12, file = "fitted_values_12.Rdata")
save(f_values_12, file = "f_values_12.Rdata")
save(rep_values_12, file = "rep_values_12")
save(T_rep_12, file = "T_rep_12.Rdata")
save(T_obs_12, file = "T_obs_12.Rdata")


## posterior predictive check for the outcome
fitted_outcome <- list()
T_rep_12_outcome <- rep(0, 500)
for(j in 1:500){
  fit_j <- rep(0, HR.data$N)
  for(i in 1:HR.data$N){
    eta <- alpha[j,1]+alpha[j,2]*log(sigma_wp[j,i]*10)+
      alpha[j,3]*log(sigma_wm1[j,i]*10)+
      alpha[j,4]*log(sigma_wt[j,i]*10)+
      alpha[j,5]*log(sigma_wm2[j,i]*10)+
      alpha[j,6]*log(sigma_wr[j,i]*10)+
      alpha[j,7]*log(sigma_bp[j,i]/10)+
      alpha[j,8]*log(sigma_bm1[j,i]/10)+
      alpha[j,9]*log(sigma_bt[j,i]/10)+
      alpha[j,10]*log(sigma_bm2[j,i]/10)+
      alpha[j,11]*log(sigma_br[j,i]/10)
    prob <- pnorm(q = eta, mean = 0, sd = 1)
    fit_j[i] <- rbinom(n = 1, size = 1, prob = prob)
  }
  fitted_outcome <- append(fitted_outcome, list(fit_j))
  T_rep_12_outcome[j] <- sum(fit_j)
}

mean(sum(HR.data$group) < T_rep_12_outcome)

## check the tail of the posterior predictive distribution
set.seed(500)
Tk_rep_12_outcome <- matrix(data = NA, nrow = 10, ncol = 500)

Tk_obs_12_outcome <- matrix(data = NA, nrow = 10, ncol = 500)

for(j in 1:500){
  eta_j <- rep(0, HR.data$N)
  fit_j <- rep(0, HR.data$N)
  for(i in 1:HR.data$N){
    eta_j[i] <- alpha[j,1]+alpha[j,2]*log(sigma_wp[j,i]*10)+
      alpha[j,3]*log(sigma_wm1[j,i]*10)+
      alpha[j,4]*log(sigma_wt[j,i]*10)+
      alpha[j,5]*log(sigma_wm2[j,i]*10)+
      alpha[j,6]*log(sigma_wr[j,i]*10)+
      alpha[j,7]*log(sigma_bp[j,i]/10)+
      alpha[j,8]*log(sigma_bm1[j,i]/10)+
      alpha[j,9]*log(sigma_bt[j,i]/10)+
      alpha[j,10]*log(sigma_bm2[j,i]/10)+
      alpha[j,11]*log(sigma_br[j,i]/10)
    prob <- pnorm(q = eta_j[i], mean = 0, sd = 1)
    fit_j[i] <- rbinom(n = 1, size = 1, prob = prob)
  }
  for(k in 1:10){
    decile <- ntile(eta_j, n = 10)
    Tk_obs_12_outcome[k,j] <- sum(HR.data$group[which(decile==k)])
    Tk_rep_12_outcome[k,j] <- sum(fit_j[which(decile==k)])
  }
}

for(k in 1:10){
  print(mean(Tk_obs_12_outcome[k,] < Tk_rep_12_outcome[k,]))
}

save(fitted_outcome, file = "fitted_outcome.Rdata")
save(T_rep_12_outcome, file = "T_rep_12_outcome.Rdata")
save(Tk_obs_12_outcome, file = "Tk_obs_12_outcome.Rdata")
save(Tk_rep_12_outcome, file = "Tk_rep_12_outcome.Rdata")






