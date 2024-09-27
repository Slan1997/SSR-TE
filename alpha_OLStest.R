library(pacman)
p_load(mvtnorm, # generate multivariate normal
       #PwrGSD, # Lan-DeMets bounds with both efficacy and futility
       gsDesign, # Lan-DeMets bounds with efficacy only
       tidyr,dplyr,readr,magrittr
)
### this file was run 72 times for 72 scenarios.

# if use cluster to run
# Sys.getenv("SLURM_ARRAY_TASK_ID") extracts array id from slurm file: "#SBATCH --array=1-72"
index=Sys.getenv("SLURM_ARRAY_TASK_ID")
scena = as.numeric(index)   # scena = 1:72

nsim = 1e6 # 1 million
sub_sim_chunk = 1:nsim 

### initial scenarios investigated
scenarios = expand.grid(n_control = c(16,30,60,120,1000,2000),
                        alloc_ratio = c(1,2),
                        SSR = c('no', 'yes'),
                        beta = c(.1,.2)
)  %>%
  as.data.frame
### added more scenarios of large sample sizes to see the trend better
scenarios1 = expand.grid(n_control = c(1500,2500,3000),
                         alloc_ratio = c(1,2),
                         SSR = c('no', 'yes'),
                         beta = c(.1,.2)
)  %>%
  as.data.frame

scenarios = bind_rows(scenarios,scenarios1) 
scenarios
scenarios[scena,]

# number of endpoints
n_endpts = 6
# target alpha for Lan-DeMets boundaries
alpha0 = .025
# target beta
beta0 =scenarios[scena,'beta']
beta0
# two groups, treatment vs control. allocation ratio r = n_trt/n_ct
alloc_rate = scenarios[scena,'alloc_ratio']
alloc_rate
# control group sample size for two stages
N_ct = scenarios[scena,'n_control']
# control group sample size for each stage 
n_ct = N_ct/2 
# treatment group sample size for two stages
N_trt = N_ct*alloc_rate
# treatment group sample size for each stage 
n_trt = N_trt/2  

# prefixed maximal SSR sample size for control group, 2-fold of original planned sample size
N_ct_max = 2*N_ct

# SSR scenarios
SSR = scenarios[scena,'SSR']
SSR

##### Get critical values from O'Brien-Fleming spending function
alpha <- alpha0 # 1-sided Type I error
k <- 2 # Number of planned analyses
test.type <- 1 # Asymmetric bound design with non-binding futility bound
timing <- .5 # information fraction at interim analyses
sfu <- sfLDOF # O'Brien-Fleming spending function for alpha-spending
x <- gsDesign(
  k = k,
  test.type = test.type,
  alpha = alpha,
  timing = timing,
  sfu = sfu
)
C1 = x$upper$bound[1]
C2 = x$upper$bound[2]
C1;C2

# expected cohen's d (theta*) for each endpoint (= expected theta for each endpoint since assume sigma = 1)
# under H0: mean theta=0 
exp_cohen_d = rep(0,n_endpts) ## assume for each endpoint, the effect size is 0

## Simulate 6 endpoints multivariate normal with correlation rho
sigma = 1
sd_trt = rep(sigma,n_endpts) # the sd_trt vector for 6 endpoints
sd_ct = rep(sigma,n_endpts) # the sd_ct vector for 6 endpoints
rho = .3

# build the covariance matrix
covariance = rho*sd_trt*sd_ct
cov_mat_trt = diag(sd_trt-rho) + matrix(covariance,ncol=n_endpts,nrow=n_endpts) 
cov_mat_ct = diag(sd_ct-rho) + matrix(covariance,ncol=n_endpts,nrow=n_endpts) 
cov_mat_trt
cov_mat_ct

mu_trt =  c(4,6,4,5,7,6)
mu_ct = c(4,6,4,5,7,6)

## set seed for reproducibility
seed1 = 1:1000000


# rejection at final analysis indicator:
sub_nsim = length(sub_sim_chunk)
rej_final_ind =rep(NA,sub_nsim)  # indicator of rejection, either interim or final
stop_at_interim_ind = rep(NA,sub_nsim) # if stop early 

# exact T statistics (T)
T_stat_interim =  rep(NA,sub_nsim)  
T_stat_stage2 =  rep(NA,sub_nsim) 
T_stat_realrho_interim =  rep(NA,sub_nsim) 
T_stat_realrho_stage2 =  rep(NA,sub_nsim) 

# normal test statistics (Z*) converted from T
test_stat_interim = rep(NA,sub_nsim) 
test_stat_stage2 = rep(NA,sub_nsim) 
test_stat_final = rep(NA,sub_nsim) 
test_stat_realrho_final = rep(NA,sub_nsim) 
test_stat_realrho_interim =  rep(NA,sub_nsim) 
test_stat_realrho_stage2=  rep(NA,sub_nsim) 

new_sample_size_trt = rep(NA,sub_nsim) 
new_sample_size_ct = rep(NA,sub_nsim) 

## for endpoint selection (not included anymore)
#endpoints_removed = as.data.frame(matrix(NA,nrow=sub_nsim,ncol=n_endpts))
#colnames(endpoints_removed) = paste0("end_pt_rm",1:n_endpts)

start = Sys.time()
for (s in (1:sub_nsim)){
  if (s%%50000 == 0) print(s)
  
  set.seed(seed1[s])
  
  ############# Data Generation - Stage 1
  y_trt_mat = rmvnorm(n_trt, mean = mu_trt,sigma=cov_mat_trt ) # covariance matrix
  y_ct_mat = rmvnorm(n_ct, mean = mu_ct, sigma=cov_mat_ct)  # covariance matrix
  
  # if endpoint selection, simulate baseline data for control group
  # if (end_select == 'yes'){
  #   # true percentage change
  #   percent_change = -c(.5, .5, .5, .5, 0.03, 0.01) # endpoint 5,6 are supposed to be removed.
  #   # assume the larger the endpoint value is better (baseline should be better than interim)
  #   # (interim-baseline)/baseline = pc <=> interim/baseline = 1+pc => baseline = interim/(1+pc)
  #   mu_ct_baseline = mu_ct/(1+percent_change)
  #   y_ct_mat0 = rmvnorm(n_ct, mean = mu_ct_baseline,
  #                       sigma=cov_mat_ct)
  #   ## colMeans((y_ct_mat-y_ct_mat0)/y_ct_mat0)
  # }
  # 
  #################################################### 
  ############# Interim Analysis
  
  # ####### 1. Endpoint Selection
  # if (end_select == 'yes'){
  #   ### compare observed means between interim and baseline for control group
  #   obs_perc_change = colMeans( (y_ct_mat-y_ct_mat0) / y_ct_mat0)
  #   ### selection with threshold = -20% percentage change
  #   if (any(obs_perc_change >= -.2)){
  #     var_keep_idx = which(obs_perc_change < -.2)
  #     var_remove_idx = obs_perc_change >= -.2
  #     endpoints_removed[s, var_remove_idx] = 1
  #     # remove the endpoints didn't get worse than baseline.
  #     y_trt1 = y_trt_mat[,var_keep_idx]
  #     y_ct1 = y_ct_mat[,var_keep_idx]
  #   }else{
  #     y_trt1 = y_trt_mat
  #     y_ct1 = y_ct_mat
  #     var_keep_idx = 1:n_endpts
  #   }
  # }else{ #  no endpoint selection 
  
  y_trt1 = y_trt_mat
  y_ct1 = y_ct_mat
  var_keep_idx = 1:n_endpts
  
  #}
  
  # remaining number of endpoints
  K = ncol(y_ct1)
  # expected effect size for remaining endpoints
  exp_cohen_d_K = exp_cohen_d[var_keep_idx] 
  mean_exp_cohen_d_K = mean(exp_cohen_d_K) # i.e., \bar{theta}^*, which equals 0 here.
  
  # The true SD of z_bar (mean z-score)
  A_exp = sqrt( 1/K^2 * (K + K*(K-1)*rho) )
  
  ## pooled sample variance matrix 
  S1_pool = (cov(y_trt1)*(n_trt-1) + cov(y_ct1)*(n_ct-1) )/ (n_trt+n_ct-2)
  ## get terms with sigma_hat and rho_hat from pooled sample variance
  sigma_hat = sqrt(diag(S1_pool))
  rho_ij = (S1_pool - diag(S1_pool)*diag(1,K)) / (matrix(sigma_hat,K,1)%*%sigma_hat)
  sum_rho = sum(rho_ij)# K*(K-1)*.3#sum(rho_ij)
  
  # The (estimated) SD of z_bar (mean z-score)
  A_obs = sqrt( 1/K^2 * (K + sum_rho) )
  
  ## cohen_d = (bar_Y_trt1_k - bar_Y_ct1_k)/sigma_hat1_k
  cohen_d =  (colMeans(y_trt1) - colMeans(y_ct1))/sigma_hat
  mean_cohen_d = mean(cohen_d) 
  
  # Standardized mean observed Z-score across all endpoints (T_1)
  bar_T_star1_realrho = sqrt(n_ct/(1/alloc_rate+1))*mean_cohen_d / A_exp # T_1 using real SD
  bar_T_star1 = sqrt(n_ct/(1/alloc_rate+1))*mean_cohen_d / A_obs # T_1 using estimated SD
  T_stat_realrho_interim[s] = bar_T_star1_realrho 
  T_stat_interim[s] = bar_T_star1
  # df proposed by logan and Tamhane 2024
  new_df = 0.5*(n_ct*(1+alloc_rate )-2)*(1+1/K^2)
  #bar_Z_star1 = qnorm(1-pt(bar_T_star1,df=n_ct*(1+alloc_rate )-2,lower.tail = F ))
  
  ## corresponding critical value in the standard normal distribution
  bar_Z_star1 = qnorm(1-pt(bar_T_star1,df=new_df,lower.tail = F )) 
  
  test_stat_interim[s] = bar_Z_star1
  test_stat_realrho_interim[s] = qnorm(1-pt(bar_T_star1_realrho,df=new_df,lower.tail = F )) # bar_Z_star1 under real SD
  
  #if (bar_Z_star1 < Cfut1){ # do not reject H0 at interim, stop trial
  # rej_final_ind[s] = 0 
  # stop_at_interim_ind[s] = 1 # stop b/c futility
  
  #}else 
  if (bar_Z_star1 > C1){ # reject H0 at interim, stop trial
    rej_final_ind[s] = 1
    stop_at_interim_ind[s] = 2 # stop b/c efficacy
    
  }else{ # bar_Z_star1 <= C1, do not reject H0
    stop_at_interim_ind[s] = 0
    
    ####### 4. Reestimate the sample size

    # when mean cohen_d <= 0, adjust it to be 0.01 to avoid denominator = 0 in SSR formula
    mean_cohen_d_adj = max(mean_cohen_d,.001)
    
    if (SSR=="yes"){
      n_hat = ceiling(A_obs^2 * ((C1+qnorm(1-beta0))*sqrt(1/alloc_rate+1)/mean_cohen_d_adj )^2 )
      # Or equivalently (as in the manuscript), 
      ## n_tot = A_obs^2 * ((C1+qnorm(1-beta0))*sqrt((1/alloc_rate+1)*(1+alloc_rate))/mean_cohen_d_adj )^2 
      ## n_hat = ceiling(n_tot/(1+alloc_rate) )
      M_ct = min(max(2*n_hat,N_ct),N_ct_max)
    }else{
      M_ct = N_ct
    }
    m_ct =  M_ct-n_ct
    m_trt = alloc_rate*m_ct 
    
    # record the new sample sizes
    new_sample_size_ct[s] = M_ct
    new_sample_size_trt[s] = alloc_rate*M_ct
     
    # #################################################### 
    # ############# Data Generation - Stage 2
    mu_trt_K = mu_trt[var_keep_idx]
    mu_ct_K = mu_ct[var_keep_idx]
    cov_mat_trt_K = cov_mat_trt[var_keep_idx,var_keep_idx]
    cov_mat_ct_K = cov_mat_ct[var_keep_idx,var_keep_idx]
    
    y_trt2 = rmvnorm(m_trt, mean = mu_trt_K,sigma=cov_mat_trt_K) # covariance matrix
    y_ct2 = rmvnorm(m_ct, mean = mu_ct_K, sigma=cov_mat_ct_K)  # covariance matrix
    
    ####################################################
    ############# Final Analysis

    ####### 1. Calculate the stage-2 standardized mean observed Z-score across all endpoints (T_2)
    ## pooled sample variance matrix
    S2_pool = (cov(y_trt2)*(m_trt-1) + cov(y_ct2)*(m_ct-1) )/ (m_trt+m_ct-2)
    ## get terms with sigma_hat and rho_hat
    sigma_hat2 = sqrt(diag(S2_pool))
    rho_ij2 = (S2_pool - diag(S2_pool)*diag(1,K)) / (matrix(sigma_hat2,K,1)%*%sigma_hat2)
    sum_rho2 = sum(rho_ij2)
    A_obs2 = sqrt( 1/K^2 * (K + sum_rho2) )
    
    ## cohen_d = (bar_Y_trt1_k - bar_Y_ct1_k)/sigma_hat1_k
    cohen_d2 =  (colMeans(y_trt2) - colMeans(y_ct2))/sigma_hat2
    mean_cohen_d2 = mean(cohen_d2) 
    ## mean_cohen_d2 can be negative and don't need adjust it by max(mean_cohen_d,.001) here, 
    # since no sample size re-estimation followed
    
    # # Standardized mean observed Z-score across all endpoints (T_2)
    bar_T_star2 = sqrt(m_ct/(1/alloc_rate+1))*mean_cohen_d2 / A_obs2
    T_stat_stage2[s] = bar_T_star2
    T_stat_realrho_stage2[s] = sqrt(m_ct/(1/alloc_rate+1))*mean_cohen_d2/A_exp 
    
    #pt(bar_Z_star2,df=m_ct*(1+alloc_rate )-2,lower.tail = F )
    #bar_Z_star2 = qnorm(1- pt(bar_T_star2,df=m_ct*(1+alloc_rate )-2,lower.tail = F ) )
    
    # df proposed by logan and Tamhane 2024
    new_df2 = 0.5*(m_ct*(1+alloc_rate )-2)*(1+1/K^2)
    #bar_Z_star1 = qnorm(1-pt(bar_T_star1,df=n_ct*(1+alloc_rate )-2,lower.tail = F ))
   
    ## corresponding critical value in the standard normal distribution
    bar_Z_star2 = qnorm(1-pt(bar_T_star2,df=new_df2,lower.tail = F ))
    
    test_stat_stage2[s] = bar_Z_star2
    test_stat_realrho_stage2[s] = qnorm(1- pt( T_stat_realrho_stage2[s],df=new_df2,lower.tail = F ) )
    
    ####### 2. Get final test statistics resulting from the inverse normal method of combining independent normal test statistics
    bar_Z_star_final = sqrt(1/2)*(bar_Z_star1+bar_Z_star2  )
    test_stat_final[s] = bar_Z_star_final
    test_stat_realrho_final[s] = sqrt(1/2)*(test_stat_realrho_interim[s]+test_stat_realrho_stage2[s] )
    # reject at final analysis: if bar_Z_star_final > C2=Cfut2
    if (bar_Z_star_final>C2){
      rej_final_ind[s] = 1
    }else{  # T_final < C2=Cfut2 do not reject at final analysis
      rej_final_ind[s] = 0
    }
  }  
}  
end = Sys.time()

running_time = difftime(end,start,units = "hours")

output = tibble(sub_sim_chunk,
                rej_final_ind,
                stop_at_interim_ind,
                #endpoints_removed,
                new_sample_size_trt,
                new_sample_size_ct,
                test_stat_interim,
                T_stat_interim,
                test_stat_realrho_interim,
                T_stat_realrho_interim,
                
                test_stat_stage2,
                T_stat_stage2,
                test_stat_realrho_stage2,
                T_stat_realrho_stage2,
                
                test_stat_final,
                test_stat_realrho_final,
                
                running_time,
                seed1=seed1[sub_sim_chunk]
)

mean(output$rej_final_ind)
# save the whole table
write_csv(output,paste0(getwd(),"/result_alpha729/",'output_debug715_1e6newdf',scena,'.csv'))

alpha = mean(output$rej_final_ind)
# save the empirical type I error
write_csv(as.data.frame(alpha),paste0(getwd(),"/result_alpha729/",'debug715_1e6newdf',scena,'.csv'))

### Later these individual results will be combined by running "combine_results.R".
## The combined results are already saved in the "Result datasets" folder 
