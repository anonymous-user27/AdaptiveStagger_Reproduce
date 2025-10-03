#----------------------------------------------------------------------------------------------#
####                           MCMC Bayesian Model -- regress on lambda                     ####
#----------------------------------------------------------------------------------------------#
joint_model_lambda="
model{
  #likelihood
  for(p in 1:length(n.low.tox)){
    
    ## Poisson likelihood for low-grade toxicities
    n.low.tox[p] ~ dpois(lambda[dose_each_patient[p]] * t.low.tox[p])
    
    ## multinomial likelihood for tox and eff
    cat_obs_data[p,] ~ dmulti(weighted_pi[p,], 1)
    weighted_pi[p,1] = omega.tox[p]*omega.eff[p]*pi_subject[p,1]
    weighted_pi[p,2] = (1-omega.tox[p])*omega.eff[p]*pi_subject[p,1] +  omega.eff[p]*pi_subject[p,2]
    weighted_pi[p,3] = omega.tox[p]*(1-omega.eff[p])*pi_subject[p,1] +  omega.tox[p]*pi_subject[p,3]
    weighted_pi[p,4]= 1-weighted_pi[p,1]-weighted_pi[p,2] - weighted_pi[p,3]
  }
  
  ### parameter related to Poisson distribution
  for(d in 1:ndose){
  lambda[d] <-  exp(gamma0 + gamma1 * dose[d])
  }
  
  ### parameter related to multinomial distribution
  ## probabilities for each subject 
  for(p in 1:length(n.low.tox)){
  
    mu_T_dose[p] = ((mu_T[dose_each_patient[p]] >=0) -(mu_T[dose_each_patient[p]]<=0))*(min(abs(mu_T[dose_each_patient[p]]),5))
    
    x[p] = max(-mu_T_dose[p],-mu_E_subject[p])
    y[p] = min(-mu_T_dose[p],-mu_E_subject[p])
    mu_a[p] = dnorm(x[p],0,1)/pnorm(-x[p],0,1)
    zeta_abrho[p] = (rho*mu_a[p] - y[p])/sqrt(1-rho^2)
    sigma_2a[p] = 1+x[p]*mu_a[p]-mu_a[p]^2
    pi_subject[p,1] = pnorm(-x[p],0,1)*(pnorm(zeta_abrho[p],0,1) - 1/2*rho^2/(1-rho^2)*zeta_abrho[p]*dnorm(zeta_abrho[p],0,1)*sigma_2a[p])
    pi_subject[p,2] = pnorm(mu_E_subject[p],0,1) - pi_subject[p,1]
    pi_subject[p,3] = pnorm(mu_T_dose[p],0,1) - pi_subject[p,1]
    pi_subject[p,4] = pi_subject[p,1] -1 +pnorm(x[p],0,1) +pnorm(y[p],0,1)
  }
  
  
  # probabilities for each dose without any delay
  ## parameter related to multinomial distribution
  # for(d in 1:ndose){
  #   x_dose[d] = max(-mu_T[d],-mu_E[d])
  #   y_dose[d] = min(-mu_T[d],-mu_E[d])
  #   mu_a_dose[d] = dnorm(x_dose[d],0,1)/pnorm(max(-x_dose[d],-5),0,1)
  #   zeta_abrho_dose[d] = (rho*mu_a_dose[d] - y_dose[d])/sqrt(1-rho^2)
  #   sigma_2a_dose[d] = 1+x_dose[d]*mu_a_dose[d]-mu_a_dose[d]^2
  #   pi[d,1] = pnorm(-x_dose[d],0,1)*(pnorm(zeta_abrho_dose[d],0,1) - 1/2*rho^2/(1-rho^2)*zeta_abrho_dose[d]*dnorm(zeta_abrho_dose[d],0,1)*sigma_2a_dose[d])
  #   pi[d,2] = pnorm(mu_E[d],0,1) - pi[d,1]
  #   pi[d,3] = pnorm(mu_T[d],0,1) - pi[d,1]
  #   pi[d,4] = pi[d,1] -1 +pnorm(x_dose[d],0,1) +pnorm(y_dose[d],0,1)
  # }
    
  ## calculate the muE for each subject accounting the delay 
  for(p in 1:length(n.low.tox)){
    mu_E_pts[p] = ((mu_E[dose_each_patient[p]] >=0) -(mu_E[dose_each_patient[p]]<=0))*(min(abs(mu_E[dose_each_patient[p]]),5))
    mu_E_subject[p] <-qnorm(eta_elicited_piE[p]*pnorm(mu_E_pts[p],0,1),0,1) 
  }
  
  ## latent mean 
  for(d in 1:ndose){
    mu_T[d]= beta0 + beta1*dose[d]+beta2*lambda[d]
    mu_E[d] = alpha0 + alpha2/(1+exp(-alpha1*(10*dose[d]-alpha3)))
  }
  
  ### Priors 
  gamma0 ~ dnorm(0, 1/1.25^2)
  gamma1 ~ dnorm(0, 1/1.25^2)T(0,)
  
  beta0 ~ dnorm(beta0_bar, tau_beta0)
  beta1 ~ dnorm(beta1_bar, tau_beta1)T(0,)
  beta2 ~ dnorm(beta2_bar, tau_beta2)T(0,)
  
  alpha0 ~dnorm(alpha0_bar, tau_alpha0)
  alpha1 ~dnorm(alpha1_bar, tau_alpha1)T(0,)
  alpha2 ~dnorm(alpha2_bar, tau_alpha2)T(0,)
  alpha3 ~dnorm(alpha3_bar, tau_alpha3)
  
  sigma2_beta0 = 1/tau_beta0
  sigma2_beta1 = 1/tau_beta1
  sigma2_beta2 = 1/tau_beta2
  
  sigma2_alpha0 = 1/tau_alpha0
  sigma2_alpha1 = 1/tau_alpha1
  sigma2_alpha2 = 1/tau_alpha2
  sigma2_alpha3 = 1/tau_alpha3
  
  rho ~ dunif(0,1)
}
"
#--------------------------------------------------------------------------------------------#
####                           MCMC Bayesian Model-regress on X(t)/t                      ####
#--------------------------------------------------------------------------------------------#
joint_model ="
model{
  #likelihood
  for(p in 1:length(n.low.tox)){
    
    ## Poisson likelihood for low-grade toxicities
    n.low.tox[p] ~ dpois(lambda[dose_each_patient[p]] * t.low.tox[p])
    
    ## multinomial likelihood for tox and eff
    cat_obs_data[p,] ~ dmulti(weighted_pi[p,], 1)
    weighted_pi[p,1] = omega.tox[p]*omega.eff[p]*pi_subject[p,1]
    weighted_pi[p,2] = (1-omega.tox[p])*omega.eff[p]*pi_subject[p,1] +  omega.eff[p]*pi_subject[p,2]
    weighted_pi[p,3] = omega.tox[p]*(1-omega.eff[p])*pi_subject[p,1] +  omega.tox[p]*pi_subject[p,3]
    weighted_pi[p,4]= 1-weighted_pi[p,1]-weighted_pi[p,2] - weighted_pi[p,3]
  }
  
  ### parameter related to Poisson distribution
  for(d in 1:ndose){
  lambda[d] <-  exp(gamma0 + gamma1 * dose[d])
  }
  
  ### parameter related to multinomial distribution
  ## probabilities for each subject 
  for(p in 1:length(n.low.tox)){
  
    mu_T_subject[p] = ((mu_T_pts[p] >=0) -(mu_T_pts[p]<=0))*(min(abs(mu_T_pts[p]),5))
    
    x[p] = max(-mu_T_subject[p],-mu_E_subject[p])
    y[p] = min(-mu_T_subject[p],-mu_E_subject[p])
    mu_a[p] = dnorm(x[p],0,1)/pnorm(-x[p],0,1)
    zeta_abrho[p] = (rho*mu_a[p] - y[p])/sqrt(1-rho^2)
    sigma_2a[p] = 1+x[p]*mu_a[p]-mu_a[p]^2
    pi_subject[p,1] = pnorm(-x[p],0,1)*(pnorm(zeta_abrho[p],0,1) - 1/2*rho^2/(1-rho^2)*zeta_abrho[p]*dnorm(zeta_abrho[p],0,1)*sigma_2a[p])
    pi_subject[p,2] = pnorm(mu_E_subject[p],0,1) - pi_subject[p,1]
    pi_subject[p,3] = pnorm(mu_T_subject[p],0,1) - pi_subject[p,1]
    pi_subject[p,4] = pi_subject[p,1] -1 +pnorm(x[p],0,1) +pnorm(y[p],0,1)
  }
  
  
  # probabilities for each dose without any delay
  ## parameter related to multinomial distribution
  # for(d in 1:ndose){
  #   x_dose[d] = max(-mu_T[d],-mu_E[d])
  #   y_dose[d] = min(-mu_T[d],-mu_E[d])
  #   mu_a_dose[d] = dnorm(x_dose[d],0,1)/pnorm(max(-x_dose[d],-5),0,1)
  #   zeta_abrho_dose[d] = (rho*mu_a_dose[d] - y_dose[d])/sqrt(1-rho^2)
  #   sigma_2a_dose[d] = 1+x_dose[d]*mu_a_dose[d]-mu_a_dose[d]^2
  #   pi[d,1] = pnorm(-x_dose[d],0,1)*(pnorm(zeta_abrho_dose[d],0,1) - 1/2*rho^2/(1-rho^2)*zeta_abrho_dose[d]*dnorm(zeta_abrho_dose[d],0,1)*sigma_2a_dose[d])
  #   pi[d,2] = pnorm(mu_E[d],0,1) - pi[d,1]
  #   pi[d,3] = pnorm(mu_T[d],0,1) - pi[d,1]
  #   pi[d,4] = pi[d,1] -1 +pnorm(x_dose[d],0,1) +pnorm(y_dose[d],0,1)
  # }
    
  ## calculate the muE for each subject accounting the delay 
  for(p in 1:length(n.low.tox)){
    mu_E_pts[p] = ((mu_E[dose_each_patient[p]] >=0) -(mu_E[dose_each_patient[p]]<=0))*(min(abs(mu_E[dose_each_patient[p]]),5))
    mu_E_subject[p] <-qnorm(eta_elicited_piE[p]*pnorm(mu_E_pts[p],0,1),0,1)
  }
  
  ## calculate the muT for each subject accounting low grade toxicity
  for(p in 1:length(n.low.tox)){
    mu_T_pts[p] = beta0 + beta1*dose[dose_each_patient[p]]+beta2*n.low.tox[p]/t.low.tox[p]
  }
  ## latent mean 
  for(d in 1:ndose){
    mu_T[d]= beta0 + beta1*dose[d]+beta2*lambda[d]
    mu_E[d] = alpha0 + alpha2/(1+exp(-alpha1*(10*dose[d]-alpha3)))
  }
  
  ### Priors 
  gamma0 ~ dnorm(0, 1)
  gamma1 ~ dnorm(0, 1)T(0,)
  
  beta0 ~ dnorm(beta0_bar, tau_beta0)
  beta1 ~ dnorm(beta1_bar, tau_beta1)T(0,)
  beta2 ~ dnorm(beta2_bar, tau_beta2)T(0,)
  
  alpha0 ~dnorm(alpha0_bar, tau_alpha0)
  alpha1 ~dnorm(alpha1_bar, tau_alpha1)T(0,)
  alpha2 ~dnorm(alpha2_bar, tau_alpha2)T(0,)
  alpha3 ~dnorm(alpha3_bar, tau_alpha3)
  
  sigma2_beta0 = 1/tau_beta0
  sigma2_beta1 = 1/tau_beta1
  sigma2_beta2 = 1/tau_beta2
  
  sigma2_alpha0 = 1/tau_alpha0
  sigma2_alpha1 = 1/tau_alpha1
  sigma2_alpha2 = 1/tau_alpha2
  sigma2_alpha3 = 1/tau_alpha3
  
  rho ~ dunif(0,1)
}
"

#--------------------------------------------------------------------------#
####                  Model without low-grade toxicity                  ####
#--------------------------------------------------------------------------#
eff_tox_model="
model{
  #likelihood
  for(p in 1:length(n.low.tox)){
    ## multinomial likelihood for tox and eff
    cat_obs_data[p,] ~ dmulti(weighted_pi[p,], 1)
    weighted_pi[p,1] = omega.tox[p]*omega.eff[p]*pi_subject[p,1]
    weighted_pi[p,2] = (1-omega.tox[p])*omega.eff[p]*pi_subject[p,1] +  omega.eff[p]*pi_subject[p,2]
    weighted_pi[p,3] = omega.tox[p]*(1-omega.eff[p])*pi_subject[p,1] +  omega.tox[p]*pi_subject[p,3]
    weighted_pi[p,4]= 1-weighted_pi[p,1]-weighted_pi[p,2] - weighted_pi[p,3]
  }
  
  
  ### parameter related to multinomial distribution
  ## probabilities for each subject 
  for(p in 1:length(n.low.tox)){
    
    mu_T_dose[p] = ((mu_T[dose_each_patient[p]] >=0) -(mu_T[dose_each_patient[p]]<=0))*(min(abs(mu_T[dose_each_patient[p]]),5))
     
    x[p] = max(-mu_T_dose[p],-mu_E_subject[p])
    y[p] = min(-mu_T_dose[p],-mu_E_subject[p])
    
    
    mu_a[p] = dnorm(x[p],0,1)/pnorm(-x[p],0,1)
    zeta_abrho[p] = (rho*mu_a[p] - y[p])/sqrt(1-rho^2)
    sigma_2a[p] = 1+x[p]*mu_a[p]-mu_a[p]^2
    pi_subject[p,1] = pnorm(-x[p],0,1)*(pnorm(zeta_abrho[p],0,1) - 1/2*rho^2/(1-rho^2)*zeta_abrho[p]*dnorm(zeta_abrho[p],0,1)*sigma_2a[p])
    pi_subject[p,2] = pnorm(mu_E_subject[p],0,1) - pi_subject[p,1]
    pi_subject[p,3] = pnorm(mu_T_dose[p] ,0,1) - pi_subject[p,1]
    pi_subject[p,4] = pi_subject[p,1] -1 +pnorm(x[p],0,1) +pnorm(y[p],0,1)
  }
  
  
  ## probabilities for each dose without any delay
      ### parameter related to multinomial distribution
    # for(d in 1:ndose){
    #   x_dose[d] = max(-mu_T[d],-mu_E[d])
    #   y_dose[d] = min(-mu_T[d],-mu_E[d])
    #   mu_a_dose[d] = dnorm(x_dose[d],0,1)/pnorm(-x_dose[d],0,1)
    #   zeta_abrho_dose[d] = (rho*mu_a_dose[d] - y_dose[d])/sqrt(1-rho^2)
    #   sigma_2a_dose[d] = 1+x_dose[d]*mu_a_dose[d]-mu_a_dose[d]^2
    #   pi[d,1] = pnorm(-x_dose[d],0,1)*(pnorm(zeta_abrho_dose[d],0,1) - 1/2*rho^2/(1-rho^2)*zeta_abrho_dose[d]*dnorm(zeta_abrho_dose[d],0,1)*sigma_2a_dose[d])
    #   pi[d,2] = pnorm(-y_dose[d],0,1) - pi[d,1]
    #   pi[d,3] = pnorm(-x_dose[d],0,1) - pi[d,1]
    #   pi[d,4] = pi[d,1] -1 +pnorm(x_dose[d],0,1) +pnorm(y_dose[d],0,1)
    # }
    
  ## calculate the muE for each subject accounting the delay 
  for(p in 1:length(n.low.tox)){
    mu_E_pts[p] = ((mu_E[dose_each_patient[p]] >=0) -(mu_E[dose_each_patient[p]]<=0))*(min(abs(mu_E[dose_each_patient[p]]),5))
    mu_E_subject[p] <-qnorm(eta_elicited_piE[p]*pnorm(mu_E_pts[p],0,1),0,1) 
  }
  
  ## latent mean 
  for(d in 1:ndose){
    mu_T[d]= beta0 + beta1*dose[d]
    mu_E[d] = alpha0 + alpha2/(1+exp(-alpha1*(10*dose[d]-alpha3)))
  }
  
  ### Priors 
  
  beta0 ~ dnorm(beta0_bar, tau_beta0)
  beta1 ~ dnorm(beta1_bar, tau_beta1)T(0,)
  
  
  alpha0 ~dnorm(alpha0_bar, tau_alpha0)
  alpha1 ~dnorm(alpha1_bar, tau_alpha1)T(0,)
  alpha2 ~dnorm(alpha2_bar, tau_alpha2)T(0,)
  alpha3 ~dnorm(alpha3_bar, tau_alpha3)
  
  sigma2_beta0 = 1/tau_beta0
  sigma2_beta1 = 1/tau_beta1
  
  sigma2_alpha0 = 1/tau_alpha0
  sigma2_alpha1 = 1/tau_alpha1
  sigma2_alpha2 = 1/tau_alpha2
  sigma2_alpha3 = 1/tau_alpha3
  
  rho ~ dunif(0,1)
}
"

#--------------------------------------------------------------------------#
#### Function to substract delayed pi11,pi01,pi10,p00 from MCMC samples ####
#--------------------------------------------------------------------------#
# delayed_joint_prob <- function(mcmc.sample, ndose, t_delay){
#   mcmc.rho <- mcmc.sample[,"rho"]
#   post_mean_pi_delay <- matrix(NA, nrow=ndose,ncol=4)
#   
#   for(d in 1:ndose){
#     
#     pi.delay.dose <- matrix(0,nrow = dim(mcmc.sample)[1], ncol=4 ) 
#     muT_dose_ind <- paste("mu_T", "[",d,"]", sep = "")
#     muE_dose_ind <- paste("mu_E", "[",d,"]", sep = "")
#     
#     mu_T = mcmc.sample[,muT_dose_ind]
#     mu_E_delay = qnorm(eta_elicited_func(eta_elicited_rate, t_delay)*pnorm(mcmc.sample[,muE_dose_ind]))
#     
#     x_dose = apply(cbind(-mu_T,-mu_E_delay), 1, max)
#     y_dose = apply(cbind(-mu_T,-mu_E_delay), 1, min)
#     
#     mu_a_dose = dnorm(x_dose,0,1)/pnorm(-x_dose,0,1)
#     zeta_abrho_dose = (mcmc.rho*mu_a_dose - y_dose)/sqrt(1-mcmc.rho^2)
#     sigma_2a_dose = 1+x_dose*mu_a_dose-mu_a_dose^2
#     
#     pi.delay.dose[,1] = pnorm(-x_dose,0,1)*(pnorm(zeta_abrho_dose,0,1) - 1/2*mcmc.rho^2/(1-mcmc.rho^2)*zeta_abrho_dose*dnorm(zeta_abrho_dose,0,1)*sigma_2a_dose)
#     pi.delay.dose[,2] = pnorm(mu_E_delay,0,1) - pi.delay.dose[,1]
#     pi.delay.dose[,3] = pnorm(mu_T,0,1) - pi.delay.dose[,1]
#     pi.delay.dose[,4] = pi.delay.dose[,1] -1 +pnorm(x_dose,0,1) +pnorm(y_dose,0,1)
#     
#     post_mean_pi_delay[d,] <- apply(pi.delay.dose,2, function(x){mean(x,na.rm = TRUE)})
#   }
#   return(post_mean_pi_delay)
# }

delayed_joint_prob <- function(post_mean_muT.stag, post_mean_muE.stag, post.rho.stag, ndose, t_delay){
  
  post_mean_pi_delay <- matrix(NA, nrow=ndose,ncol=4)
  
  for(d in 1:ndose){
    
    mu_T_ini = post_mean_muT.stag[d]
    mu_T = sign(mu_T_ini)*min(abs(mu_T_ini),5)
    
    mu_E_ini0 = post_mean_muE.stag[d]
    mu_E_ini = sign(mu_E_ini0)*min(abs(mu_E_ini0),5)
    mu_E_delay = qnorm(eta_elicited_func(eta_elicited_rate, t_delay)*pnorm(mu_E_ini))
    
    x_dose = apply(cbind(-mu_T,-mu_E_delay), 1, max)
    y_dose = apply(cbind(-mu_T,-mu_E_delay), 1, min)
    
    
    
    mu_a_dose = dnorm(x_dose,0,1)/pnorm(-x_dose,0,1)
    zeta_abrho_dose = (post.rho.stag*mu_a_dose - y_dose)/sqrt(1-post.rho.stag^2)
    sigma_2a_dose = 1+x_dose*mu_a_dose-mu_a_dose^2
    
    post_mean_pi_delay[d,1] = pnorm(-x_dose,0,1)*(pnorm(zeta_abrho_dose,0,1) - 1/2*post.rho.stag^2/(1-post.rho.stag^2)*zeta_abrho_dose*dnorm(zeta_abrho_dose,0,1)*sigma_2a_dose)
    post_mean_pi_delay[d,2] = pnorm(mu_E_delay,0,1) - post_mean_pi_delay[d,1]
    post_mean_pi_delay[d,3] = pnorm(mu_T,0,1) - post_mean_pi_delay[d,1]
    post_mean_pi_delay[d,4] = post_mean_pi_delay[d,1] -1 +pnorm(x_dose,0,1) +pnorm(y_dose,0,1)
    
  }
  return(post_mean_pi_delay)
}
# 
# 
# as.numeric(pmvnorm(lower=c(0,0), upper=c(Inf,Inf),mean=c(mu_T, mu_E_delay), corr = post.Sigma.stag))
# as.numeric(pmvnorm(lower=c(-Inf,0), upper=c(0,Inf),mean=c(mu_T, mu_E_delay), corr = post.Sigma.stag ))
# as.numeric(pmvnorm(lower=c(0,-Inf), upper=c(Inf,0),mean=c(mu_T, mu_E_delay), corr = post.Sigma.stag ))
# as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(0,0),mean=c(mu_T, mu_E_delay), corr = post.Sigma.stag ))
# 




#delayed_joint_prob(mcmc.sample, 5, 30)
#--------------------------------------------------------------------------#
####     Function to generate possible outcomes after staggering        ####
#--------------------------------------------------------------------------#

### c_status = (tox, eff); 0 or 1; must be length 2
stagger_info <- function(c_status, pi_t, pi_e, delta_t, delta_e){
  if(delta_t==1){
    t_vec=c_status[1];pi_t_s=0
  }else{
    t_vec <- (c_status[1]):1
    if(length(t_vec)==2){pi_t_s=c(pi_t,1-pi_t)}else{pi_t_s=0}
  }
  if(delta_e==1){
    e_vec=c_status[2];pi_e_s=0
  }else{
    e_vec <- (c_status[2]):1
    if(length(e_vec)==2){pi_e_s=c(pi_e,1-pi_e)}else{pi_e_s=0}
  }
  all_status <- expand.grid(tox=t_vec,eff=e_vec)
  all_p <- as.vector(outer(1-pi_t_s,1-pi_e_s))
  result <- cbind(all_status,all_p)
  return(result)
}
### Input No tox = 0, No eff = 0
### Pi_t=0.8, Pi_e=0.4
### Results:
### 1. tox = 0 & eff = 0 => pi=(1-0.8)*(1-0.4) = 0.12
### 2. tox = 1 & eff = 0 => pi=0.8*(1-0.4)     = 0.48
### 3. tox = 0 & eff = 1 => pi=(1-0.8)*0.4     = 0.08
### 4. tox = 1 & eff = 1 => pi=0.8*0.4         = 0.32
#stagger_info(c_status = c(0,0), pi_t = 0.8, pi_e = 0.4, delta_t=1, delta_e=0)

mult_pts_stagger <- function(pts_outcome){
  ind_info <- lapply(1:nrow(pts_outcome), 
                     function(x){
                       stagger_info(
                         c_status = c(pts_outcome$tox[x], pts_outcome$eff[x]),
                         pi_t = pts_outcome$pi_t[x], pi_e = pts_outcome$pi_e[x],
                         delta_t = pts_outcome$delta_t[x], delta_e = pts_outcome$delta_e[x])})
  ind_info_prob <- lapply(ind_info, function(df) df$all_p)
  ind_info_outcomes <- lapply(seq_along(ind_info), function(i){
    df <- ind_info[[i]][,c("tox","eff")]
    names(df) <- paste(names(df), paste0("pts", i), sep = ".")
    df
  })
  
  combined_outcome <- Reduce(function(x, y) do.call(merge, list(x, y, by = NULL)), ind_info_outcomes)
  all_prob <- as.vector(Reduce(function(x,y) do.call(outer, list(x,y)), ind_info_prob))
  result <- data.frame(combined_outcome, all_prob = all_prob)
  return(result)
}
### Pts outcome data
# pts_outcome <- data.frame(
#   ID = c(1,2,3),
#   tox = c(1,0,0),
#   eff = c(0,0,1),
#   pi_t = c(0.2,0.8,0.2),
#   pi_e = c(0.6,0.6,0.6),
#   delta_t = c(0,0,1),
#   delta_e = c(0,1,0)
# )
# pts_out_stagger <- mult_pts_stagger(pts_outcome)

# Sigma = matrix(c(1,post.rho.stag,post.rho.stag,1),nrow=2)
# pmvnorm(lower=c(0,0), upper = c(Inf,Inf), mean =c(post_mean_muT.stag[d], post_mean_muE.stag[d]), corr = Sigma)
# pmvnorm(lower=c(-Inf,0), upper = c(0,Inf), mean =c(post_mean_muT.stag[d], post_mean_muE.stag[d]), corr = Sigma)
# pmvnorm(lower=c(0,-Inf), upper = c(Inf,0), mean =c(post_mean_muT.stag[d], post_mean_muE.stag[d]), corr = Sigma)
# pmvnorm(lower=c(-Inf,-Inf), upper = c(0,0), mean =c(post_mean_muT.stag[d], post_mean_muE.stag[d]), corr = Sigma)

#-------------------------------------------------------------------------------------------#
#   isotonic transformation using the pool adjacent violator algorithm (PAVA)           ####
#-------------------------------------------------------------------------------------------#
pava <- function(x, wt = rep(1, length(x))) {
  n <- length(x)
  if (n <= 1)
    return(x)
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol)))
      break
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}

#-------------------------------------------------------------------------------------------#
#   Function to get BOIN12 Boundaries                    ####
#-------------------------------------------------------------------------------------------#
get.boundary.BOIN12 <- function(target, targetE, ncohort, cohortsize=3, p.saf=NA, p.tox=NA,  cutoff.eli=0.95,
                         cutoff.eli.E=0.90)
{
  if(is.na(p.saf)) p.saf=0.6*target;
  if(is.na(p.tox)) p.tox=1.4*target;
  
  ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
  npts = ncohort*cohortsize;
  ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;elimE=NULL;
  for(n in (1:ncohort)*cohortsize)
  {
    error.min=3;
    for(m1 in 0:(n-1))
    {
      for(m2 in (m1+1):n)
      {
        
        error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
        error2 = 1-pbinom(m1, n, p.saf);
        error3 = pbinom(m2-1, n, p.tox);
        
        error=error1+error2+error3;
        if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
      }
    }
    ntrt = c(ntrt, n);
    b.e = c(b.e, cutoff1);
    b.d = c(b.d, cutoff2);
    
    elimineed=0; # indicating whether elimination is needed
    elimineedE=0
    if(n<3) { elim = c(elim, NA); elimE = c(elimE,NA)}  # require treating at least 3 patients before eliminating a dose
    else
    {
      for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
      {
        if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
      }
      if(elimineed==1) { elim = c(elim, ntox); }
      else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
      
      for(neff in n:0){
        if(pbeta(targetE,neff+1,n-neff+1)>cutoff.eli.E){elimineedE=1; break;}
      }
      if(elimineedE==1){elimE=c(elimE,neff)} else {elimE=c(elimE,NA)}
    }
  }
  for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
  boundaries = rbind(ntrt, elim, b.d, b.e,elimE);
  rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=",
                           "Deescalate if # of DLT >=",  "Escalate if # of DLT <=", "Eliminate if # of Eff <=");
  colnames(boundaries) = rep("", ncohort);
  
  return(boundaries);
}
