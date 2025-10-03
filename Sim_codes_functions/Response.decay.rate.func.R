### Response rate decay function 
eta_elicited_func <- function(eta_elicited_rate, t_delay){
  eta_elicited_rate_all <- c(0,eta_elicited_rate)
  multiplier <- 1 - eta_elicited_rate_all  
  t_delay_week <-  t_delay/7
  
  h <- ceiling(t_delay_week)
  if(h==0){
    delay_latent_multiplier <- 1
  }else if(h > (length(eta_elicited_rate))){
    delay_latent_multiplier <- multiplier[h]
  }else{
    delay_latent_multiplier <- multiplier[h+1] + (multiplier[h+1] - multiplier[h])*(t_delay_week - h)
  }
  return(delay_latent_multiplier)
}

# y<- NULL
# for(i in 0:28){
#   y <- c(y,eta_elicited_func(eta_elicited_rate, i))
# }
# plot(0:28, y*0.4, type="l",ylim=c(0,0.5))
# dose <- c(0.1, 0.3, 0.5, 0.7, 0.9)
# ## number of doses
# ndose <- length(dose)
# ## standardize the dose to be mean 0 and standard deviation 0.5
# # Step 1: Subtract the mean
# centered_dose <- dose - mean(dose)
# # Step 2: Standardize to standard deviation of 1
# standardized_dose <- centered_dose / sd(centered_dose)* 0.5
# 
# ## dose is multipied by 10 to stable numeric search, later can try without 10 version 
# mu_E <- -1.3 + (1.09)/(1+exp((-1.13)*(10*standardized_dose-(-2.89)))) 
# # marginal efficacy rate is pnorm(coefficients[1,1] + (coefficients[3,1])/(1+exp((-coefficients[2,1])*(10*standardized_dose-(coefficients[4,1])))))
# # 0.1006302 0.2009631 0.4012644 0.4163788 0.4168211
# 
# ### response rate decay elicited from phycisions unit: weeks
# eta_elicited_rate <- c(0.05,0.1,0.18,0.3)
# 
# 
# ### transform the rate decay to latent scale unit: weeks 
# ## mu_E needs to be estimated from the model
# ## t_delay is the potential stagger time 
# eta_elicited_latent <- 1/mu_E*qnorm(eta_elicited_func(eta_elicited_rate, t_delay)*pnorm(mu_E))
# # mu_E
# # eta_elicited_latent
# # pnorm(mu_E)
# # pnorm(eta_elicited_latent*mu_E) = reponse rate * multiplier
# 
# 
# 
