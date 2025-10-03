
##### Function of outcome 2*2 probability #####
outcome_prob_func <- function(rho, pt, pe){
  Delta <- rho*sqrt(pt*(1-pt)*pe*(1-pe))
  p1 <- -Delta+(1-pt)*pe
  p2 <- Delta+(1-pt)*(1-pe)
  p3 <- Delta+(pt*pe)
  p4 <- -Delta+pt*(1-pe)
  p_cell <- c(p1,p2,p3,p4)
  return(p_cell)
}

utility.table1<- c(100, 40, 60, 0)
####



Scenarios <- list()
#------------------------------------------------------------------------
## Scenarios when low-grade toxicities are not predictive of DLT   ##
## For results in Table S4
#----------------------------------------------------------------------
#### Scenario 1 ####
lambda_d <- c(20,20,20,20,20)/30

prob_T <- c(0.1, 0.25, 0.5, 0.65, 0.7)
mu_T <- qnorm(prob_T)

prob_E <- c(0.05, 0.1, 0.15, 0.3, 0.45)
mu_E <- qnorm(prob_E)

utility <- c()
for(i in 1:length(prob_T)){
  utility[i] <- sum(outcome_prob_func(0.3,prob_T[i],prob_E[i]) * utility.table1)
}
Scenarios[[1]] <- list(lambda_d = lambda_d, mu_T = mu_T, mu_E = mu_E, 
                       prob_T=prob_T, prob_E = prob_E, utility = utility)


#### Scenario 2  ####
lambda_d <- c(20,20,20,20,20)/30
prob_T <- c(0.15, 0.3, 0.5, 0.6, 0.7)
mu_T <- qnorm(prob_T)
prob_E <- c(0.4, 0.42, 0.45, 0.55, 0.65)
mu_E <- qnorm(prob_E)
utility <- c()
for(i in 1:length(prob_T)){
  utility[i] <- sum(outcome_prob_func(0.3,prob_T[i],prob_E[i]) * utility.table1)
}
Scenarios[[2]] <- list(lambda_d = lambda_d, mu_T = mu_T, mu_E = mu_E,
                       prob_T=prob_T, prob_E = prob_E, utility = utility)




#### Scenario 3 ####
lambda_d <- c(20,20,20,20,20)/30
prob_T <-c(0.02,0.05,0.3,0.5,0.9)
mu_T <- qnorm(prob_T)

prob_E <- c(0.08, 0.1, 0.4, 0.42, 0.45)
mu_E <- qnorm(prob_E)
utility <- c()
for(i in 1:length(prob_T)){
  utility[i] <- sum(outcome_prob_func(0.3,prob_T[i],prob_E[i]) * utility.table1)
}
Scenarios[[3]] <- list(lambda_d = lambda_d, mu_T = mu_T, mu_E = mu_E, 
                       prob_T=prob_T, prob_E = prob_E, utility = utility)




# #### Scenario 4####
lambda_d <-  c(10,10,10,10,10)/30


prob_T <- c(0.02,0.04,0.08,0.14,0.35)
mu_T <- qnorm(prob_T)

prob_E <- c(0.05, 0.1,0.25,0.45,0.46)
mu_E <- qnorm(prob_E)

utility <- c()
for(i in 1:length(prob_T)){
  utility[i] <- sum(outcome_prob_func(0.3,prob_T[i],prob_E[i]) * utility.table1)
}
Scenarios[[4]] <- list(lambda_d = lambda_d, mu_T = mu_T, mu_E = mu_E,
                        prob_T=prob_T, prob_E = prob_E, utility = utility)



# #### Scenario 5  ####

lambda_d <-  c(10, 10, 10, 10, 10)/30

prob_T <- c(0.03, 0.05, 0.15, 0.25, 0.5)
mu_T <- qnorm(prob_T)

prob_E <- c(0.05, 0.15, 0.3, 0.5, 0.55)
mu_E <- qnorm(prob_E)

utility <- c()
for(i in 1:length(prob_T)){
  utility[i] <- sum(outcome_prob_func(0.3,prob_T[i],prob_E[i]) * utility.table1)
}
Scenarios[[5]] <- list(lambda_d = lambda_d, mu_T = mu_T, mu_E = mu_E,
                        prob_T=prob_T, prob_E = prob_E, utility = utility)


# #### Scenario 6 ####

lambda_d <- c(5,5,5,5,5)/30

prob_T <- c(0.05, 0.08, 0.12,0.18,0.25)
mu_T <- qnorm(prob_T)

prob_E <- c(0.1,0.15,0.2,0.4,0.6)
mu_E <- qnorm(prob_E)

utility <- c()
for(i in 1:length(prob_T)){
  utility[i] <- sum(outcome_prob_func(0.3,prob_T[i],prob_E[i]) * utility.table1)
}
Scenarios[[6]] <- list(lambda_d = lambda_d, mu_T = mu_T, mu_E = mu_E,
                        prob_T=prob_T, prob_E = prob_E, utility = utility)

