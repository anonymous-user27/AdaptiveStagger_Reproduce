###################################################################################################################
#                      Bayesian Phase I-II Designs with Adaptive Rules for Staggering Patient Entry
#             
#  This file R code to replicate Simulation results for "Strict Stagger" design in Table S8
###################################################################################################################

####  Package loading     ####
library(MASS)
library(truncdist)
library(rjags)
library(R2jags)
library(runjags)
library(matrixcalc) ## check whether a matrix is positive definite
library(mvtnorm)
###############################################

setwd("/rsrch5/scratch/biostatistics/swang23/Staggering")

source("Sim_codes_functions/Sim.func.Strict.Stagger.R")

##scenarios
source("Sim_codes_functions/Sim.scenarios.R")


args <- commandArgs(trailingOnly = TRUE)
id = as.numeric(args[1])



####       Data generation          ####
dose <- c( 4e+7, 8e+7, 4e+8, 8e+8,4e+9)


## start dose
startdose <- 1
## maximum sample size
max_n <- 45

## cohort size 
cohortsize=3

# ## number of cohorts
# ncohort = max_n/cohortsize

## toxicity upper bound 
pi_T_bar <- 0.35

## efficacy lower bound
pi_E_bar <- 0.25

## toxicity/efficacy Assessment window
DLT_window <- 30
efficacy_window <- 30

## required staggering window
stagger_window <- 30

## patients accrual rate
accrual.rate <- 3   ## average 3 patients per month
## distribution used to generate arrival time: accept "Uniform" and "Exponential"
arrive.dist =  "Uniform"
# distribution used to generate time to event: accept "Uniform" and "Weibull"
event.dist = "Uniform"


### utility 
utility.table <- c(60,100, 0, 40)  ## U(1,1), U(0,1), U(1,0), U(0,0)

##Parameters for prior distribution 
beta0_bar =-1.2
tau_beta0 = 1/(4*beta0_bar)^2
beta1_bar = 0.69
tau_beta1 = 1/(4*beta1_bar)^2 ## prior parameter for mu_T
alpha0_bar = -1.64
tau_alpha0 = 1/(4*alpha0_bar)^2
alpha1_bar = 0.2
tau_alpha1 = 1/(4*alpha1_bar)^2
alpha2_bar = 3.6
tau_alpha2 = 1/(4*alpha2_bar)^2
alpha3_bar = -1
tau_alpha3 = 1/(4*alpha3_bar)^2 ## prior parameter for mu_E


## correlation between two latent outcomes
rho=0.3
Sigma <- matrix(c(1,rho,rho,1), ncol=2)   #covariance matrix of the variables.

## Number of simulation 
Nsim = 50

Sce.list <- length(Scenarios)
scen_num <- seq(from = 1, to=Sce.list,by=1)

### response rate decay elicited from physicians unit: weeks
eta_elicited_rate <- c(0.05,0.1,0.15,0.2)


cut.tox <- c(0.1)
cut.eff <- c(0.05)
cut.escalation <- c(0.5)
cut.tox.stagger <- c(0.4)
  
subset.try <- c(1:50)

all_para <- expand.grid(scen_num, subset.try)
colnames(all_para) <- c("scenario", "subset")


sce = all_para[id,]$`scenario`
subset = all_para[id,]$`subset`


##
## Low-grade toxicity intensity with dose 
lambda_d <- Scenarios[[sce]]$lambda_d

##true marginal DLT and efficacy
pT.true <- Scenarios[[sce]]$prob_T
pE.true <- Scenarios[[sce]]$prob_E

File_name="Strict_Stagger_AR3"

s1 <- Sys.time()
Strict_Stagger_sce <- Strict_Stagger_design(dose, startdose, max_n, cohortsize,
                                              DLT_window, efficacy_window,  stagger_window,
                                              pi_T_bar, pi_E_bar, eta_elicited_rate, utility.table,
                                              accrual.rate, arrive.dist,event.dist,
                                              lambda_d, rho, pT.true, pE.true,
                                              cut.tox,cut.eff,cut.tox.stagger, cut.escalation,Nsim,
                                              beta0_bar, tau_beta0, beta1_bar, tau_beta1,  ## prior parameter for mu_T
                                              alpha0_bar,tau_alpha0,alpha1_bar,tau_alpha1,alpha2_bar,tau_alpha2, ## prior parameter for mu_E
                                              alpha3_bar, tau_alpha3, sce, subset,File_name)
#Strict_Stagger_sce
s2 <- Sys.time()
s2-s1
filename <-sprintf("Results/%s_ctox_%.2f_ceff_%.2f_ces_%.2f_cstag_%.2f_scen_%d_subset_%d.RData", File_name, cut.tox,cut.eff, cut.escalation,cut.tox.stagger,sce,subset)


save(Strict_Stagger_sce,file=filename)

