###################################################################################################################
#                      Bayesian Phase I-II Designs with Adaptive Rules for Staggering Patient Entry
#             
#  This file R code to replicate Simulation results for "BOIN12" design in Table S2 
###################################################################################################################

source("Sim_codes_functions/Sim.func.BOIN12.R")
##scenarios
source("Sim_codes_functions/Sim.scenarios.051825.R")

####       Data generation          ####
ndose = 5
## start dose
startdose <- 1
## maximum sample size
max_n <- 45

## cohort size 
cohortsize=3

## toxicity upper bound 
pi_T_bar <- 0.35

## efficacy lower bound
pi_E_bar <- 0.25

## decision boundaries
p.saf=0.6*pi_T_bar
p.tox=1.4*pi_T_bar

## toxicity/efficacy Assessment window
DLT_window <- 30
efficacy_window <- 30

## required staggering window
stagger_window <- 30

## patients accrual rate
accrual.rate <- 5   ## average 5 patients per month
## distribution used to generate arrival time: accept "Uniform" and "Exponential"
arrive.dist =  "Uniform"
# distribution used to generate time to event: accept "Uniform" and "Weibull"
event.dist = "Uniform"

### utility 
utility.table <- c(60,100, 0, 40)  ## U(1,1), U(0,1), U(1,0), U(0,0)


## correlation between two latent outcomes
rho=0.3
Sigma <- matrix(c(1,rho,rho,1), ncol=2)   #covariance matrix of the variables.

## Early stiop if n.earlystop patients has been treated in a dose
n.earlystop = max_n
## Additional rules
N1=6
N2=9

## Number of simulation 
Nsim = 5000


## admissible cutoffs
Sce.list <- length(Scenarios)
scen_num <- seq(from = 1, to=Sce.list,by=1)

cut.tox <- c(0.9)
cut.eff <- c(0.95)

eta_elicited_rate <- c(0,0,0,0)

Final_result <- c()
for(sce in 1:Sce.list){
  ## Low-grade toxicity intensity with dose 
  lambda_d <- Scenarios[[sce]]$lambda_d
  
  ## latent variable mean 
  pT.true <- Scenarios[[sce]]$prob_T

  ## dose is multipied by 10 to stable numeric search, later can try without 10 version 
  pE.true <- Scenarios[[sce]]$prob_E
  
  BOIN12_design_sce <- BOIN12_Design(pT.true, pE.true, max_n,cohortsize,startdose,
                                     DLT_window, efficacy_window,stagger_window,
                                     pi_T_bar,pi_E_bar,utility.table, eta_elicited_rate,
                                     accrual.rate, arrive.dist, event.dist,
                                     lambda_d, rho,p.saf,p.tox,cut.tox,cut.eff,
                                      N1, N2,Nsim,sce)
 
  Final_result <- rbind(Final_result, BOIN12_design_sce)
}

Final_result1 <- data.frame(Final_result)
colnames(Final_result1) <- c("Scenario",paste("Dose ", c(0:ndose), sep = ""), "Trial Duration", "# pts treated","# Stagger" ,
                             "# stagger decision",
                             paste("# pts in Dose ", seq(1:ndose), sep = ""), "# DLT", "# Res","% stagger",
                             paste("% pts in Dose ", seq(1:ndose), sep = ""),"actual.Nsim")

Final_result1$Method = "BOIN12_nodecay"



Final_result_combine <- Final_result1[,c("Scenario","Method",paste("Dose ", c(0:ndose), sep = ""), "Trial Duration", "# pts treated",
                                         "# Stagger" ,"# stagger decision",
                                         paste("# pts in Dose ", seq(1:ndose), sep = ""),"% stagger",
                                         paste("% pts in Dose ", seq(1:ndose), sep = ""),"actual.Nsim","# DLT", "# Res")]

filename <-sprintf("Results/BOIN12_nodecay_ctox_%.2f_ceff_%.2f.csv", cut.tox,cut.eff)
write.csv(Final_result_combine,file=filename)










