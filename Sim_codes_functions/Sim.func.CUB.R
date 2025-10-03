#############################################################
#  Function to generate operating characteristics of 
#  Conventional Utility-based Design (CUB)  
#############################################################

####  Package loading     ####
library(MASS)
library(truncdist)
library(rjags)
library(R2jags)
library(runjags)
library(matrixcalc) ## check whether a matrix is positive definite
library(mvtnorm)
###############################################

source("Sim_codes_functions/Get.patients.tox.eff.outcome.R")
source("Sim_codes_functions/Response.decay.rate.func.R")
source("Sim_codes_functions/Functions.R")


#' @param dose the vector of  doses under investigation 
#' @param startdose the starting dose level for the trial
#' @param max_n maximum sample size
#' @param cohortsize the cohort size
#' @param DLT_window DLT Assessment window  
#' @param efficacy_window efficacy Assessment window
#' @param stagger_window required staggering window 
#' @param pi_T_bar the highest acceptable DLT rate
#' @param pi_E_bar the lowest acceptable efficacy rate
#' @param utility.table vector of utility of four outcomes 
#'                      (Y_T,Y_E) = c( (1,1), (0,1), (1,0), (0,1)  ) 
#' @param accrual.rate patients accrual rate, e.g., 
#'                     =5 means average 5 patients per month
#' @param arrive.dist  distribution used to generate arrival time: 
#'                     accept "Uniform" and "Exponential" 
#' @param event.dist  distribution used to generate time to event: 
#'                    accept "Uniform" and "Weibull" 
#' @param lambda_d  a vector of number of low grade toxicity intensity
#' @param rho the association between the latent bivariate normal random variables
#' @param p.saf the highest DLT probability that is deemed subtherapeutic
#'              (i.e. below the MTD) such that dose escalation should be undertaken.
#'              The default value is \code{p.saf=0.6*target}.
#' @param pT.true the true marginal toxicity rate vector.
#' @param pE.true the true marginal efficacy rate vector.
#' @param cut.tox the cutoff to eliminate an overly toxic dose for safety.
#'                  We recommend the default value of (\code{cut.tox=0.1}) for general use.
#'
#' @param cut.eff the cutoff to eliminate a futile dose.
#'                  We recommend the default value of (\code{cut.eff=0.05}) for general use.
#' @param cut.escalation Escalation cutoff, 
#'                       We recommend the default value of (\code{cut.escalation=0.5}) for general use   
#' @param Nsim the number of simulated trials.
#' @param sce indicator of scenario, in our simulations, we considered six scenarios.
#' @param subset We partition 5000 simulations to 40 subsets, this is an indicator of which subset is running.
#' @param File_name Prefix of the result file name to specify which design is considered.
#-------------------------------------
#' Models Prior parameters
# Prior specification details can be found in section 2.2 and S3
#' @Examples 
# beta0_bar =-1.2 
# tau_beta0 = 1/(4*beta0_bar)^2
# beta1_bar = 0.69
# tau_beta1 = 1/(4*beta1_bar)^2   ## prior parameter for mu_T
# alpha0_bar = -1.64
# tau_alpha0 = 1/(2*alpha0_bar)^2
# alpha1_bar = 0.2
# tau_alpha1 = 1/(2*alpha1_bar)^2
# alpha2_bar = 3.6
# tau_alpha2 = 1/(2*alpha2_bar)^2
# alpha3_bar = -1
# tau_alpha3 = 1/(2*alpha3_bar)^2  ## prior parameter for mu_E
#-------------------------------------

CUB_design <- function(dose, startdose, max_n, cohortsize, 
                           DLT_window, efficacy_window,stagger_window,
                           pi_T_bar, pi_E_bar, utility.table,
                           accrual.rate, arrive.dist, event.dist,
                           lambda_d, rho, pT.true, pE.true,
                           cut.tox,cut.eff, cut.escalation, Nsim,
                           beta0_bar, tau_beta0, beta1_bar, tau_beta1,  ## prior parameter for mu_T
                           alpha0_bar,tau_alpha0,alpha1_bar,tau_alpha1,alpha2_bar,tau_alpha2, ## prior parameter for mu_E
                           alpha3_bar, tau_alpha3,sce=1,subset=1,File_name="CUB"){
  
  s1<- Sys.time()
  #number of cohorts 
  ncohort = max_n/cohortsize
  ## latent variable mean
  mu_T <- qnorm(pT.true)
  ## dose is multipied by 10 to stable numeric search, later can try without 10 version
  mu_E <- qnorm(pE.true)

  Sigma <- matrix(c(1,rho,rho,1), ncol=2)   #covariance matrix of the variables.
  model.unsucc <- 0

  #### dose standerdization
  logdose <- log(dose)
  # ## number of doses
  ndose <- length(dose)
  # ## standardize the dose to be mean 0 and standard deviation 0.5
  # # Step 1: Subtract the mean
  centered_dose <- logdose - mean(logdose)
  # # Step 2: Standardize to standard deviation of 1
  standardized_dose <- centered_dose / sd(centered_dose)* 0.5
  dose = standardized_dose
  report_result <- c()
  result_all <- list()
  
  
  for(trial in 1:Nsim){
   # print(model.unsucc.idx)
   # if(trial>1){
   #   if(model.unsucc.idx == 1){
   #     trial = trial - 1
    #  }
   # }
    model.unsucc.idx <- 0
    set.seed(2024 + trial + (subset-1)*Nsim)
    
    ## patients enrolled in each dose
    npat_dose <- rep(0, ndose)  
    
    yT.true=yE.true=NULL;  #toxicity/efficacy indicator for each subject
    X.true = NULL;         #number of observed low-grade toxicity for each subject
    dose_each_patient <- NULL   #record dose for each subject
    
    t.enter=NULL; # time enter the study
    t.eff.event=t.tox.event=NULL; # time to event
    t.low.tox = NULL
    t.decision = 0; # decision making time
    
    n.low.tox <-NULL
    t.low.tox <-NULL
    
    ## indicator of whether the patients need to make stagger decision 
    idx_stagger <- NULL
    
    early.stop = 0
    ####  For each cohort   ####
    for(i in 1:ncohort){
      if(i == 1){
        current_dose = startdose
      }else{
        current_dose = next.dose
      }
      
 
      ## vectors to record information within each cohort 
      cohort.yT <- NULL
      cohort.t.dlt <- NULL
      cohort.yE <- NULL
      cohort.t.eff <- NULL
      cohort.pts.entry <- NULL
      cohort.t.low.tox <- NULL
      
      for(j in 1:cohortsize){
        if(j==1) { cohort.pts.entry = c(cohort.pts.entry, t.decision); } 
        else {
          if(arrive.dist=="Uniform"){ cohort.pts.entry = c(cohort.pts.entry, cohort.pts.entry[length(cohort.pts.entry)] + 30*runif(1, 0, 2/accrual.rate))}
          if(arrive.dist=="Exponential"){ cohort.pts.entry = c(cohort.pts.entry, cohort.pts.entry[length(cohort.pts.entry)] + 30*rexp(1, rate=accrual.rate))}
        }
      }
      t.enter = c(t.enter, cohort.pts.entry)
      
      
      ##Need to get this patient's outcome 
      pts.toxeff.out <- get.pts.tox.eff.outcome(d = current_dose, npts=cohortsize, mu_E, mu_T, Sigma, DLT_window, efficacy_window, dist=event.dist)
      
      ## Record patients tox and eff outcomes 
      yT.true = c(yT.true, pts.toxeff.out$tox)
      t.tox.event=c(t.tox.event, pts.toxeff.out$t.tox)
      
      yE.true = c(yE.true, pts.toxeff.out$eff)
      t.eff.event = c(t.eff.event, pts.toxeff.out$t.eff)
      
      
      ### # of stagger decisions 
      if(npat_dose[current_dose] ==0){
        pts.toxeff.out$t.tox
        
        idx_stagger <- c(idx_stagger,0)  ## first patient, no stagger 
        time.stagger.decition <- pts.toxeff.out$t.tox + cohort.pts.entry
        
        ## following patients 
        for(p in 2:cohortsize){
          idx_stagger <- c(idx_stagger,ifelse(any(cohort.pts.entry[p] < time.stagger.decition[1:(p-1)]), 1,0))
        }
        
      }else{idx_stagger <- c(idx_stagger, rep(0, cohortsize))}
      
      
      npat_dose[current_dose]<-npat_dose[current_dose] + cohortsize
      
      dose_each_patient = c(dose_each_patient, rep(current_dose, cohortsize))
      
      ## we need to wait until all patients' results are observed then make the escalation decision. 
      t.decision = max (cohort.pts.entry + pts.toxeff.out$t.tox)
      
      
      #follow up time for enrolled patients 
      followup.time <- t.decision - t.enter
      
      ## number of low-grade toxicity
      n.new.low.tox <- NULL
      
      ## The goal of low-grade toxicity is to predict DLT
      ## we only care about the number of low-grade toxicity before DLT assessment
      ## between cohorts always stagger, so no pending DLT assessments
      ## previous enrolled patients' number of low-grade toxicity won't update 
      n.low.tox <- c(n.low.tox, rpois(cohortsize,t.tox.event[(length(dose_each_patient) - cohortsize +1):length(dose_each_patient)]*lambda_d[current_dose]))
      t.low.tox <- c(t.low.tox, t.tox.event[(length(dose_each_patient) - cohortsize +1):length(dose_each_patient)])
      
      
      ### make the escalation decision 
      yT.true
      t.tox.event
      yE.true
      t.eff.event
      dose_each_patient
      t.low.tox
      n.low.tox
      omega.tox <-NULL
      omega.eff <- NULL
      yT.observe <- NULL
      yE.observe <- NULL
      cat_obs_data <- matrix(0,nrow = length(yT.true),ncol = 4)
      eta_elicited_piE <- NULL
      # actual_dose <- NULL 
      if(length(yT.true) == max_n){
        t.decision =  max(t.enter + pmax(t.tox.event, t.eff.event))
        followup.time = t.decision - t.enter
      }else{
        t.decision = t.decision
        followup.time = t.decision - t.enter
      }
      for(k in 1:length(yT.true)){
        
        ## observed outcomes at the decision point 
        yT.observe[k]  <- ifelse(followup.time[k] < t.tox.event[k], 0, yT.true[k])
        yE.observe[k] <- ifelse(followup.time[k]  < t.eff.event[k], 0, yE.true[k])
        
        ### create the multimomal dataframe: 1 <-> Pr(yT.observe =1, yE.observe=1)
        ###                                  2 <-> Pr(yT.observe =0, yE.observe=1)
        ###                                  3 <-> Pr(yT.observe =1, yE.observe=0)
        ###                                  4 <-> Pr(yT.observe =0, yE.observe=0)
        if(yT.observe[k] ==1 & yE.observe[k] ==1){
          cat_obs_data[k,1]=1
        }else if(yT.observe[k] ==0 & yE.observe[k] ==1){
          cat_obs_data[k,2]=1
        }else if(yT.observe[k] ==1 & yE.observe[k] ==0){
          cat_obs_data[k,3]=1
        }else{ cat_obs_data[k,4]=1}
        
        ## calculate the weights
        omega.tox[k] = ifelse(yT.observe[k]==1 | followup.time[k]>=DLT_window,1,followup.time[k]/DLT_window)
        omega.eff[k] = ifelse(yE.observe[k]==1 | followup.time[k]>=efficacy_window,1,followup.time[k]/efficacy_window)
        
        ## calculate the response decay 
        eta_elicited_piE[k] <- 1
        
        ## actual dose for each patient instead of dose level 
        # actual_dose[k] <- dose[dose_each_patient[k]] 
      }
      
        ### Run the Bayesian model to get the 
        jags.params<-c("rho","mu_T", "mu_E")#  parameters of interest
        jags.data<-list(dose = dose, ndose=ndose,dose_each_patient=dose_each_patient,  ## dose info
                        n.low.tox=n.low.tox,       ## low grade toxicity info
                        cat_obs_data=cat_obs_data, omega.tox=omega.tox, 
                        omega.eff=omega.eff,  eta_elicited_piE = eta_elicited_piE,     ## tox and eff info
                        beta0_bar=beta0_bar, tau_beta0=tau_beta0, beta1_bar=beta1_bar, tau_beta1=tau_beta1, ## prior parameter for mu_T
                        alpha0_bar=alpha0_bar, tau_alpha0 = tau_alpha0,alpha1_bar=alpha1_bar, tau_alpha1 = tau_alpha1,
                        alpha2_bar=alpha2_bar, tau_alpha2 = tau_alpha2,alpha3_bar=alpha3_bar, tau_alpha3 = tau_alpha3
        )
        
        jags.inits <- list(
          list("beta0"=0,"beta1"=0, "alpha0"=0,"alpha1"=0,"alpha2"=0,"alpha3"=0),
          list("beta0"=0,"beta1"=0, "alpha0"=0,"alpha1"=0,"alpha2"=0,"alpha3"=0),
          list("beta0"=0,"beta1"=0, "alpha0"=0,"alpha1"=0,"alpha2"=0,"alpha3"=0),
          list("beta0"=0,"beta1"=0, "alpha0"=0,"alpha1"=0,"alpha2"=0,"alpha3"=0)
        )
        model.unsucc.idx = 0
        runjagsmodel<- tryCatch({run.jags(model =eff_tox_model, monitor = jags.params,data = jags.data,
                                          n.chains = 4,adapt = 2000,burnin = 5000,inits = jags.inits,
                                          sample = 10000,summarise = FALSE,thin = 1,method = 'rjparallel',
                                          plots=FALSE,silent.jags = T)}, error=function(e){NULL})
        if(is.null(runjagsmodel)){
          model.unsucc.idx <- 1
          model.unsucc <- model.unsucc + 1
          break
        }
        
      
      codasamples<-as.mcmc.list(runjagsmodel)
      summary<-summary(codasamples)
      summary
      mcmc.sample <- as.matrix(codasamples)
      
      ### subtract the posterior mean for pi11,pi01,p10,p00
      post_mean_muT <- summary$statistics[substr(row.names(summary$statistics),1,4) =="mu_T","Mean"]
      post_mean_muE <- summary$statistics[substr(row.names(summary$statistics),1,4) =="mu_E","Mean"]
      
      post.rho <- summary$statistics[row.names(summary$statistics) =="rho","Mean"]
      post.Sigma <- matrix(c(1,post.rho,post.rho,1), ncol=2)   #covariance matrix of the variables.
      
      ## post joint tox and eff --> to calculate utility 
      post_mean_pi <- matrix(NA, nrow = ndose, ncol=4)
      for(u in 1:ndose){
        post_mean_pi[u,1] <- as.numeric(pmvnorm(lower=c(0,0), upper=c(Inf,Inf),mean=c(post_mean_muT[u], post_mean_muE[u]), corr = post.Sigma))
        post_mean_pi[u,2] <- as.numeric(pmvnorm(lower=c(-Inf,0), upper=c(0,Inf),mean=c(post_mean_muT[u], post_mean_muE[u]), corr = post.Sigma ))
        post_mean_pi[u,3] <- as.numeric(pmvnorm(lower=c(0,-Inf), upper=c(Inf,0),mean=c(post_mean_muT[u], post_mean_muE[u]), corr = post.Sigma ))
        post_mean_pi[u,4] <- as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(0,0),mean=c(post_mean_muT[u], post_mean_muE[u]), corr = post.Sigma ))
      }
      
      ## posterior mean of DLT rate and eff 
      post_mean_tox <-  pnorm(post_mean_muT)
      post_mean_eff  <- pnorm(post_mean_muE)
      
      ### subtract the posterior mean for pi11, pi01,p10,p00 given stagger time 
      # t_delay <- 10
      # post_mean_pi_delay <- delayed_joint_prob(mcmc.sample, 5, 30)
      # 
      
      ## first calculate the utility of all doses 
      post_utility <-  apply(post_mean_pi,1, function(x){
        sum(x*utility.table)
      } )
      
      
      
      ### first evaluate the safety of highest tried dose, if Pr(pi_T <= bar_pi_T) > Ces, 
      ##         and the dose is not the highest dose, then escalate
      #   Otherwise, identify the dose that maximizes mean utility among admissible doses 
      #             that are no higher than current dose 
      ### Admissible dose 
      ### The dose is admissible if (1) Pr(pi_T <= bar_pi_T) > cut.tox
      ###                       and (2) Pr(pi_E >= underlie_pi_E) > cut.eff
      highest_tried_dose <- max(dose_each_patient) 
      
      # pi.current.dose.ind <- c(paste("pi", "[",current_dose, ",1","]", sep = ""),
      #                          paste("pi", "[",current_dose, ",2","]", sep = ""),
      #                          paste("pi", "[",current_dose, ",3","]", sep = ""),
      #                          paste("pi", "[",current_dose, ",4","]", sep = ""))
      
      Prob_piT_leq_bar <- NULL
      Prob_piE_geq_bar <- NULL
      
      for(d in 1:ndose){
        
        pi.dose.ind <- c(paste("mu_T", "[",d,"]", sep = ""),
                         paste("mu_E", "[",d,"]", sep = ""))
        
        Prob_piT_leq_bar[d]  <- mean(pnorm(mcmc.sample[,pi.dose.ind[1]])<= pi_T_bar)
        Prob_piE_geq_bar[d]  <- mean(pnorm(mcmc.sample[,pi.dose.ind[2]]) >= pi_E_bar)
        
        # Prob_piT_leq_bar[d] <- mean(mcmc.sample[,pi.dose.ind[1]] + mcmc.sample[,pi.dose.ind[3]] <= pi_T_bar)
        # Prob_piE_geq_bar[d] <- mean(mcmc.sample[,pi.dose.ind[1]] + mcmc.sample[,pi.dose.ind[2]] >= pi_E_bar)
      }
      
      ## Pr(pi_T <= bar_pi_T|Data) 
      Prob_piT_leq_bar_highest <- Prob_piT_leq_bar[highest_tried_dose]
      
      if(Prob_piT_leq_bar_highest > cut.escalation & highest_tried_dose!= ndose){
        next.dose = highest_tried_dose + 1 
      }else{
        # prob_tox <- post_mean_pi[,1]+post_mean_pi[,3]
        # prob_eff <- post_mean_pi[,1]+post_mean_pi[,2]
        
        Admiss_ind <- Prob_piT_leq_bar>cut.tox & Prob_piE_geq_bar>cut.eff
        Admiss_set <- which(Admiss_ind == TRUE)
        
        # Admiss_set<- c(2,3,4)
        Candidate_set <- Admiss_set[Admiss_set <= highest_tried_dose]
        
        if(length(Candidate_set)==0){
          early.stop = 1 
          next.dose = 0
          print("No admissible dose and the trial is stopped.")
          break
        }
        
        ## adaptive randomization to lower admissible dose 
        de_cohort_AR_prob <- post_utility[Candidate_set]/sum(post_utility[Candidate_set])
        ind_dose_cohort <- rmultinom(1, 1, de_cohort_AR_prob)
        next.dose <- Candidate_set[which(ind_dose_cohort[,1]==1)]
        
        
      } 
      
      
    }  ## loop for cohorts 
    
    if(model.unsucc.idx == 0){
      ### dose selection 
      sim.data <- data.frame(dose_each_patient=dose_each_patient, yT.true = yT.true, yE.true = yE.true, n.low.tox=n.low.tox, t.low.tox=t.low.tox,idx_stagger=idx_stagger)
      
      est_tox <- post_mean_tox
      est_eff <- post_mean_eff
      
      
      sim.dose.info <- data.frame(est_tox=est_tox, est_eff = est_eff, post_utility=post_utility)
      
      Admiss_ind <- Prob_piT_leq_bar>cut.tox & Prob_piE_geq_bar>cut.eff
      Admiss_set <- which(Admiss_ind == TRUE & npat_dose!=0)
      ## OBD of this trial
      if(length(Admiss_set)==0){
        OBD = 0
      }else{
        OBD <- Admiss_set[which.max(post_utility[Admiss_set])]
      }
      
      ##trial length
      trial.length <- t.decision
      
      ## Number of patients staggerred 
      nstagger <- 0
      ntotal_stagger <- sum(idx_stagger)
      
      ntotal_pts <- length(yT.true)
        
      report_result <- rbind(report_result, c(OBD, trial.length, ntotal_pts,nstagger, ntotal_stagger,npat_dose))
      
      result_all[[trial]] <- list(sim.data, sim.dose.info, report_result)
    }
    
    print(paste("Simulation time: ", trial, sep=""))
  
  } 
  s2<- Sys.time()
  s2-s1
  report_result <- data.frame(report_result)
  colnames(report_result) <- c("OBD", "Trial Duration", "# pts", "# Stagger","# stagger decision" , paste("# pts in Dose ", seq(1:ndose), sep = "") )
  
  percent_select <- c()
  for(d in 1:(1+ndose)){
    percent_select[d] <- sum(report_result$OBD==(d-1))/(Nsim - model.unsucc)*100 
  }
  
  percent_select <- matrix(percent_select, byrow = TRUE, nrow=1)
  colnames(percent_select)<- paste("Dose ", seq(0:ndose)-1, sep = "") 
  print("Percentage of selection")
  print(percent_select)
  
  Average_result <- apply(report_result[,-1],2,mean)
  col_name <- names(Average_result)
  Average_result <- matrix(Average_result, byrow = TRUE, nrow=1)
  colnames(Average_result) <- col_name
  
  print(Average_result)
  
  pct_stagger <- mean(report_result[,4])/mean(report_result$`# stagger decision`)*100
  
  pct_pts_allo <- apply(report_result[,-c(1,2,3,4,5)],2,mean)/mean(report_result$`# pts`)*100
  pct_pts_allo <- matrix(pct_pts_allo, byrow = TRUE, nrow=1)
  colnames(pct_pts_allo) <- c(paste("% pts in Dose", 1:ndose,sep = ""))
  print(pct_pts_allo)
  
  result.out <- cbind(scenario = sce, subset=subset, percent_select, Average_result, `% stagger` = pct_stagger, pct_pts_allo, actual.Nsim = Nsim - model.unsucc)
  
  f1 <- sprintf("Results/%s_ctox_%.2f_ceff_%.2f_cesc_%.2f_scen_%d_subset_%d.csv", File_name,cut.tox,cut.eff, cut.escalation,sce,subset)
  
  write.csv(result.out,file=f1)
  
  return(result_all)
}

# ---------------------------------
# ####      Example           ####
#-----------------------------------
# # 
# dose <- c( 4e+7, 8e+7, 4e+8, 8e+8,4e+9)
# 
# ## True tox, eff,low-grade toxicity
# pT.true = c(0.1, 0.25, 0.5, 0.65, 0.7)
# pE.true =c(0.05, 0.1, 0.15, 0.3, 0.45)
# lambda_d = c(10,20,40,48,60)/30
# 
# ## start dose
# startdose <- 1
# ## maximum sample size
# max_n <- 45
# 
# ## cohort size
# cohortsize=3
# 
# # ## number of cohorts
# # ncohort = max_n/cohortsize
# 
# ## toxicity upper bound
# pi_T_bar <- 0.35
# 
# ## efficacy lower bound
# pi_E_bar <- 0.25
# 
# ## toxicity/efficacy Assessment window
# DLT_window <- 30
# efficacy_window <- 30
# 
# ## required staggering window
# stagger_window <- 30
# 
# ## patients accrual rate
# accrual.rate <- 5   ## average 5 patients per month
# ## distribution used to generate arrival time: accept "Uniform" and "Exponential"
# arrive.dist =  "Uniform"
# # distribution used to generate time to event: accept "Uniform" and "Weibull"
# event.dist = "Uniform"
# 
# ### utility
# utility.table <- c(60,100, 0, 40)  ## U(1,1), U(0,1), U(1,0), U(0,0)
# 
# ##Parameters for prior distribution
# beta0_bar =-1.2
# tau_beta0 = 1/(4*beta0_bar)^2
# beta1_bar = 0.69
# tau_beta1 = 1/(4*beta1_bar)^2
# alpha0_bar = -1.64
# tau_alpha0 = 1/(4*alpha0_bar)^2
# alpha1_bar = 0.2
# tau_alpha1 = 1/(4*alpha1_bar)^2
# alpha2_bar = 3.6
# tau_alpha2 = 1/(4*alpha2_bar)^2
# alpha3_bar = -1
# tau_alpha3 = 1/(4*alpha3_bar)^2
# 
# ## correlation between two latent outcomes
# rho=0.3
# Sigma <- matrix(c(1,rho,rho,1), ncol=2)   #covariance matrix of the variables.
# 
# ## Number of simulation
# Nsim = 5
# 
# 
# ## admissible cutoffs
# cut.tox <- c(0.1)
# cut.eff <- c(0.05)
# cut.escalation <- c(0.5)
# 
# 
# s1<- Sys.time()
# CUB_design_result <- CUB_design(dose, startdose, max_n, cohortsize,
#                                      DLT_window, efficacy_window,stagger_window,
#                                      pi_T_bar, pi_E_bar, utility.table,
#                                      accrual.rate, arrive.dist, event.dist,
#                                      lambda_d, rho, pT.true, pE.true,
#                                      cut.tox,cut.eff, cut.escalation, Nsim,
#                                      beta0_bar, tau_beta0, beta1_bar, tau_beta1,   ## prior parameter for mu_T
#                                      alpha0_bar,tau_alpha0,alpha1_bar,tau_alpha1,alpha2_bar,tau_alpha2, ## prior parameter for mu_E
#                                      alpha3_bar, tau_alpha3)
# 
# 
# 
# CUB_design_result
