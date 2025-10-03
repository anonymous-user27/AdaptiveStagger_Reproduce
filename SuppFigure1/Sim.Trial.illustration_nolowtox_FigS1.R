#-----------------------------------------------------------------------------------------------------------------#
#  This file R code to replicate Trial illustration in Supplementary Material S15, Figure S1 
#-----------------------------------------------------------------------------------------------------------------#

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


Adaptive_Stagger_onetime_nolowtox <- function(dose, startdose, max_n, cohortsize, 
                                    DLT_window, efficacy_window,  stagger_window,
                                    pi_T_bar, pi_E_bar, eta_elicited_rate, utility.table,
                                    accrual.rate, arrive.dist,event.dist,
                                    lambda_d, rho, pT.true, pE.true,
                                    cut.tox,cut.eff,cut.tox.stagger, cut.escalation, Nsim,
                                    beta0_bar, tau_beta0, beta1_bar, tau_beta1, beta2_bar, tau_beta2,  ## prior parameter for mu_T
                                    alpha0_bar,tau_alpha0,alpha1_bar,tau_alpha1,alpha2_bar,tau_alpha2, ## prior parameter for mu_E
                                    alpha3_bar, tau_alpha3){
  
  #number of cohorts 
  ncohort = max_n/cohortsize
  
  Sigma <- matrix(c(1,rho,rho,1), ncol=2)   #covariance matrix of the variables.
  model.unsucc <- 0
  #model.unsucc.idx <- 0
  
  ## latent variable mean
  mu_T <- qnorm(pT.true)
  ## dose is multipied by 10 to stable numeric search, later can try without 10 version
  mu_E <- qnorm(pE.true)
  
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

    
    ## patients enrolled in each dose
    npat_dose <- rep(0, ndose)  
    
    ##Count the number of stagger
    nstagger <- 0
    nstagger_dose <- rep(0,ndose)
    staggertime <-NULL
    
    yT.true=yE.true=NULL;  #toxicity/efficacy indicator for each subject
    X.true = NULL;         #number of observed low-grade toxicity for each subject
    dose_each_patient <- NULL   #record dose for each subject
    
    t.enter=NULL; # time enter the study
    t.eff.event=t.tox.event=NULL; # time to event
    t.low.tox = NULL
    t.decision = 0; # decision making time
    
    n.low.tox <-NULL
    t.low.tox <-NULL
    low.tox.ind <- NULL
    
    early.stop = 0
    
    ## indicator of whether the patients need to make stagger decision 
    idx_stagger <- NULL
    npat_initial_dose <- NULL
    
    
    ####  For each cohort   ####
    for(i in 1:ncohort){
      if(i == 1){
        current_dose = startdose
      }else{
        current_dose = next.dose
      }
      
      # if(early.stop == 1){
      #   break
      # }
      
      ## Each cohort would be assigned to a dose, if this dose has never been tried, then, we would apply the within cohort staggering rule
      if(npat_dose[current_dose] == 0){
        
        ## vectors to record information within each cohort 
        cohort.yT <- NULL
        cohort.t.dlt <- NULL
        cohort.yE <- NULL
        cohort.t.eff <- NULL
        cohort.pts.entry <- NULL
        cohort.follow <- NULL
        n.cohort.low.tox <- NULL
        cohort.t.low.tox <- NULL
        cohort.low.tox.ind <- NULL
        cohort.dose.info <- NULL
        
        for(j in 1:cohortsize){
          ## first patient in cohort is treated once arrive
          if(j==1){
            
            ## After every patient comes, we need to make a staggering decision, so the entry time for each subject is the decision time
            t.enter = c(t.enter, t.decision)
            t.pt.enter <- t.decision
            cohort.pts.entry = c(cohort.pts.entry,t.pt.enter) 
            
            staggertime = c(staggertime, 0)
            ##Need to get this patient's outcome 
            pts.toxeff.out <- get.pts.tox.eff.outcome(d = current_dose, npts=1, mu_E, mu_T, Sigma, DLT_window, efficacy_window, dist=event.dist)
            
            # DLT info
            pts.dlt <- pts.toxeff.out$tox
            pts.dlt.time <- pts.toxeff.out$t.tox
            
            # Response info
            pts.eff <- pts.toxeff.out$eff
            pts.eff.time <-  pts.toxeff.out$t.eff
            
            ## Record patients tox and eff outcomes 
            yT.true = c(yT.true, pts.dlt)
            t.tox.event=c(t.tox.event, pts.dlt.time)
            
            yE.true = c(yE.true, pts.eff)
            t.eff.event = c(t.eff.event, pts.eff.time)
            
            cohort.yT <- c(cohort.yT, pts.dlt)
            cohort.t.dlt <- c(cohort.t.dlt, pts.dlt.time)
            cohort.yE <- c(cohort.yE, pts.eff)
            cohort.t.eff <- c(cohort.t.eff, pts.eff.time)
            
            npat_dose[current_dose]<-npat_dose[current_dose] + 1
            
            dose_each_patient = c(dose_each_patient, current_dose)
            npat_initial_dose = c(npat_initial_dose, current_dose)
            idx_stagger = c(idx_stagger,0)
            
            ### when the following patients come, staggering rule will be applied
          }else{
            ## time of next patient come
            if(arrive.dist=="Uniform"){t.next.pt.enter = 30*runif(1, 0, 2/accrual.rate)}
            if(arrive.dist=="Exponential"){t.next.pt.enter = 30*rexp(1, rate=accrual.rate)}
            
            ##making staggering decision and the decision will be reapplied once DLT data is updated in this cohort 
            
            # when the new patient arrive, we need to make the decision of whether to stagger
            t.decision_stagger.first = t.pt.enter + t.next.pt.enter
            
            #follow up time for enrolled patients 
            followup.time <- t.decision_stagger.first - t.enter
            
            ## follow up time for patients in this cohort
            cohort.follow <- followup.time[((i-1)*cohortsize+1):length(followup.time)]
            
            ## dose assigned to patients in this cohort
            cohort.dose.info <- dose_each_patient[((i-1)*cohortsize+1):length(yT.true)]
            
            # indicator for whether or not toxicity and efficacy is observed in this cohort --  o means pending 
            delta.tox <- as.numeric((cohort.pts.entry + cohort.t.dlt) <= t.decision_stagger.first)
            delta.eff <- as.numeric((cohort.pts.entry + cohort.t.eff) <= t.decision_stagger.first)
            
            # count the number of applying staggering rule 
            num_stag_rule <- sum(delta.tox==0)+1
            
            pend.pts <- which(delta.tox==0)
            stag.decision.time <- c(t.decision_stagger.first,sort((cohort.pts.entry + cohort.t.dlt)[which(delta.tox==0)], decreasing = FALSE))
            pend.pts.order <- pend.pts[order((cohort.pts.entry + cohort.t.dlt)[which(delta.tox==0)],decreasing = FALSE)]
            
            
            ## stagger time for this patient 
            pt.max.stag.time <- max(DLT_window + 0.001 - cohort.follow[j-1],0)
            already.stag.time <- 0
            
            #### Need to write the loop to update stagger decision 
            for(r in 1:num_stag_rule){
              if(r==1){
                t.decision_stagger = t.decision_stagger.first
                followup.time = followup.time
                cohort.follow = cohort.follow
                
                if(is.null(n.cohort.low.tox)){
                  cohort.t.low.tox <- c(cohort.t.low.tox, min(cohort.follow, cohort.t.dlt))
                  n.cohort.low.tox <- c(n.cohort.low.tox,rpois(1, min(cohort.follow, cohort.t.dlt)*lambda_d[cohort.dose.info]))
                  cohort.low.tox.ind <- c(cohort.low.tox.ind, delta.tox)
                }else{
                  ## update the low grade toxicity for 1:(j-2) patient and get the lowgrade tox info for j-1 patient 
                  for(k in 1:(j-2)){
                    if(cohort.low.tox.ind[k] == 0  ){
                      cohort.t.low.tox[k] <- min(cohort.follow[k], cohort.t.dlt[k])
                      n.new.low.tox <- rpois(1, cohort.t.low.tox[k] * lambda_d[cohort.dose.info[k]])
                      n.cohort.low.tox[k] <- n.cohort.low.tox[k] + max(n.new.low.tox - n.cohort.low.tox[k],0)
                      cohort.low.tox.ind[k] <- ifelse(cohort.follow[k] >= cohort.t.dlt[k], 1,0 )
                    }else{n.cohort.low.tox[k] = n.cohort.low.tox[k]
                    cohort.t.low.tox[k] =  cohort.t.dlt[k]
                    cohort.low.tox.ind[k] = cohort.low.tox.ind[k]}
                  }
                  
                  cohort.t.low.tox <- c(cohort.t.low.tox, min(cohort.follow[j-1], cohort.t.dlt[j-1]))
                  n.cohort.low.tox <- c(n.cohort.low.tox,rpois(1, min(cohort.follow[j-1], cohort.t.dlt[j-1])*lambda_d[cohort.dose.info[j-1]]))
                  cohort.low.tox.ind <- c(cohort.low.tox.ind, ifelse(cohort.follow[j-1] >= cohort.t.dlt[j-1], 1,0 ))
                }
                
                if(i==1){
                  n.low.tox <- n.cohort.low.tox
                  t.low.tox <-  cohort.t.low.tox
                  low.tox.ind <- cohort.low.tox.ind
                }else{
                  n.low.tox <- c(n.low.tox[1:((i-1)*cohortsize)], n.cohort.low.tox)
                  t.low.tox <-  c(t.low.tox[1:((i-1)*cohortsize)], cohort.t.low.tox)
                  low.tox.ind <- c(low.tox.ind[1:((i-1)*cohortsize)],cohort.low.tox.ind)
                }
                already.stag.time <- 0
                
              }else{
                
                t.decision_stagger = stag.decision.time[r]
                followup.time = t.decision_stagger - t.enter
                
                cohort.follow <- followup.time[((i-1)*cohortsize+1):length(followup.time)]
                
                delta.tox <- as.numeric((cohort.pts.entry + cohort.t.dlt) <= t.decision_stagger)
                delta.eff <- as.numeric((cohort.pts.entry + cohort.t.eff) <= t.decision_stagger)
                
                already.stag.time <- t.decision_stagger - stag.decision.time[1]
                
                ## Update number of low-grade toxicity 
                
                for(k in 1:(j-1)){
                  if(cohort.low.tox.ind[k] == 0  ){
                    cohort.t.low.tox[k] <-  min(cohort.follow[k], cohort.t.dlt[k])
                    n.new.low.tox <- rpois(1, cohort.t.low.tox[k] * lambda_d[cohort.dose.info[k]])
                    n.cohort.low.tox[k] <- n.cohort.low.tox[k] + max(n.new.low.tox - n.cohort.low.tox[k],0)
                    cohort.low.tox.ind[k] <- ifelse(cohort.follow[k] >= cohort.t.dlt[k], 1,0 )
                  }else{n.cohort.low.tox[k] = n.cohort.low.tox[k]
                  cohort.t.low.tox[k] =  cohort.t.dlt[k]
                  cohort.low.tox.ind[k] = cohort.low.tox.ind[k]}
                }
                n.low.tox[((i-1)*cohortsize+1):length(n.low.tox)] <-  n.cohort.low.tox
                t.low.tox[((i-1)*cohortsize+1):length(t.low.tox)] <-  cohort.t.low.tox
                low.tox.ind[((i-1)*cohortsize+1):length(low.tox.ind)] <- cohort.low.tox.ind
              }
              
              ## evaluate whether the dose is admissible
              omega.tox <-NULL
              omega.eff <- NULL
              yT.observe <- NULL
              yE.observe <- NULL
              cat_obs_data <- matrix(0,nrow = length(yT.true),ncol = 4)
              eta_elicited_piE <- NULL
              # actual_dose <- NULL 
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
                eta_elicited_piE[k] <- eta_elicited_func(eta_elicited_rate, staggertime[k])
                
                ## actual dose for each patient instead of dose level 
                # actual_dose[k] <- dose[dose_each_patient[k]] 
              }
              
              ### Run the Bayesian model to get the 
              jags.params<-c( "rho","mu_T", "mu_E")#  parameters of interest
              jags.data<-list(dose = dose, ndose=ndose,dose_each_patient=dose_each_patient,  ## dose info
                              n.low.tox=n.low.tox,#t.low.tox=t.low.tox,       ## low grade toxicity info
                              cat_obs_data=cat_obs_data, omega.tox=omega.tox, 
                              omega.eff=omega.eff,  eta_elicited_piE = eta_elicited_piE,     ## tox and eff info
                              beta0_bar=beta0_bar, tau_beta0=tau_beta0, beta1_bar=beta1_bar, tau_beta1=tau_beta1,
                              # beta2_bar=beta2_bar, tau_beta2=tau_beta2,  ## prior parameter for mu_T
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
              runjagsmodel <- tryCatch({run.jags(model =eff_tox_model, monitor = jags.params,data = jags.data,
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
              
              # post_mean_pi <- summary$statistics[substr(row.names(summary$statistics),1,2) =="pi","Mean"]
              # post_mean_pi <- matrix(post_mean_pi, nrow=ndose,ncol=4) ## matrix: ndose * pi
              # post_mean_tox  <- post_mean_pi[,1] +  post_mean_pi[,3]
              # post_mean_eff  <- post_mean_pi[,1] +  post_mean_pi[,2]
              
              
              ####
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
              post_mean_eff <- pnorm(post_mean_muE)
              # post_mean_tox  <- post_mean_pi[,1] +  post_mean_pi[,3]
              # post_mean_eff  <- post_mean_pi[,1] +  post_mean_pi[,2]
              
              pi.current.dose.ind <- c(paste("mu_T", "[",current_dose,"]", sep = ""), paste("mu_E", "[",current_dose,"]", sep = ""))
              
              
              ## first calculate the utility of all doses 
              post_utility <-  apply(post_mean_pi,1, function(x){
                sum(x*utility.table)
              } )
              
              ### if DLT assessments are completed for patients in this cohort, 
              ##         (1) current dose is admissible, this patient is treated immediately at current dose
              ##         (2) not admissible, treated with dose lower dose with maximum mean utility or adaptive randomization 
              if(sum(delta.tox) == length(delta.tox)){
                ## Pr(pi_T <= bar_pi_T|Data) 
                
                Prob_piT_leq_bar_current  <- mean(pnorm(mcmc.sample[,pi.current.dose.ind[1]])<= pi_T_bar)
                #Prob_piE_geq_bar[d]  <- mean(pnorm(mcmc.sample[,pi.current.dose.ind[2]]) >= pi_E_bar)
                
                if(Prob_piT_leq_bar_current > cut.tox.stagger){
                  assign.dose = current_dose
                }else{
                  # prob_tox <- post_mean_pi[,1]+post_mean_pi[,3]
                  # prob_eff <- post_mean_pi[,1]+post_mean_pi[,2]
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
                  
                  Admiss_ind <- Prob_piT_leq_bar>cut.tox & Prob_piE_geq_bar>cut.eff
                  Admiss_set <- which(Admiss_ind == TRUE)
                  
                  # Admiss_set<- c(2,3,4)
                  Candidate_set <- Admiss_set[Admiss_set <= current_dose]
                  
                  ### treat at least one cohort of patients before terminate the trial
                  if(i==1){
                    assign.dose = current_dose
                  }else if(length(Candidate_set)==0){
                    early.stop = 1
                    staggertime <-c(staggertime, already.stag.time)
                    t.decision = t.decision_stagger 
                    print("No admissible dose and the trial is stopped.")
                    break
                  }else{
                    
                    ## de-escalate to the dose with maximum utility
                    assign.dose <- Candidate_set[which.max(post_utility[Candidate_set])]
                    
                  }
                  
                  
                }  
                
                
                ### Decision for this patient: enroll immediately, dose = assign.dose
                
                t.decision = t.decision_stagger
                
                t.enter = c(t.enter, t.decision)
                
                t.pt.enter <- t.decision
                
                cohort.pts.entry = c(cohort.pts.entry,t.pt.enter) 
                
                
                
                staggertime = c(staggertime, already.stag.time )
                
                mu_E_delay <- qnorm(eta_elicited_func(eta_elicited_rate, already.stag.time)*pnorm(mu_E))
                
                pts.toxeff.out <- get.pts.tox.eff.outcome(d = assign.dose, npts=1, mu_E = mu_E_delay, mu_T, Sigma, DLT_window, efficacy_window, dist=event.dist)
                
                # DLT info
                pts.dlt <- pts.toxeff.out$tox
                pts.dlt.time <- pts.toxeff.out$t.tox
                
                # Response info
                pts.eff <- pts.toxeff.out$eff
                pts.eff.time <-  pts.toxeff.out$t.eff
                
                ## Record patients tox and eff outcomes 
                yT.true = c(yT.true, pts.dlt)
                t.tox.event=c(t.tox.event, pts.dlt.time)
                
                yE.true = c(yE.true, pts.eff)
                t.eff.event = c(t.eff.event, pts.eff.time)
                
                cohort.yT <- c(cohort.yT, pts.dlt)
                cohort.t.dlt <- c(cohort.t.dlt, pts.dlt.time)
                cohort.yE <- c(cohort.yE, pts.eff)
                cohort.t.eff <- c(cohort.t.eff, pts.eff.time)
                
                npat_dose[assign.dose]<-npat_dose[assign.dose] + 1
                
                dose_each_patient = c(dose_each_patient, assign.dose)
                npat_initial_dose = c(npat_initial_dose, current_dose)
                idx_stagger = c(idx_stagger,ifelse(r==1,0,1))
                
              }else{
                
                ## if there are pending DLT, calculate expected utility 
                ## indicator of patients treated with  
                
                # cohort.t.dlt
                # cohort.follow
                # cohort.t.eff
                ## Lists possible outcome combination after staggering 
                cohort.dose <- dose_each_patient[((i-1)*cohortsize+1):length(yT.true)]
                
                pts_outcome <- data.frame(
                  ID = seq(1:length(cohort.yT)),
                  tox = ifelse(cohort.follow < cohort.t.dlt, 0, cohort.yT), ## observed outcome at the decision point
                  eff = ifelse(cohort.follow < cohort.t.eff, 0, cohort.yE), ## observed outcome at the decision point
                  pi_t = post_mean_tox[cohort.dose],
                  pi_e = post_mean_eff[cohort.dose],
                  delta_t = delta.tox,  ## 1: DLT assessment complete
                  delta_e =  delta.eff   ## 1: eff assessment complete 
                )
                
                pts_out_stagger <- mult_pts_stagger(pts_outcome)
                pts_out_stagger
                
                stagger_tox_name <- NULL
                stagger_eff_name <- NULL
                
                for(pd in pts_outcome$ID){
                  stagger_tox_name[pd] <-  paste("tox.pts",pd, sep="")
                  stagger_eff_name[pd] <-  paste("eff.pts",pd, sep="")
                }
                
                stag.utility.per.outcome <-NULL
                ##for each possible outcome combination, calculated the utility
                
                
                
                for(o in 1:dim(pts_out_stagger)[1]){
                  ##D_{i'-1}^S(a)
                  stagger_tox_outcome <- pts_out_stagger[o, stagger_tox_name]
                  stagger_eff_outcome <- pts_out_stagger[o, stagger_eff_name]
                  
                  ##D_{n+i'-1}^S(a)
                  if(i==1){
                    stagger_all_tox_out <- stagger_tox_outcome
                    stagger_all_eff_out <- stagger_eff_outcome
                  }else{
                    stagger_all_tox_out <- c(yT.true[1:((i-1)*cohortsize)], stagger_tox_outcome)
                    stagger_all_eff_out <- c(yE.true[1:((i-1)*cohortsize)], stagger_eff_outcome)
                  }
                  ### Fit the model using the data after staggering 
                  # n.low.tox 
                  # t.low.tox 
                  ## potential  stagger time for this patient
                  # pt.max.stag.time
                  
                  #follow up time for enrolled patients 
                  # followup.time 
                  
                  ## evaluate whether the dose is admissible
                  omega.tox.stag <-NULL
                  omega.eff.stag <- NULL
                  yT.observe.stag <- NULL
                  yE.observe.stag <- NULL
                  cat_obs_data.stag <- matrix(0,nrow = length(stagger_all_tox_out),ncol = 4)
                  eta_elicited_piE.stag <- NULL
                  
                  
                  # actual_dose <- NULL 
                  for(k in 1:length(stagger_all_tox_out)){
                    
                    ## observed outcomes at the decision point 
                    yT.observe.stag[k]  <-  stagger_all_tox_out[k]
                    yE.observe.stag[k]  <-  stagger_all_eff_out[k]
                    
                    ### create the multimomal dataframe: 1 <-> Pr(yT.observe =1, yE.observe=1)
                    ###                                  2 <-> Pr(yT.observe =0, yE.observe=1)
                    ###                                  3 <-> Pr(yT.observe =1, yE.observe=0)
                    ###                                  4 <-> Pr(yT.observe =0, yE.observe=0)
                    if( yT.observe.stag[k] == 1 & yE.observe.stag[k] == 1){
                      cat_obs_data.stag[k,1] = 1
                    }else if(yT.observe.stag[k] == 0 & yE.observe.stag[k] == 1){
                      cat_obs_data.stag[k,2] = 1
                    }else if(yT.observe.stag[k] == 1 & yE.observe.stag[k] == 0){
                      cat_obs_data.stag[k,3] = 1
                    }else{ cat_obs_data.stag[k,4] = 1}
                    
                    ## calculate the weights
                    omega.tox.stag[k] = 1
                    omega.eff.stag[k] = ifelse(stagger_all_eff_out[k]==1 | followup.time[k]+ pt.max.stag.time - already.stag.time >= efficacy_window,1,(followup.time[k]+ pt.max.stag.time - already.stag.time)/efficacy_window)
                    
                    ## calculate the response decay 
                    eta_elicited_piE.stag[k] <- eta_elicited_func(eta_elicited_rate, staggertime[k])
                    
                    ## actual dose for each patient instead of dose level 
                    # actual_dose[k] <- dose[dose_each_patient[k]] 
                  }
                  
                  ### Run the Bayesian model to get the 
                  jags.params <- c( "rho","mu_T", "mu_E")#  parameters of interest
                  jags.data <-list(dose = dose, ndose=ndose,dose_each_patient=dose_each_patient,  ## dose info
                                   n.low.tox=n.low.tox,#t.low.tox=t.low.tox,       ## low grade toxicity info
                                   cat_obs_data=cat_obs_data.stag, omega.tox=omega.tox.stag, 
                                   omega.eff=omega.eff.stag,  eta_elicited_piE = eta_elicited_piE.stag,     ## tox and eff info
                                   beta0_bar=beta0_bar, tau_beta0=tau_beta0, beta1_bar=beta1_bar, tau_beta1=tau_beta1,
                                   # beta2_bar=beta2_bar, tau_beta2=tau_beta2,  ## prior parameter for mu_T
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
                  runjagsmodel<-tryCatch({run.jags(model =eff_tox_model, monitor = jags.params,data = jags.data,
                                                   n.chains = 4,adapt = 2000,burnin = 5000,inits = jags.inits,
                                                   sample = 10000,summarise = FALSE,thin = 1,method = 'rjparallel',
                                                   plots=FALSE,silent.jags = T)}, error=function(e){NULL})
                  if(is.null(runjagsmodel)){
                    model.unsucc.idx <- 1
                    model.unsucc <- model.unsucc + 1
                    break
                  }
                  
                  codasamples.stag <-as.mcmc.list(runjagsmodel)
                  summary.stag<-summary(codasamples.stag)
                  summary.stag
                  mcmc.sample.stag <- as.matrix(codasamples.stag)
                  
                  ### subtract the posterior mean for pi11,pi01,p10,p00
                  
                  post_mean_muT.stag <- summary.stag$statistics[substr(row.names(summary.stag$statistics),1,4) =="mu_T","Mean"]
                  post_mean_muE.stag <- summary.stag$statistics[substr(row.names(summary.stag$statistics),1,4) =="mu_E","Mean"]
                  
                  post.rho.stag <- summary.stag$statistics[row.names(summary.stag$statistics) =="rho","Mean"]
                  post.Sigma.stag <- matrix(c(1,post.rho.stag,post.rho.stag,1), ncol=2)   #covariance matrix of the variables.
                  
                  
                  post_mean_tox.stag <-  pnorm(post_mean_muT.stag)
                  post_mean_eff.stag <- pnorm(post_mean_muE.stag)
                  
                  
                  
                  ## post joint tox and eff --> to calculate utility 
                  post_mean_pi.stag <- matrix(NA, nrow = ndose, ncol=4)
                  for(u in 1:ndose){
                    post_mean_pi.stag[u,1] <- as.numeric(pmvnorm(lower=c(0,0), upper=c(Inf,Inf),mean=c(post_mean_muT.stag[u], post_mean_muE.stag[u]), corr = post.Sigma.stag))
                    post_mean_pi.stag[u,2] <- as.numeric(pmvnorm(lower=c(-Inf,0), upper=c(0,Inf),mean=c(post_mean_muT.stag[u], post_mean_muE.stag[u]), corr = post.Sigma.stag ))
                    post_mean_pi.stag[u,3] <- as.numeric(pmvnorm(lower=c(0,-Inf), upper=c(Inf,0),mean=c(post_mean_muT.stag[u], post_mean_muE.stag[u]), corr = post.Sigma.stag ))
                    post_mean_pi.stag[u,4] <- as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(0,0),mean=c(post_mean_muT.stag[u], post_mean_muE.stag[u]), corr = post.Sigma.stag ))
                  }
                  
                  pi.current.dose.ind.stag <- c(paste("mu_T", "[",current_dose,"]", sep = ""),
                                                paste("mu_E", "[",current_dose,"]", sep = ""))
                  
                  
                  ### probability of de-escalation 
                  Prob_piT_leq_bar_stagger  <- mean(pnorm(mcmc.sample.stag[,pi.current.dose.ind.stag[1]])<= pi_T_bar)
                  
                  ### calculate the utility for each possible outcome 
                  ## subtract the posterior mean for pi11,pi01,p10,p00
                  pi_stagger <- delayed_joint_prob(post_mean_muT.stag, post_mean_muE.stag, post.rho.stag, ndose, t_delay = pt.max.stag.time)
                  
                  stagger_post_utility <- apply(pi_stagger,1, function(x){
                    sum(x*utility.table)
                  })
                  
                  
                  ## deescalation indicator, 1--de-escalate, 0-- stay
                  dees_ind <- as.numeric(Prob_piT_leq_bar_stagger <= cut.tox.stagger)
                  
                  
                  if(dees_ind == 1){
                    ###de-escalate to lower dose with highest utility in admissible set
                    ## identify admissible set 
                    Prob_piT_leq_bar.stag.all <- NULL
                    Prob_piE_geq_bar.stag.all <- NULL
                    
                    for(d in 1:ndose){
                      
                      pi.dose.ind <- c(paste("mu_T", "[",d,"]", sep = ""),
                                       paste("mu_E", "[",d,"]", sep = ""))
                      
                      Prob_piT_leq_bar.stag.all[d]  <- mean(pnorm(mcmc.sample.stag[,pi.dose.ind[1]])<= pi_T_bar)
                      Prob_piE_geq_bar.stag.all[d]  <- mean(pnorm(mcmc.sample.stag[,pi.dose.ind[2]]) >= pi_E_bar)
                      
                      # Prob_piT_leq_bar[d] <- mean(mcmc.sample[,pi.dose.ind[1]] + mcmc.sample[,pi.dose.ind[3]] <= pi_T_bar)
                      # Prob_piE_geq_bar[d] <- mean(mcmc.sample[,pi.dose.ind[1]] + mcmc.sample[,pi.dose.ind[2]] >= pi_E_bar)
                    }
                    
                    #Admiss_ind <- Prob_piT_leq_bar.stag.all>cut.tox & Prob_piE_geq_bar.stag.all>cut.eff
                    ## should consider both eff and tox
                    Admiss_safety_ind <- Prob_piT_leq_bar.stag.all > cut.tox & Prob_piE_geq_bar.stag.all>cut.eff
                    Admiss_safety_set <- which(Admiss_safety_ind == TRUE)
                    
                    ##### no admissible set, what to do 
                    Candidate_safety_set <- Admiss_safety_set[Admiss_safety_set <= current_dose]
                    
                    
                    if(length(Candidate_safety_set)==0){
                      ## no safe dose, the utility is the case y_T=1 and y_E=0
                      stag.utility.per.outcome[o] <- utility.table[3]*pts_out_stagger$all_prob[o]
                    }else{
                      ## de-escalate to the dose with maximum utility
                      opt_dose <- Candidate_safety_set[which.max(stagger_post_utility[Candidate_safety_set])]
                      stag.utility.per.outcome[o] <- stagger_post_utility[opt_dose]*pts_out_stagger$all_prob[o]
                      
                    }
                    
                    
                    
                    ### Utility of de-escalating dose 
                  } else{
                    ## treat at current dose   
                    stagger_post_utility <- apply(pi_stagger,1, function(x){
                      sum(x*utility.table)
                    })
                    stag.utility.per.outcome[o] <- stagger_post_utility[current_dose]*pts_out_stagger$all_prob[o]
                  }  
                  
                }## loop to get the utility based on each possible outcome combination
                if(model.unsucc.idx == 1){break}
                
                ## Utility after staggering 
                utility.stagger <- sum(stag.utility.per.outcome)
                
                ##compare it with the utility of non-staggering and then make stagger decision 
                ## first calculate the utility of all doses 
                pi_fur_stagger <- delayed_joint_prob(post_mean_muT, post_mean_muE, post.rho, ndose, t_delay = already.stag.time)
                post_utility_fur_stagger <-  apply(pi_fur_stagger,1, function(x){
                  sum(x*utility.table)})
                utility.non.stagger <- post_utility_fur_stagger[current_dose]
                
                stagger.ind <- ifelse(utility.stagger > utility.non.stagger + 1,1,0) ### eliminate randomness
                
                ## Prob_piT_leq_bar_current
                tox.idx  <- mean(pnorm(mcmc.sample[,pi.current.dose.ind[1]])<= pi_T_bar)
                #Prob_piE_geq_bar[d]  <- mean(pnorm(mcmc.sample[,pi.current.dose.ind[2]]) >= pi_E_bar)
                
                if(stagger.ind == 0 & tox.idx > cut.tox.stagger){
                  staggertime = c(staggertime, already.stag.time)
                  assign.dose = current_dose
                  ### Decision for this patient: enroll immediately, dose = assign.dose
                  
                  t.decision = t.decision_stagger
                  
                  t.enter = c(t.enter, t.decision)
                  
                  
                  t.pt.enter <- t.decision
                  
                  cohort.pts.entry = c(cohort.pts.entry,t.pt.enter) 
                  
                  
                  mu_E_delay <- qnorm(eta_elicited_func(eta_elicited_rate, already.stag.time)*pnorm(mu_E))
                  pts.toxeff.out <- get.pts.tox.eff.outcome(d = assign.dose, npts=1, mu_E=mu_E_delay, mu_T, Sigma, DLT_window, efficacy_window, dist=event.dist)
                  
                  # DLT info
                  pts.dlt <- pts.toxeff.out$tox
                  pts.dlt.time <- pts.toxeff.out$t.tox
                  
                  # Response info
                  pts.eff <- pts.toxeff.out$eff
                  pts.eff.time <-  pts.toxeff.out$t.eff
                  
                  ## Record patients tox and eff outcomes 
                  yT.true = c(yT.true, pts.dlt)
                  t.tox.event=c(t.tox.event, pts.dlt.time)
                  
                  yE.true = c(yE.true, pts.eff)
                  t.eff.event = c(t.eff.event, pts.eff.time)
                  
                  cohort.yT <- c(cohort.yT, pts.dlt)
                  cohort.t.dlt <- c(cohort.t.dlt, pts.dlt.time)
                  cohort.yE <- c(cohort.yE, pts.eff)
                  cohort.t.eff <- c(cohort.t.eff, pts.eff.time)
                  
                  npat_dose[assign.dose]<-npat_dose[assign.dose] + 1
                  
                  dose_each_patient = c(dose_each_patient, assign.dose)
                  npat_initial_dose = c(npat_initial_dose, current_dose)
                  idx_stagger = c(idx_stagger,1)
                  
                  ### No further staggering is needed. 
                  print("patient should be treated without further delay")
                  break
                }

                
              }### DLT pending - staggering decision
              
      
              
            }### staggering rule re-apply 
            
            if(early.stop == 1){break}  ## break the simulation of j 
            if(model.unsucc.idx == 1){break}
            
            ###   
          }
          
          
        }  
        if(model.unsucc.idx == 1){break}
        
        ### between cohort always stagger 
        if(early.stop == 1){
          t.decision = t.decision
        }else{
          t.decision = max (cohort.pts.entry + cohort.t.dlt)
        }
        
        
        #follow up time for enrolled patients 
        followup.time <- t.decision - t.enter
        
        ## follow up time for patients in this cohort
        cohort.follow <- followup.time[((i-1)*cohortsize+1):length(followup.time)]
        
        cohort.dose.info <- dose_each_patient[((i-1)*cohortsize+1):length(yT.true)]
        
        if(length(cohort.yT) == cohortsize){
          n.new.low.tox <- NULL
          for(k in 1:(cohortsize-1)){
            if(cohort.low.tox.ind[k] == 0  ){
              cohort.t.low.tox[k] <- min(cohort.follow[k], cohort.t.dlt[k])
              n.new.low.tox <- rpois(1, cohort.t.low.tox[k] * lambda_d[cohort.dose.info[k]])
              n.cohort.low.tox[k] <- n.cohort.low.tox[k] + max(n.new.low.tox - n.cohort.low.tox[k],0)
              cohort.low.tox.ind[k] <- ifelse(cohort.follow[k] >= cohort.t.dlt[k], 1,0 )
            }else{n.cohort.low.tox[k] = n.cohort.low.tox[k]
            cohort.t.low.tox[k] =  cohort.t.dlt[k]
            cohort.low.tox.ind[k] = cohort.low.tox.ind[k]}
          }
          cohort.t.low.tox <- c(cohort.t.low.tox, min(cohort.follow[cohortsize], cohort.t.dlt[cohortsize]))
          n.cohort.low.tox <- c(n.cohort.low.tox,rpois(1, min(cohort.follow[cohortsize], cohort.t.dlt[cohortsize])*lambda_d[cohort.dose.info[cohortsize]]))
          cohort.low.tox.ind <- c(cohort.low.tox.ind, ifelse(cohort.follow[cohortsize] >= cohort.t.dlt[cohortsize], 1,0 ))
          
          n.low.tox[((i-1)*cohortsize+1):length(yT.true)] <-  n.cohort.low.tox
          t.low.tox[((i-1)*cohortsize+1):length(yT.true)] <-  cohort.t.low.tox
          low.tox.ind[((i-1)*cohortsize+1):length(yT.true)] <- cohort.low.tox.ind
        }else{
          n.new.low.tox <- NULL
          for(k in 1:length(cohort.yT)){
            if(cohort.low.tox.ind[k] == 0  ){
              cohort.t.low.tox[k] <- min(cohort.follow[k], cohort.t.dlt[k])
              n.new.low.tox <- rpois(1, cohort.t.low.tox[k] * lambda_d[cohort.dose.info[k]])
              n.cohort.low.tox[k] <- n.cohort.low.tox[k] + max(n.new.low.tox - n.cohort.low.tox[k],0)
              cohort.low.tox.ind[k] <- ifelse(cohort.follow[k] >= cohort.t.dlt[k], 1,0 )
            }else{n.cohort.low.tox[k] = n.cohort.low.tox[k]
            cohort.t.low.tox[k] =  cohort.t.dlt[k]
            cohort.low.tox.ind[k] = cohort.low.tox.ind[k]}
          }
          n.low.tox[((i-1)*cohortsize+1):length(yT.true)] <-  n.cohort.low.tox
          t.low.tox[((i-1)*cohortsize+1):length(yT.true)] <-  cohort.t.low.tox
          low.tox.ind[((i-1)*cohortsize+1):length(yT.true)] <- cohort.low.tox.ind
        }
        
  
        ###### if there is no admissible dose in the lowest dose, then terminate the trial 
        if(early.stop == 1){
          next.dose = 0
          print("No admissible dose and the trial is stopped.")
          break
        }
        
      }else{
        ## if this dose has been tried, then, we assign all patients in this cohort to this dose
        
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
        staggertime = c(staggertime, rep(0, cohortsize))
        
        ##Need to get this patient's outcome 
        pts.toxeff.out <- get.pts.tox.eff.outcome(d = current_dose, npts=cohortsize, mu_E, mu_T, Sigma, DLT_window, efficacy_window, dist=event.dist)
        
        ## Record patients tox and eff outcomes 
        yT.true = c(yT.true, pts.toxeff.out$tox)
        t.tox.event=c(t.tox.event, pts.toxeff.out$t.tox)
        
        yE.true = c(yE.true, pts.toxeff.out$eff)
        t.eff.event = c(t.eff.event, pts.toxeff.out$t.eff)
        
        
        npat_dose[current_dose]<-npat_dose[current_dose] + cohortsize
        
        dose_each_patient = c(dose_each_patient, rep(current_dose, cohortsize))
        npat_initial_dose = c(npat_initial_dose, rep(current_dose, cohortsize))
        
        idx_stagger = c(idx_stagger,rep(0, cohortsize))
        
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
        low.tox.ind <- c(low.tox.ind, rep(1,cohortsize))
      }
      
      
      ### make the escalation decision 
      yT.true
      t.tox.event
      yE.true
      t.eff.event
      dose_each_patient
      staggertime
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
        eta_elicited_piE[k] <- eta_elicited_func(eta_elicited_rate, staggertime[k])
        
        ## actual dose for each patient instead of dose level 
        # actual_dose[k] <- dose[dose_each_patient[k]] 
      }
      
      ### Run the Bayesian model to get the 
      jags.params<-c( "rho","mu_T", "mu_E")#  parameters of interest
      jags.data<-list(dose = dose, ndose=ndose,dose_each_patient=dose_each_patient,  ## dose info
                      n.low.tox=n.low.tox,#t.low.tox=t.low.tox,       ## low grade toxicity info
                      cat_obs_data=cat_obs_data, omega.tox=omega.tox, 
                      omega.eff=omega.eff,  eta_elicited_piE = eta_elicited_piE,     ## tox and eff info
                      beta0_bar=beta0_bar, tau_beta0=tau_beta0, beta1_bar=beta1_bar, tau_beta1=tau_beta1,
                      # beta2_bar=beta2_bar, tau_beta2=tau_beta2,  ## prior parameter for mu_T
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
      runjagsmodel<-tryCatch({run.jags(model = eff_tox_model, monitor = jags.params,data = jags.data,
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
      post_mean_eff <- pnorm(post_mean_muE)
      
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
        
      } ##hybrid utility 2 
      
      
    }  ## loop for cohorts 
    
    
    if(model.unsucc.idx == 0){
      ### dose selection 
      sim.data <- data.frame(ID=1:length(dose_each_patient), Dose=dose_each_patient, t.enrollment=t.enter-staggertime,
                             t.stagger=staggertime, t.treated = t.enter,YT = yT.true, VT = t.tox.event,  
                             YE = yE.true, VE= t.eff.event, X=n.low.tox, 
                             idx_stagger = idx_stagger, `Initial dose`=npat_initial_dose)
      
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
      ##Number of patients in each dose 
      npat_dose
      ## Number of patients staggerred 
      nstagger_treated <- sum(staggertime[1:length(yT.true)]>0)
      
      ntotal_stagger_treated <- sum(idx_stagger)
      
      ntotal_pts_treated <- length(yT.true)
      
      
      
      report_result <- rbind(report_result, c(OBD, trial.length, ntotal_pts_treated,nstagger_treated, ntotal_stagger_treated,
                                              npat_dose))
      colnames(report_result) <- c("OBD", "Trial Duration", "# pts treated","# Stagger" ,"# stagger decision",
                                   paste("# pts in Dose ", seq(1:ndose), sep = "") )
      
      result_all <- list(sim.data, sim.dose.info, report_result,staggertime=staggertime)
    }
  
  return(result_all)
}
#------------------------------------------------------------------#
# Settings for Trial Illustration in Supplementary Material Figure S1
#--------------------------------------------------------------------#

## dose
dose <- c( 4e+7, 8e+7, 4e+8, 8e+8,4e+9)

## start dose
startdose <- 1
## maximum sample size
max_n <- 30

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
accrual.rate <- 5   ## average 5 patients per month
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
tau_alpha3 = 1/(4*alpha3_bar)^2


#pnorm(-1+0.2*dose+0.2*lambda_d)
#pnorm(-2+1/(1+exp(-1*(dose +1))))
## correlation between two latent outcomes
rho=0.3
Sigma <- matrix(c(1,rho,rho,1), ncol=2)   #covariance matrix of the variables.


### response rate decay elicited from physicians unit: weeks
eta_elicited_rate <- c(0.05,0.1,0.15,0.2)

cut.tox <- c(0.1)
cut.eff <- c(0.05)
cut.escalation <- c(0.5)
cut.tox.stagger <- c(0.4)





pT.true = c(0.03,0.05,0.15,0.25,0.50)
pE.true =c(0.05,0.15,0.30,0.50,0.55)
lambda_d = c(5,10,15,25,45)/30

set.seed(2237)




Adaptive_Stagger_onetime_nolowtox <- Adaptive_Stagger_onetime_nolowtox(dose , startdose, max_n, cohortsize,
                                                DLT_window, efficacy_window,  stagger_window,
                                                pi_T_bar, pi_E_bar, eta_elicited_rate, utility.table,
                                                accrual.rate, arrive.dist,event.dist,
                                                lambda_d, rho, pT.true, pE.true,
                                                cut.tox,cut.eff,cut.tox.stagger, cut.escalation, Nsim,
                                                beta0_bar, tau_beta0, beta1_bar, tau_beta1, beta2_bar, tau_beta2,  ## prior parameter for mu_T
                                                alpha0_bar,tau_alpha0,alpha1_bar,tau_alpha1,alpha2_bar,tau_alpha2, ## prior parameter for mu_E
                                                alpha3_bar, tau_alpha3)

Adaptive_Stagger_onetime_nolowtox

## round the numbers
df <- Adaptive_Stagger_onetime_nolowtox[[1]]
cols_to_round <- grep("^t\\.|^(VT|VE)$", names(df), value = TRUE)
df[cols_to_round] <- lapply(df[cols_to_round], round)

write.csv(df,"SuppFigure1/Trial_illustration_nolowtox_FigureS1.csv", row.names = FALSE)


