#############################################################
#  Function to generate operating characteristics of 
#  Backfill Design  
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
#' @param eta_elicited_rate the vector of response rate decay elicited 
#'                          from physicians unit: weeks   
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
#'                  We recommend the default value of (\code{cut.tox=0.05}) for general use.
#'
#' @param cut.eff the cutoff to eliminate a futile dose.
#'                  We recommend the default value of (\code{cut.eff=0.1}) for general use.
#' @param cut.escalation Escalation cutoff, 
#'                       We recommend the default value of (\code{cut.escalation=0.5}) for general use   
#' @param cut.tox.stagger Safety evaluation cutoff to decide whether stay at current dose after staggering
#'                        We recommend the default value of (\code{cut.tox.stagger=0.4}) for general use
#' @param Nsim the number of simulated trials.
#' @param max_acc maximum number of patients enrolled during staggering period, 
#'                We recommend the default value of (\code{max_acc=20}) for general use.
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



Backfill_design <- function(dose, startdose, max_n, cohortsize, 
                            DLT_window, efficacy_window,  stagger_window,
                            pi_T_bar, pi_E_bar, eta_elicited_rate, utility.table,
                            accrual.rate, arrive.dist,event.dist,
                            lambda_d, rho, pT.true, pE.true,
                            cut.tox,cut.eff,cut.tox.stagger, cut.escalation, Nsim,
                            beta0_bar, tau_beta0, beta1_bar, tau_beta1,  ## prior parameter for mu_T
                            alpha0_bar,tau_alpha0,alpha1_bar,tau_alpha1,alpha2_bar,tau_alpha2, ## prior parameter for mu_E
                            alpha3_bar, tau_alpha3,max_acc = 20, sce=1, subset=1, File_name="Backfill"){
  

  s1 <- Sys.time()
  #number of cohorts 
  ncohort = max_n/cohortsize
  
  Sigma <- matrix(c(1,rho,rho,1), ncol=2)   #covariance matrix of the variables.
  
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
  
  #model.unsucc.idx <- 0
  model.unsucc <- 0
  report_result <- c()
  result_all <- list()
  
  #fail_mark <- c()
  
  #trial = 1
  for(trial in 1:Nsim){
  
    set.seed(2024 + trial + (subset-1)*Nsim)
    
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
    
    ## indicator of whether the patients need to make stagger decision 
    idx_stagger <- NULL
    npat_initial_dose <- NULL
    
    early.stop = 0
    ####  For each cohort   ####
    for(i in 1:ncohort){
      if(i == 1){
        current_dose = startdose
        cohortsize_use = cohortsize
      }else{
        current_dose = next.dose
        cohortsize_use = cohortsize
        if(length(yT.true) + cohortsize >max_n){
          cohortsize_use = max_n - length(yT.true)
        }
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
        cohort.t.stagger <- NULL
        cohort.follow <- NULL
        n.cohort.low.tox <- NULL
        cohort.t.low.tox <- NULL
        cohort.low.tox.ind <- NULL
        cohort.dose.info <- NULL
        
        backfill.cohort.pts.entry <- NULL
        backfill.cohort.dose.info <- NULL
        
        
        for(j in 1:cohortsize_use){
          if(length(yT.true) == max_n){break}
          ## first patient in cohort is treated once arrive
          if(j==1){
            
            ## After every patient comes, we need to make a staggering decision, so the entry time for each subject is the decision time
            t.enter = c(t.enter, t.decision)
            t.pt.enter <- t.decision
            cohort.pts.entry = c(cohort.pts.entry,t.pt.enter) 
            
            staggertime = c(staggertime, 0)
            cohort.t.stagger = c(cohort.t.stagger,0)
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
            
            # cohort.t.low.tox <- c(cohort.t.low.tox, pts.dlt.time)
            # n.cohort.low.tox <- c(n.cohort.low.tox, rpois(1, pts.dlt.time * lambda_d[current_dose]))
            # 
            t.low.tox <- c(t.low.tox, pts.dlt.time)
            n.low.tox <- c(n.low.tox, rpois(1, pts.dlt.time * lambda_d[current_dose]))
            
            
            npat_dose[current_dose]<-npat_dose[current_dose] + 1
            
            dose_each_patient = c(dose_each_patient, current_dose)
            npat_initial_dose = c(npat_initial_dose, current_dose)
            cohort.dose.info <- c(cohort.dose.info, current_dose)
            idx_stagger = c(idx_stagger,0)
            
            ### when the following patients come, staggering rule will be applied
          }else{
            ## time of next patient come
            
            if(arrive.dist=="Uniform"){t.next.pt.enter = 30*runif(1, 0, 2/accrual.rate)}
            if(arrive.dist=="Exponential"){t.next.pt.enter = 30*rexp(1, rate=accrual.rate)}
            
            ##making staggering decision and the decision will be reapplied once DLT data is updated in this cohort 
            
            # when the new patient arrive, we need to make the decision of whether to stagger
            t.decision_stagger = t.pt.enter + t.next.pt.enter
            
            #follow up time for enrolled patients 
            followup.time <- t.decision_stagger - t.enter
            
            
            ## follow up time for patients in this cohort
            cohort.follow <- t.decision_stagger- cohort.pts.entry
            
            ## dose assigned to patients in this cohort
            cohort.dose.info
            
            ## stagger time for this patient 
            pt.max.stag.time <- max(min(cohort.t.dlt[j-1],DLT_window + 0.001) - cohort.follow[j-1],0)
            
            ## One dose level down rule
            ### When there is pending data, do not stagger, treat the patient one dose level down. 
            ### When there is no lower dose, the patient will need to be staggered. 
            
            if(current_dose == 1 | (current_dose > 1 & pt.max.stag.time == 0)){
              ## When there is no lower dose, the patient will need to be staggered.
              ## decision time 
              t.decision = t.decision_stagger + pt.max.stag.time
              ## decide the dose assigned to this patient 
              followup.time = t.decision - t.enter
              cohort.follow <- t.decision - cohort.pts.entry
              
              
              # delta.tox <- as.numeric((cohort.pts.entry + cohort.t.dlt) <= t.decision)
              # delta.eff <- as.numeric((cohort.pts.entry + cohort.t.eff) <= t.decision)
              # 
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
              
              ### change to depend only on observed data, nopending DLT data 
              # actual_dose <- NULL
              id.complete.data <- which(omega.tox == 1)
              dose_each_patient_comp <- dose_each_patient[id.complete.data]
              n.low.tox.comp <- n.low.tox[id.complete.data]
              if(length(id.complete.data) == 1){
                cat_obs_data_comp <- matrix(cat_obs_data[id.complete.data,],ncol=4)
              }else{cat_obs_data_comp <-cat_obs_data[id.complete.data,]}
              omega.tox.comp <- omega.tox[id.complete.data]
              omega.eff.comp <- omega.eff[id.complete.data]
              eta_elicited_piE.comp <- eta_elicited_piE[id.complete.data]
              
              
              ### Run the Bayesian model to get the 
              jags.params<-c( "rho","mu_T", "mu_E")#  parameters of interest
              jags.data<-list(dose = dose, ndose=ndose,dose_each_patient=dose_each_patient_comp,  ## dose info
                              n.low.tox=n.low.tox.comp,       ## low grade toxicity info
                              cat_obs_data=cat_obs_data_comp, omega.tox=omega.tox.comp, 
                              omega.eff=omega.eff.comp,  eta_elicited_piE = eta_elicited_piE.comp,     ## tox and eff info
                              beta0_bar=beta0_bar, tau_beta0=tau_beta0, beta1_bar=beta1_bar, tau_beta1=tau_beta1,  ## prior parameter for mu_T
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
                                               plots=FALSE,silent.jags = T)}, error=function(e){NULL;cat("error-no-one-dose-down")})
              
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
              post_mean_eff <-  pnorm(post_mean_muE)
              post_mean_tox <- pnorm(post_mean_muT)
              
              ### subtract the posterior mean for pi11, pi01,p10,p00 given stagger time 
              # t_delay <- 10
              # post_mean_pi_delay <- delayed_joint_prob(mcmc.sample, 5, 30)
              # 
              
              ## first calculate the utility of all doses 
              post_utility <-  apply(post_mean_pi,1, function(x){
                sum(x*utility.table)
              } )
              
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
              
              Prob_piT_leq_bar_current <- Prob_piT_leq_bar[current_dose]
              if(Prob_piT_leq_bar_current > cut.tox.stagger){
                assign.dose = current_dose
              }else{
                Admiss_ind <- Prob_piT_leq_bar>cut.tox & Prob_piE_geq_bar>cut.eff
                Admiss_set <- which(Admiss_ind == TRUE)
                
                # Admiss_set<- c(2,3,4)
                Candidate_set <- Admiss_set[Admiss_set <= current_dose]
                
                
                ### Treat at least one cohort before terminate the trial
                if(i==1){
                  assign.dose <- current_dose
                }else if(length(Candidate_set)==0){
                  early.stop = 1
                  staggertime = c(staggertime, pt.max.stag.time )
                  print("No admissible dose and the trial is stopped.")
                  break
                }else{
                  ## adaptive randomization to lower admissible dose 
                  de_sub_AR_prob <- post_utility[Candidate_set]/sum(post_utility[Candidate_set])
                  indicator_dose <- rmultinom(1, 1, de_sub_AR_prob)
                  assign.dose <- Candidate_set[which(indicator_dose[,1]==1)]
                }
                
              }  
              
              
              staggertime = c(staggertime, pt.max.stag.time )
              cohort.t.stagger = c(cohort.t.stagger,pt.max.stag.time)
              
              t.enter = c(t.enter, t.decision)
              
              t.pt.enter <- t.decision
              
              cohort.pts.entry = c(cohort.pts.entry,t.decision) 
              
              
              mu_E_delay <- qnorm(eta_elicited_func(eta_elicited_rate, pt.max.stag.time)*pnorm(mu_E))
              
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
              cohort.dose.info <- c(cohort.dose.info, assign.dose)
              npat_initial_dose = c(npat_initial_dose, current_dose)
              idx_stagger = c(idx_stagger,ifelse(pt.max.stag.time==0, 0,1))
              
              t.low.tox <- c(t.low.tox, pts.dlt.time)
              n.low.tox <- c(n.low.tox, rpois(1, pts.dlt.time * lambda_d[assign.dose]))
            }else if(current_dose > 1 & pt.max.stag.time > 0){
              ## When there is pending data, do not stagger, treat the patient one dose level down.
              ## decision time 
              t.decision = t.decision_stagger + pt.max.stag.time
              ## decide the dose assigned to this patient 
              followup.time = t.decision - t.enter
              cohort.follow <- t.decision - cohort.pts.entry
              
              ## Backfill cohort. 
              if(arrive.dist=="Uniform"){t.seq.pt.enter = 30*runif(max_acc, 0, 2/accrual.rate)}
              if(arrive.dist=="Exponential"){t.seq.pt.enter = 30*rexp(max_acc, rate=accrual.rate)}
              
              actual.arrive.backfill <- cumsum(t.seq.pt.enter) + t.decision_stagger
              
              ### patients come during the interval 
              backfill.cohort.pts.entry <- c(t.decision_stagger, actual.arrive.backfill[which(actual.arrive.backfill < t.decision)])
              backfill.cohort.dose.info <- rep(current_dose-1,length(backfill.cohort.pts.entry))
              
              if(length(yT.true) + length(backfill.cohort.pts.entry) > max_n){
                add_n <- max_n - length(yT.true)
                backfill.cohort.pts.entry <- backfill.cohort.pts.entry[1:add_n]
                backfill.cohort.dose.info <- backfill.cohort.dose.info[1:add_n]
              }
              
              
              if(length(backfill.cohort.pts.entry)>0){
                ##Need to get this patient's outcome 
                t.enter = c(t.enter, backfill.cohort.pts.entry)
                backfill.pts.toxeff.out <- get.pts.tox.eff.outcome(d = current_dose-1, npts=length(backfill.cohort.pts.entry), mu_E, mu_T, Sigma, DLT_window, efficacy_window, dist=event.dist)
                yT.true =  c(yT.true, backfill.pts.toxeff.out$tox)
                t.tox.event =c(t.tox.event,backfill.pts.toxeff.out$t.tox)
                yE.true =  c(yE.true, backfill.pts.toxeff.out$eff)
                t.eff.event =  c(t.eff.event,backfill.pts.toxeff.out$t.eff)
                idx_stagger = c(idx_stagger, rep(0,length(backfill.cohort.pts.entry)))
                npat_dose[current_dose-1]<-npat_dose[current_dose-1] +  length(backfill.cohort.pts.entry)
                dose_each_patient = c(dose_each_patient, backfill.cohort.dose.info)
                npat_initial_dose = c(npat_initial_dose, rep(NA, length(backfill.cohort.dose.info)))
                
                ##generate the low grade tox
                n.low.tox <- c(n.low.tox, rpois(length(backfill.cohort.pts.entry),backfill.pts.toxeff.out$t.tox*lambda_d[current_dose-1]))
                t.low.tox <- c(t.low.tox, backfill.pts.toxeff.out$t.tox)
                low.tox.ind <- c(low.tox.ind, rep(1,length(backfill.cohort.pts.entry)))
                staggertime <- c(staggertime, rep(0, length(backfill.cohort.pts.entry)))
              }
              
              if(length(yT.true) == max_n){break}
              
              ### arrival of new pt.
              if(arrive.dist=="Uniform"){t.next.pt.cohort.enter = 30*runif(1, 0, 2/accrual.rate)}
              if(arrive.dist=="Exponential"){t.next.pt.cohort.enter = 30*rexp(1, rate=accrual.rate)}
              t.decision = t.next.pt.cohort.enter + t.decision
              
              followup.time = t.decision - t.enter
              ## evaluate whether the dose is admissible
              omega.tox <-NULL
              omega.eff <- NULL
              yT.observe <- NULL
              yE.observe <- NULL
              cat_obs_data <- matrix(0,nrow = length(yT.true),ncol = 4)
              eta_elicited_piE <- NULL
              
              
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
              
              ### change to depend only on observed data, nopending DLT data 
              # actual_dose <- NULL
              id.complete.data <- which(omega.tox == 1)
              dose_each_patient_comp <- dose_each_patient[id.complete.data]
              n.low.tox.comp <- n.low.tox[id.complete.data]
              if(length(id.complete.data) == 1){
                cat_obs_data_comp <- matrix(cat_obs_data[id.complete.data,],ncol=4)
              }else{cat_obs_data_comp <-cat_obs_data[id.complete.data,]}
              omega.tox.comp <- omega.tox[id.complete.data]
              omega.eff.comp <- omega.eff[id.complete.data]
              eta_elicited_piE.comp <- eta_elicited_piE[id.complete.data]
              
              ### Run the Bayesian model to get the 
              jags.params<-c( "rho","mu_T", "mu_E")#  parameters of interest
              jags.data<-list(dose = dose, ndose=ndose,dose_each_patient=dose_each_patient_comp,  ## dose info
                              n.low.tox=n.low.tox.comp,       ## low grade toxicity info
                              cat_obs_data=cat_obs_data_comp, omega.tox=omega.tox.comp, 
                              omega.eff=omega.eff.comp,  eta_elicited_piE = eta_elicited_piE.comp,     ## tox and eff info
                              beta0_bar=beta0_bar, tau_beta0=tau_beta0, beta1_bar=beta1_bar, tau_beta1=tau_beta1,  ## prior parameter for mu_T
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
                                               plots=FALSE,silent.jags = T)}, error=function(e){NULL;cat("error-one-dose-down")})
              
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
              post_mean_eff <-  pnorm(post_mean_muE)
              post_mean_tox <- pnorm(post_mean_muT)
              
              ### subtract the posterior mean for pi11, pi01,p10,p00 given stagger time 
              # t_delay <- 10
              # post_mean_pi_delay <- delayed_joint_prob(mcmc.sample, 5, 30)
              # 
              
              ## first calculate the utility of all doses 
              post_utility <-  apply(post_mean_pi,1, function(x){
                sum(x*utility.table)
              } )
              
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
              
              Prob_piT_leq_bar_current <- Prob_piT_leq_bar[current_dose]
              if(Prob_piT_leq_bar_current > cut.tox.stagger){
                assign.dose = current_dose
              }else{
                Admiss_ind <- Prob_piT_leq_bar>cut.tox & Prob_piE_geq_bar>cut.eff
                Admiss_set <- which(Admiss_ind == TRUE)
                
                # Admiss_set<- c(2,3,4)
                Candidate_set <- Admiss_set[Admiss_set <= current_dose]
                
                
                ### Treat at least one cohort before terminate the trial
                if(i==1){
                  assign.dose <- current_dose
                }else if(length(Candidate_set)==0){
                  early.stop = 1
                  staggertime = c(staggertime, pt.max.stag.time )
                  print("No admissible dose and the trial is stopped.")
                  break
                }else{
                  ## adaptive randomization to lower admissible dose 
                  de_sub_AR_prob <- post_utility[Candidate_set]/sum(post_utility[Candidate_set])
                  indicator_dose <- rmultinom(1, 1, de_sub_AR_prob)
                  assign.dose <- Candidate_set[which(indicator_dose[,1]==1)]
                }
                
              }  
              
             
              ##making staggering decision and the decision will be reapplied once DLT data is updated in this cohort 
      
              # pt.max.stag.time = 0
              staggertime = c(staggertime, 0 )
              cohort.t.stagger = c(cohort.t.stagger,0)
              
              t.enter = c(t.enter, t.decision)
              
              t.pt.enter <- t.decision
              
              cohort.pts.entry = c(cohort.pts.entry,t.decision) 
              
              mu_E_delay <- qnorm(eta_elicited_func(eta_elicited_rate, 0)*pnorm(mu_E))
              
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
              cohort.dose.info <- c(cohort.dose.info, assign.dose)
              npat_initial_dose = c(npat_initial_dose, current_dose)
              idx_stagger = c(idx_stagger,ifelse(pt.max.stag.time==0, 0,1))
              
              t.low.tox <- c(t.low.tox, pts.dlt.time)
              n.low.tox <- c(n.low.tox, rpois(1, pts.dlt.time * lambda_d[assign.dose]))
              
              ###   
              
            }
            
          } ##### J=2/3 in first cohort assign to a dose 
        }  
        
        
        if(model.unsucc.idx == 1){break}
        ### between cohort always stagger
        ### between cohort always stagger 
        if(early.stop == 1){
          t.decision = t.decision
        }else{
          t.decision = t.decision + pts.dlt.time 
        }
        
        followup.time <- t.decision - t.enter
        ###### if there is no admissible dose in the lowest dose, then terminate the trial 
        if(early.stop == 1){
          next.dose = 0
          print("No admissible dose and the trial is stopped.")
          break
        }
        # if(length(yT.true) == max_n){break}
      }else{
        ## if this dose has been tried, then, we assign all patients in this cohort to this dose
        
        ## vectors to record information within each cohort 
        cohort.yT <- NULL
        cohort.t.dlt <- NULL
        cohort.yE <- NULL
        cohort.t.eff <- NULL
        cohort.pts.entry <- NULL
        cohort.t.low.tox <- NULL
        
        
        for(j in 1:cohortsize_use){
          if(j==1) { cohort.pts.entry = c(cohort.pts.entry, t.decision); } 
          else {
            if(arrive.dist=="Uniform"){ cohort.pts.entry = c(cohort.pts.entry, cohort.pts.entry[length(cohort.pts.entry)] + 30*runif(1, 0, 2/accrual.rate))}
            if(arrive.dist=="Exponential"){ cohort.pts.entry = c(cohort.pts.entry, cohort.pts.entry[length(cohort.pts.entry)] + 30*rexp(1, rate=accrual.rate))}
          }
        }
        t.enter = c(t.enter, cohort.pts.entry)
        staggertime = c(staggertime, rep(0, cohortsize_use))
        
        ##Need to get this patient's outcome 
        pts.toxeff.out <- get.pts.tox.eff.outcome(d = current_dose, npts=cohortsize_use, mu_E, mu_T, Sigma, DLT_window, efficacy_window, dist=event.dist)
        
        ## Record patients tox and eff outcomes 
        yT.true = c(yT.true, pts.toxeff.out$tox)
        t.tox.event=c(t.tox.event, pts.toxeff.out$t.tox)
        
        yE.true = c(yE.true, pts.toxeff.out$eff)
        t.eff.event = c(t.eff.event, pts.toxeff.out$t.eff)
        
        
        npat_dose[current_dose]<-npat_dose[current_dose] + cohortsize_use
        
        dose_each_patient = c(dose_each_patient, rep(current_dose, cohortsize_use))
        npat_initial_dose = c(npat_initial_dose, rep(current_dose, cohortsize_use))
        
        idx_stagger <- c(idx_stagger, rep(0, cohortsize_use))
        
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
        n.low.tox <- c(n.low.tox, rpois(cohortsize_use,t.tox.event[(length(dose_each_patient) - cohortsize_use +1):length(dose_each_patient)]*lambda_d[current_dose]))
        t.low.tox <- c(t.low.tox, t.tox.event[(length(dose_each_patient) - cohortsize_use +1):length(dose_each_patient)])
        low.tox.ind <- c(low.tox.ind, rep(1,cohortsize_use))
      }
      
      
      ### make the escalation decision 
  
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
      
      ### change to depend only on observed data, nopending DLT data 
      # actual_dose <- NULL
      id.complete.data <- which(omega.tox == 1)
      dose_each_patient_comp <- dose_each_patient[id.complete.data]
      n.low.tox.comp <- n.low.tox[id.complete.data]
      if(length(id.complete.data) == 1){
        cat_obs_data_comp <- matrix(cat_obs_data[id.complete.data,],ncol=4)
      }else{cat_obs_data_comp <-cat_obs_data[id.complete.data,]}
      omega.tox.comp <- omega.tox[id.complete.data]
      omega.eff.comp <- omega.eff[id.complete.data]
      eta_elicited_piE.comp <- eta_elicited_piE[id.complete.data]
      
      
      ### Run the Bayesian model to get the 
      jags.params<-c("rho","mu_T", "mu_E")#  parameters of interest
      jags.data<-list(dose = dose, ndose=ndose,dose_each_patient=dose_each_patient_comp,  ## dose info
                      n.low.tox=n.low.tox.comp,       ## low grade toxicity info
                      cat_obs_data=cat_obs_data_comp, omega.tox=omega.tox.comp, 
                      omega.eff=omega.eff.comp,  eta_elicited_piE = eta_elicited_piE.comp,     ## tox and eff info
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
      runjagsmodel<-tryCatch({run.jags(model =eff_tox_model, monitor = jags.params,data = jags.data,
                                       n.chains = 4,adapt = 2000,burnin = 5000,inits = jags.inits,
                                       sample = 10000,summarise = FALSE,thin = 1,method = 'rjparallel',
                                       plots=FALSE,silent.jags = T)}, error=function(e){NULL;cat("error-dose-es")})
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
      post_mean_eff <-  pnorm(post_mean_muE)
      post_mean_tox <- pnorm(post_mean_muT)
      
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
      
      if(length(yT.true) == max_n){break} ## if maximum sample size is reached. stop enroll patient, select OBD
    }  ## loop for cohorts 
    
    
    if(model.unsucc.idx == 0){
      ### dose selection 
      sim.data <- data.frame(dose_each_patient=dose_each_patient, yT.true = yT.true, yE.true = yE.true, 
                             n.low.tox=n.low.tox, t.low.tox=t.low.tox,idx_stagger = idx_stagger,npat_initial_dose=npat_initial_dose)
      
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
      ## Number of patients staggerred 
      nstagger_treated <- sum(staggertime[1:length(yT.true)]>0)
      
      ntotal_stagger_treated <- sum(idx_stagger)
      
      ntotal_pts_treated <- length(yT.true)
      
      
      ##trial length
      trial.length <- t.decision
     
      
      report_result <- rbind(report_result, c(OBD, trial.length, ntotal_pts_treated,nstagger_treated, ntotal_stagger_treated,
                                               npat_dose))
      
      result_all[[trial]] <- list(sim.data, sim.dose.info, report_result,staggertime=staggertime)
      
    }  
    #fail_mark <- c(fail_mark, model.unsucc.idx)
    print(paste("Simulation time: ", trial, sep=""))
    # trial = trial +1
  }
  
  s2 <- Sys.time()
  s2-s1
  report_result <- data.frame(report_result)
  colnames(report_result) <- c("OBD", "Trial Duration", "# pts treated", "# Stagger" ,"# stagger decision",
                              paste("# pts in Dose ", seq(1:ndose), sep = "") )
  
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
  
  pct_stagger_treated <- mean(report_result[,4])/mean(report_result$`# stagger decision`)*100

  
  pct_pts_allo <- apply(report_result[,-c(1,2,3,4,5)],2,mean)/mean(report_result$`# pts treated`)*100
  pct_pts_allo <- matrix(pct_pts_allo, byrow = TRUE, nrow=1)
  colnames(pct_pts_allo) <- c(paste("% pts in Dose", 1:ndose,sep = ""))
  print(pct_pts_allo)
  
  #result.out <- cbind(percent_select,Average_result, pct_pts_allo)
  result.out <- cbind(scenario = sce, subset=subset, percent_select,Average_result, `% stagger` = pct_stagger_treated,
                      pct_pts_allo, actual.Nsim = Nsim - model.unsucc)

  
  f1 <- sprintf("Results/%s_ctox_%.2f_ceff_%.2f_ces_%.2f_cstag_%.2f_scen_%d_subset_%d.csv", File_name, cut.tox,cut.eff, cut.escalation,cut.tox.stagger,sce,subset)
  
  write.csv(result.out,file=f1)
  
  return(result_all)
}

#####################################
# ---------------------------------
# ####      Example           ####
#-----------------------------------
#  dose <- c( 4e+7, 8e+7, 4e+8, 8e+8,4e+9)
# 
#  ## True tox, eff,low-grade toxicity
#  pT.true = c(0.1, 0.25, 0.5, 0.65, 0.7)
#  pE.true =c(0.05, 0.1, 0.15, 0.3, 0.45)
#  lambda_d = c(10,20,40,48,60)/30
# 
#  ## start dose
#  startdose <- 1
#  ## maximum sample size
#  max_n <- 45
# 
#  ## cohort size
#  cohortsize=3
# 
# 
#  ## toxicity upper bound
#  pi_T_bar <- 0.35
# 
#  ## efficacy lower bound
#  pi_E_bar <- 0.25
# 
#  ## toxicity/efficacy Assessment window
#  DLT_window <- 30
#  efficacy_window <- 30
# 
#  ## required staggering window
#  stagger_window <- 30
# 
#  ## patients accrual rate
#  accrual.rate <- 5   ## average 5 patients per month
#  ## distribution used to generate arrival time: accept "Uniform" and "Exponential"
#  arrive.dist =  "Uniform"
#  # distribution used to generate time to event: accept "Uniform" and "Weibull"
#  event.dist = "Uniform"
# 
# 
#  ### utility
#  utility.table <- c(60,100, 0, 40)  ## U(1,1), U(0,1), U(1,0), U(0,0)
# 
#  ##Parameters for prior distribution
#  beta0_bar =-1.2
#  tau_beta0 = 1/(4*beta0_bar)^2
#  beta1_bar = 0.69
#  tau_beta1 = 1/(4*beta1_bar)^2  ## prior parameter for mu_T
#  alpha0_bar = -1.64
#  tau_alpha0 = 1/(4*alpha0_bar)^2
#  alpha1_bar = 0.2
#  tau_alpha1 = 1/(4*alpha1_bar)^2
#  alpha2_bar = 3.6
#  tau_alpha2 = 1/(4*alpha2_bar)^2
#  alpha3_bar = -1
#  tau_alpha3 = 1/(4*alpha3_bar)^2
# 
# 
#  ## correlation between two latent outcomes
#  rho=0.3
#  Sigma <- matrix(c(1,rho,rho,1), ncol=2)   #covariance matrix of the variables.
# 
#  ## Number of simulation
#  Nsim = 5
# 
# 
#  ### response rate decay elicited from physicians unit: weeks
#  eta_elicited_rate <- c(0.05,0.1,0.15,0.2)
# 
#  ## admissible cutoffs
#  cut.tox <- c(0.1)
#  cut.eff <- c(0.05)
#  cut.escalation <- c(0.5)
#  cut.tox.stagger <- c(0.4)
# 
# 
# s1 <- Sys.time()
# Backfill_design_sce <- Backfill_design(dose, startdose, max_n, cohortsize,
#                                               DLT_window, efficacy_window,  stagger_window,
#                                               pi_T_bar, pi_E_bar, eta_elicited_rate, utility.table,
#                                               accrual.rate, arrive.dist,event.dist,
#                                               lambda_d, rho, pT.true, pE.true,
#                                               cut.tox,cut.eff,cut.tox.stagger, cut.escalation,Nsim,
#                                               beta0_bar, tau_beta0, beta1_bar, tau_beta1,  ## prior parameter for mu_T
#                                               alpha0_bar,tau_alpha0,alpha1_bar,tau_alpha1,alpha2_bar,tau_alpha2, ## prior parameter for mu_E
#                                               alpha3_bar, tau_alpha3)
# Backfill_design_sce


