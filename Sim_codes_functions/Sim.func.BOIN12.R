#############################################################
#  Function to generate operating characteristics of BOIN12  #
#############################################################

source("Sim_codes_functions/Get.patients.tox.eff.outcome.R")
source("Sim_codes_functions/Response.decay.rate.func.R")
source("Sim_codes_functions/Functions.R")

library(MASS)
library(matrixcalc) ## check whether a matrix is positive definite
library(mvtnorm)

#' @param pT.true the true marginal toxicity rate vector
#' @param pE.true the true marginal efficacy rate vector
#' @param max_n maximum sample size
#' @param cohortsize the cohort size
#' @param startdose the starting dose level for the trial
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
#' @param p.tox the lowest DLT probability that is deemed overly toxic such
#'              that deescalation is required. The default value is \code{p.tox=1.4*target}).
#' @param N1 if number of patients at the current dose level >= N1, then the next higher dose
#'           level won't be considered when the toxicity rate at the current dose is greater than escalation boundary
#' @param N2 if the current dose has number of patients  >= N2 and is safe, and the next higher dose level has not been
#'           tested, then we will assign the next cohort to the next cohort.
#' @param cut.tox the cutoff to eliminate an overly toxic dose for safety.
#'                  We recommend the default value of (\code{cut.tox=0.9}) for general use.
#'
#' @param cut.eff the cutoff to eliminate a futile dose.
#'                  We recommend the default value of (\code{cut.eff=0.95}) for general use.
#' @param Nsim the number of simulated trials
#' @param sce indicator of scenario, in our simulations, we considered six scenarios



BOIN12_Design <- function(pT.true, pE.true, max_n,cohortsize,startdose,
                          DLT_window, efficacy_window,stagger_window,
                          pi_T_bar,pi_E_bar,utility.table, eta_elicited_rate,
                          accrual.rate, arrive.dist, event.dist,
                          lambda_d, rho,p.saf,p.tox,cut.tox,cut.eff,
                           N1, N2,Nsim,sce=1){
  #set.seed(rseed)
  
  s1<- Sys.time()
  
  ## number of doses 
  ndose<-length(pT.true)
  #number of cohorts 
  ncohort = max_n/cohortsize
  
  mu_E =   qnorm(pE.true)
  mu_T =   qnorm(pT.true)
  
  Sigma <- matrix(c(1,rho,rho,1), ncol=2)   #covariance matrix of the variables.
  
  report_result <- c()
  result_all <- list()
  
  
  temp= get.boundary.BOIN12(pi_T_bar, pi_E_bar, max_n, 1,cutoff.eli=cut.tox, cutoff.eli.E = cut.eff)
  p.tox=1.4*pi_T_bar
  b.e=temp[4,];   # escalation boundary
  b.d=temp[3,];   # deescalation boundary
  
  b.elim=temp[2,];  # elimination boundary
  b.elimE=temp[5,]
  u01=utility.table[2]
  u10=utility.table[3]
  u11 = utility.table[1]
  u00 = utility.table[4]
  # utility.table <- c(60,100, 0, 40)  ## U(1,1), U(0,1), U(1,0), U(0,0)
  
  utility=c(u11,u10,u01,u00)
  # Assume independence between toxicity and efficacy
  targetP<-c(pi_E_bar*pi_T_bar,pi_T_bar*(1-pi_E_bar),(1-pi_T_bar)*pi_E_bar,(1-pi_T_bar)*(1-pi_E_bar))
  
  # Calculate the benchmark utility
  uu = sum(targetP*utility) #highest unacceptable utility
  uu = uu+(100-uu)/2        # benchmark utility (i.e., desirable utility)
  
  
  for(trial in 1:Nsim){
    
    set.seed(2024 + trial)
    
    ##Count the number of stagger
    nstagger <- 0
    nstagger_dose <- rep(0,ndose)
    staggertime <-NULL
    
    ## patients enrolled in each dose
    npat_dose <- rep(0, ndose) 
    
    yT.true=yE.true=NULL;  #toxicity/efficacy indicator for each subject
    X.true = NULL;         #number of observed low-grade toxicity for each subject
    
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
    npat_initial_dose <- NULL
    
    early.stop=0;              ## indicate whether the trial terminates early
    d=startdose;               ## starting dose level
    elimi = rep(0, ndose);     ## indicate whether doses are eliminated due to toxicity
    elimiE=  rep(0,ndose);     ## indicate whether doses are eliminated due to efficacy
    safe = 0
    posH<-rep(1-uu/100,ndose)
    
    
    yT <- yE <- rep(0, ndose);    ## number of DLT/efficacy at each dose level
    y01 <- y10 <- y11 <- y00 <-rep(0,ndose); ## number of different outcomes at each dose level
    
    
    
    for(i in 1:ncohort){
      if(i == 1){
        current_dose = startdose
      }else{
        current_dose = next.dose
      }
      if(npat_dose[current_dose] == 0){
        
        ## vectors to record information within each cohort 
        cohort.yT <- NULL
        cohort.t.dlt <- NULL
        cohort.yE <- NULL
        cohort.t.eff <- NULL
        cohort.pts.entry <- NULL
        cohort.t.stagger <- NULL
        n.cohort.low.tox <- NULL
        cohort.t.low.tox <- NULL
        cohort.dose.info <- NULL
        
        for(j in 1:cohortsize){
          if(j == 1){
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
            
            y01[current_dose]<-y01[current_dose]+sum(pts.toxeff.out$tox == 0 &pts.toxeff.out$eff == 1)
            y10[current_dose]<-y10[current_dose]+sum(pts.toxeff.out$tox == 1 &pts.toxeff.out$eff == 0)
            y11[current_dose]<-y11[current_dose]+sum(pts.toxeff.out$tox == 1 &pts.toxeff.out$eff == 1)
            y00[current_dose]<-y00[current_dose]+sum(pts.toxeff.out$tox == 0 &pts.toxeff.out$eff == 0)
            
            yT<-y10+y11
            yE<-y01+y11
            
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
            idx_stagger = c(idx_stagger,0)
            
          }else{
            if(arrive.dist=="Uniform"){t.next.pt.enter = 30*runif(1, 0, 2/accrual.rate)}
            if(arrive.dist=="Exponential"){t.next.pt.enter = 30*rexp(1, rate=accrual.rate)}
            ##making staggering decision and the decision will be reapplied once DLT data is updated in this cohort 
            
            # when the new patient arrive, we need to make the decision of whether to stagger
            t.decision_stagger = t.pt.enter + t.next.pt.enter
            
            #follow up time for enrolled patients 
            followup.time <- t.decision_stagger - t.enter
            
            ## follow up time for patients in this cohort
            cohort.follow <- followup.time[((i-1)*cohortsize+1):length(followup.time)]
            
            ## dose assigned to patients in this cohort
            cohort.dose.info <- dose_each_patient[((i-1)*cohortsize+1):length(yT.true)]
            
            ## stagger time for this patient 
            pt.max.stag.time <- max(min(cohort.t.dlt[j-1],DLT_window + 0.001) - cohort.follow[j-1],0)
            
            ## decision time 
            t.decision = t.decision_stagger + pt.max.stag.time
            
            ## decide the dose assigned to this patient 
            followup.time = t.decision - t.enter
            
            cohort.follow <- followup.time[((i-1)*cohortsize+1):length(followup.time)]
            
            ## After waiting, evaluate whether the dose is admissible and design which dose should be assigned to this patient
            assign.dose <- current_dose
            
            within_nc <- npat_dose[assign.dose]/1
            
            # determine whether current dose level is overly toxic
            if(!is.na(b.elim[within_nc]))
            {
              if(yT[current_dose]>=b.elim[within_nc])
              {
                elimi[current_dose:ndose]=1;
                if(current_dose==1) {
                  early.stop=1
                  assign.dose = 0
                  print("No admissible dose and the trial is stopped.")
                  break
                }
              }
            }
            
            if(!is.na(b.elimE[within_nc]))
            {
              if(yE[current_dose]<=b.elimE[within_nc])
              {
                elimi[current_dose]=1;
              }
            }
            if(sum(elimi==1)==ndose) {
              early.stop=1
              assign.dose = 0
              print("No admissible dose and the trial is stopped.")
              break
            }
            
            ### if current dose is safe, then stay, if not safe, determine using utility 
            ptox.est.current <- yT[current_dose]/npat_dose[current_dose]
            if(ptox.est.current < p.tox){
              assign.dose= current_dose
            }else{
              u_curr<-(u01*y01[current_dose]+u10*y10[current_dose]+u11*y11[current_dose]+u00*y00[current_dose])/100
              posH[current_dose] = 1-pbeta(uu/100,1+u_curr,npat_dose[current_dose]-u_curr+1)
              posH <- posH*(1-elimi);
              if(npat_dose[current_dose]>=N1){safe=1} else{safe=0}
              if(npat_dose[current_dose]>=max_n){break}
              ## choose the dose 
              ## find admissible set first 
              Candidate_dose <- 1:current_dose
              Candidate_set <- Candidate_dose[which(elimi[Candidate_dose]==0)]
              if(i == 1){
                assign.dose <- current_dose
              }else if(length(Candidate_set)==0){
                early.stop = 1
                staggertime = c(staggertime, pt.max.stag.time )
                print("No admissible dose and the trial is stopped.")
              }else{
                de_sub_AR_prob <- posH[Candidate_set]/sum(posH[Candidate_set])
                indicator_dose <- rmultinom(1, 1, de_sub_AR_prob)
                assign.dose <- Candidate_set[which(indicator_dose[,1]==1)]
              }
              
            }
            
            ###
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
            
            y01[assign.dose]<-y01[assign.dose]+sum(pts.toxeff.out$tox == 0 &pts.toxeff.out$eff == 1)
            y10[assign.dose]<-y10[assign.dose]+sum(pts.toxeff.out$tox == 1 &pts.toxeff.out$eff == 0)
            y11[assign.dose]<-y11[assign.dose]+sum(pts.toxeff.out$tox == 1 &pts.toxeff.out$eff == 1)
            y00[assign.dose]<-y00[assign.dose]+sum(pts.toxeff.out$tox == 0 &pts.toxeff.out$eff == 0)
            
            
            yT<-y10+y11
            yE<-y01+y11
            
            cohort.yT <- c(cohort.yT, pts.dlt)
            cohort.t.dlt <- c(cohort.t.dlt, pts.dlt.time)
            cohort.yE <- c(cohort.yE, pts.eff)
            cohort.t.eff <- c(cohort.t.eff, pts.eff.time)
            
            npat_dose[assign.dose]<-npat_dose[assign.dose] + 1
            
            dose_each_patient = c(dose_each_patient, assign.dose)
            npat_initial_dose = c(npat_initial_dose, current_dose)
            idx_stagger = c(idx_stagger,ifelse(pt.max.stag.time==0, 0,1))
            
            t.low.tox <- c(t.low.tox, pts.dlt.time)
            n.low.tox <- c(n.low.tox, rpois(1, pts.dlt.time * lambda_d[assign.dose]))
            ###
            
          }
        }### j
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
        
      }else{
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
        
        t.decision = max (cohort.pts.entry + pts.toxeff.out$t.tox)
        
        
        #follow up time for enrolled patients 
        followup.time <- t.decision - t.enter
        
        ### # of stagger decisions 
        idx_stagger <- c(idx_stagger, rep(0, cohortsize))
        
        y01[current_dose]<-y01[current_dose]+sum(pts.toxeff.out$tox == 0 &pts.toxeff.out$eff == 1)
        y10[current_dose]<-y10[current_dose]+sum(pts.toxeff.out$tox == 1 &pts.toxeff.out$eff == 0)
        y11[current_dose]<-y11[current_dose]+sum(pts.toxeff.out$tox == 1 &pts.toxeff.out$eff == 1)
        y00[current_dose]<-y00[current_dose]+sum(pts.toxeff.out$tox == 0 &pts.toxeff.out$eff == 0)
        npat_dose[current_dose]<-npat_dose[current_dose]+cohortsize
        
        
        yT<-y10+y11
        yE<-y01+y11
        
        dose_each_patient = c(dose_each_patient, rep(current_dose, cohortsize))
        
        n.low.tox <- c(n.low.tox, rpois(cohortsize,t.tox.event[(length(dose_each_patient) - cohortsize +1):length(dose_each_patient)]*lambda_d[current_dose]))
        t.low.tox <- c(t.low.tox, t.tox.event[(length(dose_each_patient) - cohortsize +1):length(dose_each_patient)])
        
      }
      if(length(yT.true) == max_n){
        t.decision =  max(t.enter + pmax(t.tox.event, t.eff.event))
        followup.time = t.decision - t.enter
      }else{
        t.decision = t.decision
        followup.time = t.decision - t.enter
      }
      
      nc<-npat_dose[current_dose]/1   #nc-number of cohorts treated
      
      ## we need to wait until all patients' results are observed then make the escalation decision. 
      
      next.dose <- current_dose
      
      # determine whether current dose level is overly toxic
      if(!is.na(b.elim[nc]))
      {
        if(yT[current_dose]>=b.elim[nc])
        {
          elimi[current_dose:ndose]=1;
          if(current_dose==1) {
            early.stop=1
            next.dose = 0
            print("No admissible dose and the trial is stopped.")
            break
          }
        }
      }
      
      if(!is.na(b.elimE[nc]))
      {
        if(yE[current_dose]<=b.elimE[nc])
        {
          elimi[current_dose]=1;
        }
      }
      if(sum(elimi==1)==ndose) {
        early.stop=1
        next.dose = 0
        print("No admissible dose and the trial is stopped.")
        break
      }
      
      u_curr<-(u01*y01[current_dose]+u10*y10[current_dose]+u11*y11[current_dose]+u00*y00[current_dose])/100
      posH[current_dose] = 1-pbeta(uu/100,1+u_curr,npat_dose[current_dose]-u_curr+1)
      posH <- posH*(1-elimi);
      if(npat_dose[current_dose]>=N1){safe=1} else{safe=0}
      if(npat_dose[current_dose]>=max_n){break}
      
      if (yT[current_dose]>=b.d[nc] && current_dose!=1) {
        if(sum(elimi[1:(current_dose-1)]==0)>0){next.dose=max(which(elimi[1:(current_dose-1)]==0))} else {
          if(elimi[current_dose]==1){early.stop=1;break} else{next.dose=current_dose}}
      } else if (yT[current_dose]>=b.d[nc] && current_dose==1) {if(elimi[current_dose]==0){next.dose=current_dose} else{
        early.stop=1
        next.dose = 0
        print("No admissible dose and the trial is stopped.")
        break
      }
      } else{
        admi_set=current_dose;
        if(current_dose>1){
          if(sum(elimi[1:(current_dose-1)]==0)>0){admi_set<-c(admi_set,max(which(elimi[1:(current_dose-1)]==0)))}
        }
        if(current_dose<ndose){
          if(safe==0){
            if(sum(elimi[(current_dose+1):ndose]==0)>0){admi_set<-c(admi_set,current_dose+min(which(elimi[(current_dose+1):ndose]==0)))}
          } else {
            if(yT[current_dose]<=b.e[nc] & sum(elimi[(current_dose+1):ndose]==0)>0){admi_set<-c(admi_set,current_dose+min(which(elimi[(current_dose+1):ndose]==0)))}
          }
        }
        
        temp.posH<-posH[admi_set]+runif(length(admi_set))*(10^-15)
        next.dose=admi_set[which.max(temp.posH)]
        
        #####-----------------------------
        if(length(admi_set)==3 ){
          posH_max.index = which(posH[admi_set] == max(posH[admi_set]))
          next.dose = max(admi_set[posH_max.index])
        }
      }
      
      if (elimi[next.dose]==1) {
        early.stop=1
        next.dose = 0
        print("No admissible dose and the trial is stopped.")
        break
      }
      if (sum(elimi)==ndose) {
        early.stop=1
        next.dose = 0
        print("No admissible dose and the trial is stopped.")
        break
      }
      
      if (current_dose<ndose){
        if(sum(elimi[(current_dose+1):ndose]==0)>0){
          d_temp=current_dose+min(which(elimi[(current_dose+1):ndose]==0))
          if(npat_dose[current_dose]>=N2 & npat_dose[min(d_temp,ndose)]==0 & yT[current_dose]<b.d[npat_dose[current_dose]/1]){next.dose<-d_temp}
        }
      }
      
    }   ### 
    
    sim.data <- data.frame(dose_each_patient=dose_each_patient, yT.true = yT.true, yE.true = yE.true, n.low.tox=n.low.tox, t.low.tox=t.low.tox,idx_stagger=idx_stagger)
    est_tox <- yT/npat_dose
    est_eff <- yE/npat_dose
    est.utility <- (u01*y01+u10*y10+u11*y11+u00*y00)/100
    sim.dose.info <- data.frame(est_tox=est_tox, est_eff = est_eff, post_utility=est.utility)
    
    ## OBD of this trial
    if(early.stop==1){
      OBD = 0
    }else{
      pT<-rep(pi_T_bar,ndose)
      pE<-rep(pi_E_bar,ndose)
      pT[npat_dose!=0]<-(yT[npat_dose!=0])/(npat_dose[npat_dose!=0])
      pE[npat_dose!=0]<-(yE[npat_dose!=0])/(npat_dose[npat_dose!=0])
      #pE<-peestimate(yE,n)
      pT<-pava(pT,w = 1/((yT + 0.05) * (npat_dose - yT + 0.05)/((npat_dose + 0.1)^2 * (npat_dose + 0.1 + 1))))+0.001*seq(1,ndose)
      d_mtd<-which.min(abs(pT-pi_T_bar))
      pT[npat_dose!=0]<-(yT[npat_dose!=0])/(npat_dose[npat_dose!=0])
      u<-(u01*y01+u10*y10+u11*y11+u00*y00)/100
      #u1*pE+(1-pT)*u2
      #u<-u/100*n
      u<-(u+1)/(npat_dose+2)
      u[elimi==1]<--100
      u[elimiE==1]<--100
      u[npat_dose==0]<--100
      OBD<-which.max(u[1:d_mtd])
    }
    
    
    ##trial length
    trial.length <- t.decision/30
    
    
    ## Number of patients staggerred 
    nstagger_treated <- sum(staggertime[1:length(yT.true)]>0)
    
    ntotal_stagger_treated <- sum(idx_stagger)
    
    # nstagger_enrolled <- sum(staggertime>0)
    # ntotal_stagger_enrolled <- sum(idx_stagger) + ifelse(length(staggertime) > length(yT.true) & staggertime[length(staggertime)] > 0,1,0)
    # 
    ntotal_pts_treated <- length(yT.true)
    
    # ntotal_pts_enrolled <- length(staggertime)
    
    ntotal_pts <- length(yT.true)
    
    ## number of response and no.of DLT
    No.DLT <- sum(yT) 
    No.Res <- sum(yE)
    
    report_result <- rbind(report_result, c(OBD, trial.length, ntotal_pts_treated,nstagger_treated, ntotal_stagger_treated,
                                             npat_dose, No.DLT, No.Res))
    
    
    result_all[[trial]] <- list(sim.data, sim.dose.info, report_result,staggertime=staggertime)
    
    #print(paste("Simulation time: ", trial, sep=""))
    
  }
  s2<- Sys.time()
  s2-s1
  
  report_result <- data.frame(report_result)
  colnames(report_result) <- c("OBD", "Trial Duration", "# pts treated","# Stagger" ,
                               "# stagger decision",
                               paste("# pts in Dose ", seq(1:ndose), sep = ""), "# DLT", "# Res" )
  
  percent_select <- c()
  for(d in 1:(1+ndose)){
    percent_select[d] <- sum(report_result$OBD==(d-1))/(Nsim) 
  }
  
  percent_select <- matrix(percent_select, byrow = TRUE, nrow=1)*100
  colnames(percent_select)<- paste("Dose ", seq(0:ndose)-1, sep = "") 
  #print("Percentage of selection")
  #print(percent_select)
  
  Average_result <- apply(report_result[,-1],2,mean)
  col_name <- names(Average_result)
  Average_result <- matrix(Average_result, byrow = TRUE, nrow=1)
  colnames(Average_result) <- col_name
  
  #print(Average_result)
  
  pct_stagger <- mean(report_result$`# Stagger`)/mean(report_result$`# stagger decision`) *100
 
  pct_pts_allo <- apply(report_result[,-c(1,2,3,4,5,6+ndose,7+ndose)],2,mean)/mean(report_result$`# pts treated`)*100
  pct_pts_allo <- matrix(pct_pts_allo, byrow = TRUE, nrow=1)
  colnames(pct_pts_allo) <- c(paste("% pts in Dose", 1:ndose,sep = ""))
  # print(pct_pts_allo)
  
  result.out <- cbind(scenario = sce, percent_select, Average_result, `% stagger` = pct_stagger,pct_pts_allo, actual.Nsim = Nsim)
  
  #f1 <- sprintf("BOIN12_ctox_%.2f_ceff_%.2f_scen_%d.csv", cut.tox,cut.eff, sce)
  
  # f1 <- sprintf("/rsrch5/scratch/biostatistics/swang23/Staggering/Results/BOIN12_ctox_%.2f_ceff_%.2f_scen_%d.csv", cut.tox,cut.eff, sce)
  # 
  # write.csv(result.out,file=f1)
  
  return(result.out)
}

# ---------------------------------
# ####      Example           ####
#-----------------------------------
# # ## True tox, eff,low-grade toxicity
# pT.true = c(0.1, 0.25, 0.5, 0.65, 0.7)
# pE.true =c(0.05, 0.1, 0.15, 0.3, 0.45)
# lambda_d = c(10,20,40,48,60)/30
# ## start dose
# startdose <- 1
# ## maximum sample size
# max_n <- 45
# 
# ## cohort size
# cohortsize=3
# 
# ## toxicity upper bound
# pi_T_bar <- 0.35
# 
# ## efficacy lower bound
# pi_E_bar <- 0.25
# 
# ## decision boundaries
# p.saf=0.6*pi_T_bar
# p.tox=1.4*pi_T_bar
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
# ## correlation between two latent outcomes
# rho=0.3
# Sigma <- matrix(c(1,rho,rho,1), ncol=2)   #covariance matrix of the variables.
# 
# ## Additional rules
# N1=6
# N2=9
# 
# ## Number of simulation
# Nsim = 5000
# 
# cut.tox <- c(0.9)
# cut.eff <- c(0.95)
# 
# eta_elicited_rate <- c(0.05,0.1,0.15,0.2)
# 
# 
# 
# BOIN12_sim_results <- BOIN12_Design(pT.true, pE.true, max_n,cohortsize,startdose,
#                                     DLT_window, efficacy_window,stagger_window,
#                                     pi_T_bar,pi_E_bar,utility.table, eta_elicited_rate,
#                                     accrual.rate, arrive.dist, event.dist,
#                                     lambda_d, rho,p.saf,p.tox,cut.tox,cut.eff,
#                                      N1, N2,Nsim)
# 
# BOIN12_sim_results

