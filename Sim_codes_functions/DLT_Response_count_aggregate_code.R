##########################################################################
#  Function to aggregate the simulation results to get #DLT and #Response
#########################################################################
#' @description The parameters are used to identify the files you want to aggregate
#'              If your files are named in other ways, you can change the parameters 
#'              required to identify the files accordingly
#' @param Design The prefix of the file name indicating the design you conducted
#' @param cut.tox Admissible safety cutoff used in the trial
#' @param cut.eff Admissible efficacy cutoff used in the trial
#' @param cut.escalation Escalation cutoff used in the trial   
#' @param cut.tox.stagger Safety evaluation cutoff used in the trial 
#' @param scen_num The number of scenarios considered
#' @param subset  The number of subsets running for each scenario. 
#' @param Nsim    In each subset, the number of simulations were run,
#' @param ndose   The number of doses in the trial

source("Sim_codes_functions/Sim.scenarios.R")
Sce.list <- length(Scenarios)
scen_num <- seq(from = 1, to=Sce.list,by=1)

library(plyr)


aggregate_DLT_Res_results <- function(Design,cut.tox, cut.eff, cut.escalation, cut.tox.stagger, scen_num, subset,Nsim, ndose){
  

  if(Design == "CUB"|Design == "CUB_AR3"|Design == "CUB_AR3"|Design == "CUB_N60"|
     Design == "CUB_Weibull"|Design == "CUB_U70"){
    filename <-sprintf("Results/%s_ctox_%.2f_ceff_%.2f_ces_%.2f", 
                       Design, cut.tox,cut.eff, cut.escalation)

    f1_com <- NULL
    for(i in 1:length(scen_num)){
      num_DLT <- NULL
      num_Respond <- NULL
      actual.Nsim <- 0
      for(j in 1:subset){
        oname1 <- paste(filename,"_scen_",scen_num[i],"_subset_",j, sep="")
        
        tryCatch({
          f1 <-  load(paste(oname1, ".RData", sep=""))
          #sim = dim(Hybrid_Utility_sce[[length(Hybrid_Utility_sce)]][[3]])[1]
          #actual.Nsim = actual.Nsim + sim
          
          for(c in 1:length(CUB_design_sce)){
            if(!is.null(CUB_design_sce[[c]])){
              num_DLT <- c(num_DLT, sum(CUB_design_sce[[c]][[1]]$yT.true))
              num_Respond <- c(num_Respond, sum(CUB_design_sce[[c]][[1]]$yE.true))
              actual.Nsim = actual.Nsim+1
            }
          }
          
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        #f1_com <- rbind(f1_com,f1)
      }
      
      mean_num_DLT = mean(num_DLT)
      mean_num_Respond = mean(num_Respond)
      
      f1_com <- rbind(f1_com, c(i, mean_num_DLT, mean_num_Respond, actual.Nsim))
    }
    
    
    combine.result <- data.frame(f1_com)
    colnames(combine.result) <- c("scenario", "# pts DLT", "# pts Respond", "actual.Nsim") 
    
    
  }else if(Design == "Strict_Stagger"|Design == "Strict_Stagger_N60"|Design == "Backfill"|
           Design == "Strict_Stagger_AR3"|Design == "Strict_Stagger_AR7"|Design == "Strict_Stagger_quickdecay"|Design == "Strict_Stagger_nodecay"|
           Design == "Strict_Stagger_Weib"|Design == "Strict_Stagger_U70"){
    filename <-sprintf("Results/%s_ctox_%.2f_ceff_%.2f_ces_%.2f_cstag_%.2f", 
                       Design, cut.tox,cut.eff, cut.escalation, cut.tox.stagger)
    f1_com <- NULL
    for(i in 1:length(scen_num)){
      num_DLT <- NULL
      num_Respond <- NULL
      actual.Nsim <- 0
      for(j in 1:subset){
        oname1 <- paste(filename,"_scen_",scen_num[i],"_subset_",j, sep="")
        
        tryCatch({
          f1 <-  load(paste(oname1, ".RData", sep=""))
          #sim = dim(Hybrid_Utility_sce[[length(Hybrid_Utility_sce)]][[3]])[1]
          #actual.Nsim = actual.Nsim + sim
          
          for(c in 1:length(Strict_Stagger_sce)){
            if(!is.null(Strict_Stagger_sce[[c]])){
              num_DLT <- c(num_DLT, sum(Strict_Stagger_sce[[c]][[1]]$yT.true))
              num_Respond <- c(num_Respond, sum(Strict_Stagger_sce[[c]][[1]]$yE.true))
              actual.Nsim = actual.Nsim+1
            }
          }
          
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        #f1_com <- rbind(f1_com,f1)
      }
      
      mean_num_DLT = mean(num_DLT)
      mean_num_Respond = mean(num_Respond)
      
      f1_com <- rbind(f1_com, c(i, mean_num_DLT, mean_num_Respond, actual.Nsim))
    }
    
    
    combine.result <- data.frame(f1_com)
    colnames(combine.result) <- c("scenario", "# pts DLT", "# pts Respond", "actual.Nsim") 
    
  }else if(Design == "Adaptive_Stagger"|Design == "Adaptive_Stagger_N60"|
           Design == "Adaptive_Stagger_quickdecay"|Design == "Adaptive_Stagger_AR3"|
           Design == "Adaptive_Stagger_AR7"|Design == "Adaptive_Stagger_nodecay"|
           Design == "Adaptive_Stagger_Weibull"|Design == "Adaptive_Stagger_U70"|
           Design == "Adaptive_Stagger_sens_A"|Design == "Adaptive_Stagger_sens_B"|
           Design == "Adaptive_Stagger_sens_no"|Design == "Adaptive_Stagger_Cor0.3"|
           Design == "Adaptive_Stagger_Cor0.5"|Design == "Adaptive_Stagger_Cor0.8"){
    filename <-sprintf("Results/%s_ctox_%.2f_ceff_%.2f_ces_%.2f_cstag_%.2f",
                       Design, cut.tox,cut.eff, cut.escalation,cut.tox.stagger)
    
    f1_com <- NULL
    for(i in 1:length(scen_num)){
      num_DLT <- NULL
      num_Respond <- NULL
      actual.Nsim <- 0
      for(j in 1:subset){
     
        oname1 <- paste(filename,"_scen_",scen_num[i],"_subset_",j, sep="")
        tryCatch({
          f1 <-  load(paste(oname1, ".RData", sep=""))
          
          
          for(c in 1:length(Adaptive_Stagger_sce)){
            if(!is.null(Adaptive_Stagger_sce[[c]])){
              num_DLT <- c(num_DLT, sum(Adaptive_Stagger_sce[[c]][[1]]$yT.true))
              num_Respond <- c(num_Respond, sum(Adaptive_Stagger_sce[[c]][[1]]$yE.true))
              actual.Nsim = actual.Nsim+1
            }
          }
          
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        #f1_com <- rbind(f1_com,f1)
      }
      
      mean_num_DLT = mean(num_DLT)
      mean_num_Respond = mean(num_Respond)
      
      f1_com <- rbind(f1_com, c(i, mean_num_DLT, mean_num_Respond, actual.Nsim))
    }
    
    
    combine.result <- data.frame(f1_com)
    colnames(combine.result) <- c("scenario", "# pts DLT", "# pts Respond", "actual.Nsim") 
    
  }
  #col_name <- c( paste("Dose ", seq(0:ndose)-1, sep = ""),"Trial.Duration","# pts","# Stagger","# stagger decision", paste("# pts in Dose ", seq(1:ndose), sep = ""), "% stagger",paste("% pts in Dose", 1:ndose,sep = ""))
  
  
  if(Design == "CUB"|Design == "CUB_AR3"|Design == "CUB_AR3"|Design == "CUB_N60"|
     Design == "CUB_Weibull"|Design == "CUB_U70"){
    # method <- paste(Design,"ctox",cut.tox, 'ceff',cut.eff, "ces", cut.escalation,sep="_")
    method <- Design
  }else if(Design == "Strict_Stagger"|Design == "Strict_Stagger_N60"|Design == "Backfill"|
           Design == "Strict_Stagger_AR3"|Design == "Strict_Stagger_AR7"|Design == "Strict_Stagger_quickdecay"|Design == "Strict_Stagger_nodecay"|
           Design == "Strict_Stagger_Weib"|Design == "Strict_Stagger_U70"){
    # method <- paste(Design,"ctox",cut.tox, 'ceff',cut.eff, "ces", cut.escalation, "cstag",cut.tox.stagger,sep="_")
    method <- Design
  }else if(Design == "Adaptive_Stagger"|Design == "Adaptive_Stagger_N60"|
           Design == "Adaptive_Stagger_quickdecay"|Design == "Adaptive_Stagger_AR3"|
           Design == "Adaptive_Stagger_AR7"|Design == "Adaptive_Stagger_nodecay"|
           Design == "Adaptive_Stagger_Weibull"|Design == "Adaptive_Stagger_U70"|
           Design == "Adaptive_Stagger_sens_A"|Design == "Adaptive_Stagger_sens_B"|
           Design == "Adaptive_Stagger_sens_no"|Design == "Adaptive_Stagger_Cor0.3"|
           Design == "Adaptive_Stagger_Cor0.5"|Design == "Adaptive_Stagger_Cor0.8"){
    # method <- paste(Design,"ctox",cut.tox, 'ceff',cut.eff, "ces", cut.escalation, "cstag",cut.tox.stagger, sep="_")
    method <- Design
  }
  
  combine.result$Method <- method
  
 
  return(combine.result)
}


