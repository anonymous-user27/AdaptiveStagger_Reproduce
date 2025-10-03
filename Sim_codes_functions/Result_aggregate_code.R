#############################################################
#  Function to aggregate the simulation results:
#   %select #patients per dose, Trial Duration, #Stagger
#############################################################
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


aggregate_sim_OC_results <- function(Design,cut.tox, cut.eff, cut.escalation, cut.tox.stagger, scen_num, subset,Nsim, ndose){
      
  
  col_name <- c( paste("Dose ", seq(0:ndose)-1, sep = ""),"Trial.Duration","# pts treated","# Stagger" ,"# stagger decision",
                  paste("# pts in Dose ", seq(1:ndose), sep = ""), "% stagger",
                  paste("% pts in Dose", 1:ndose,sep = ""))
  
  
  if(Design == "CUB"|Design == "CUB_AR3"|Design == "CUB_AR3"|Design == "CUB_N60"|
     Design == "CUB_Weibull"|Design == "CUB_U70"){
    filename <-sprintf("Results/%s_ctox_%.2f_ceff_%.2f_cesc_%.2f", 
                       Design, cut.tox,cut.eff, cut.escalation)
    col_name1 <- c( paste("Dose ", seq(0:ndose)-1, sep = ""),"Trial.Duration","# pts treated","# Stagger","# stagger decision"
                    , paste("# pts in Dose ", seq(1:ndose), sep = ""), "% stagger",paste("% pts in Dose", 1:ndose,sep = ""))
    
    f_list1 <- list()
    for(i in 1:length(scen_num)){
      f1_com <- NULL
      for(j in 1:subset){
        oname1 <- paste(filename,"_scen_",scen_num[i],"_subset_",j, sep="")
        f1 <- NULL
        tryCatch({
          f1 <- assign(oname1, read.csv(paste(oname1, ".csv", sep="")))
          f1 <- as.matrix(f1)
          if(!is.null(f1)){f1[,4:(dim(f1)[2]-1)] <- f1[,4:(dim(f1)[2]-1)]*f1[,"actual.Nsim"]}
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        f1_com <- rbind(f1_com,f1)
      }
      if(!is.null(f1_com)){
        f1_com <- data.frame(f1_com)
        f1_com <- f1_com[!duplicated(f1_com),]
        act_Nsim = sum(f1_com$actual.Nsim)
        avg_result = matrix(apply(f1_com[,4:(dim(f1_com)[2]-1)],2,sum)/act_Nsim,nrow=1 )
        avg_result = matrix(c(scen_num[i],avg_result, act_Nsim = act_Nsim),nrow=1)
        colnames(avg_result) <- c("scenario", col_name1, "actual.Nsim")
        f_list1[[i]] <- avg_result}
    }
    com <- do.call(rbind,f_list1)           # rbind a list of matrix
    combine.result <- data.frame(com)
    colnames(combine.result) <- c("scenario", col_name1, "actual.Nsim") 
    
    combine.result <- combine.result[, c("scenario", col_name, "actual.Nsim")]
    
  }else if(Design == "Strict_Stagger"|Design == "Strict_Stagger_N60"|Design == "Backfill"|
           Design == "Strict_Stagger_AR3"|Design == "Strict_Stagger_AR7"|Design == "Strict_Stagger_quickdecay"|Design == "Strict_Stagger_nodecay"|
           Design == "Strict_Stagger_Weib"|Design == "Strict_Stagger_U70"){
    filename <-sprintf("Results/%s_ctox_%.2f_ceff_%.2f_ces_%.2f_cstag_%.2f", 
                       Design, cut.tox,cut.eff, cut.escalation, cut.tox.stagger)
    f_list1 <- list()
    for(i in 1:length(scen_num)){
      f1_com <- NULL
      for(j in 1:subset){
        oname1 <- paste(filename,"_scen_",scen_num[i],"_subset_",j, sep="")
        f1 <- NULL
        tryCatch({
          f1 <- assign(oname1, read.csv(paste(oname1, ".csv", sep="")))
          f1 <- as.matrix(f1)
          if(!is.null(f1)){f1[,4:(dim(f1)[2]-1)] <- f1[,4:(dim(f1)[2]-1)]*f1[,"actual.Nsim"]}
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        f1_com <- rbind(f1_com,f1)
      }
      if(!is.null(f1_com)){
        f1_com <- data.frame(f1_com)
        f1_com <- f1_com[!duplicated(f1_com),]
        act_Nsim = sum(f1_com$actual.Nsim)
        avg_result = matrix(apply(f1_com[,4:(dim(f1_com)[2]-1)],2,sum)/act_Nsim,nrow=1 )
        avg_result = matrix(c(scen_num[i],avg_result, act_Nsim = act_Nsim),nrow=1)
        colnames(avg_result) <- c("scenario", col_name, "actual.Nsim")
        
        
        f_list1[[i]] <- avg_result}
    }
    com <- do.call(rbind,f_list1)           # rbind a list of matrix
    combine.result <- data.frame(com)
    colnames(combine.result) <- c("scenario", col_name, "actual.Nsim") 
    
  }else if(Design == "Adaptive_Stagger"|Design == "Adaptive_Stagger_N60"|
           Design == "Adaptive_Stagger_quickdecay"|Design == "Adaptive_Stagger_AR3"|
           Design == "Adaptive_Stagger_AR7"|Design == "Adaptive_Stagger_nodecay"|
           Design == "Adaptive_Stagger_Weibull"|Design == "Adaptive_Stagger_U70"|
           Design == "Adaptive_Stagger_sens_A"|Design == "Adaptive_Stagger_sens_B"|
           Design == "Adaptive_Stagger_sens_no"|Design == "Adaptive_Stagger_Cor0.3"|
           Design == "Adaptive_Stagger_Cor0.5"|Design == "Adaptive_Stagger_Cor0.8"){
    filename <-sprintf("Results/%s_ctox_%.2f_ceff_%.2f_ces_%.2f_cstag_%.2f",
                       Design, cut.tox,cut.eff, cut.escalation,cut.tox.stagger)
    
    f_list1 <- list()
    for(i in 1:length(scen_num)){
      f1_com <- NULL
      for(j in 1:subset){
        oname1 <- paste(filename,"_scen_",scen_num[i],"_subset_",j, sep="")
        f1 <- NULL
        tryCatch({
          f1 <- assign(oname1, read.csv(paste(oname1, ".csv", sep="")))
          f1 <- as.matrix(f1)
          if(!is.null(f1)){f1[,4:(dim(f1)[2]-1)] <- f1[,4:(dim(f1)[2]-1)]*f1[,"actual.Nsim"]}
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        f1_com <- rbind(f1_com,f1)
      }
      if(!is.null(f1_com)){
        f1_com <- data.frame(f1_com)
        f1_com <- f1_com[!duplicated(f1_com),]
        act_Nsim = sum(f1_com$actual.Nsim)
        avg_result = matrix(apply(f1_com[,4:(dim(f1_com)[2]-1)],2,sum)/act_Nsim,nrow=1 )
        avg_result = matrix(c(scen_num[i],avg_result, act_Nsim = act_Nsim),nrow=1)
        colnames(avg_result) <- c("scenario", col_name, "actual.Nsim")
        f_list1[[i]] <- avg_result}
    }
    com <- do.call(rbind,f_list1)           # rbind a list of matrix
    combine.result <- data.frame(com)
    colnames(combine.result) <- c("scenario", col_name, "actual.Nsim") 
  }
 
  
  if(Design == "CUB"|Design == "CUB_AR3"|Design == "CUB_AR3"|Design == "CUB_N60"|
     Design == "CUB_Weibull"|Design == "CUB_U70"){
    #method <- paste(Design,"ctox",cut.tox, 'ceff',cut.eff, "ces", cut.escalation,sep="_")
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






