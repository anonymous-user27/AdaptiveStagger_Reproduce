###################################################################################################################
#                      Bayesian Phase I-II Designs with Adaptive Rules for Staggering Patient Entry
#             
#  This file R code to aggregate simulation results and get numbers presented in Table S13 (Time to event Weibull Distribution)
###################################################################################################################



source("Sim_codes_functions/Result_aggregate_code.R")
source("Sim_codes_functions/DLT_Response_count_aggregate_code.R")


Combined_sim_OC_results <- data.frame()
Combined_DLT_res_results <- data.frame()
############# Adaptive Stagger  ###############

Design = "Adaptive_Stagger_Weibull"
ndose = 5
Nsim=25
cut.tox <- c(0.1)
cut.eff <- c(0.05)
cut.escalation <- c(0.5)
cut.tox.stagger <- c(0.4)
subset = 200

result_sim_OC <-  aggregate_sim_OC_results(Design,cut.tox, cut.eff, cut.escalation, cut.tox.stagger,
                                           scen_num, subset,Nsim, ndose)
result_sim_OC
result_DLT_res <-  aggregate_DLT_Res_results(Design,cut.tox, cut.eff, cut.escalation, cut.tox.stagger,
                                             scen_num, subset,Nsim, ndose)
result_DLT_res
Combined_sim_OC_results <- rbind(Combined_sim_OC_results, result_sim_OC)
Combined_DLT_res_results <- rbind(Combined_DLT_res_results, result_DLT_res)

############# Strict Stagger  ###############

Design = "Strict_Stagger_Weibull"
ndose = 5
Nsim=50
cut.tox <- c(0.1)
cut.eff <- c(0.05)
cut.escalation <- c(0.5)
cut.tox.stagger <- c(0.4)
subset = 100


result_sim_OC <-  aggregate_sim_OC_results(Design,cut.tox, cut.eff, cut.escalation, cut.tox.stagger,
                                           scen_num, subset,Nsim, ndose)
result_sim_OC
result_DLT_res <-  aggregate_DLT_Res_results(Design,cut.tox, cut.eff, cut.escalation, cut.tox.stagger,
                                             scen_num, subset,Nsim, ndose)
result_DLT_res
Combined_sim_OC_results <- rbind(Combined_sim_OC_results, result_sim_OC)
Combined_DLT_res_results <- rbind(Combined_DLT_res_results, result_DLT_res)

############# Conventional Utility-Based  ###############

Design = "CUB_Weibull"
ndose = 5
Nsim=125
cut.tox <- c(0.1)
cut.eff <- c(0.05)
cut.escalation <- c(0.5)
cut.tox.stagger <- NULL
subset = 40


result_sim_OC <-  aggregate_sim_OC_results(Design,cut.tox, cut.eff, cut.escalation, cut.tox.stagger,
                                           scen_num, subset,Nsim, ndose)
result_sim_OC
result_DLT_res <-  aggregate_DLT_Res_results(Design,cut.tox, cut.eff, cut.escalation, cut.tox.stagger,
                                             scen_num, subset,Nsim, ndose)
result_DLT_res
Combined_sim_OC_results <- rbind(Combined_sim_OC_results, result_sim_OC)
Combined_DLT_res_results <- rbind(Combined_DLT_res_results, result_DLT_res)

############# BOIN12  ###############
Design = "BOIN12_Weibull"
filename <-sprintf("Results/%s_ctox_%.2f_ceff_%.2f.csv", 
                   Design, 1-cut.tox,1-cut.eff)
result_BOIN12 <- read.csv(filename)
col_name <- c("Scenario","Method", paste("Dose ", seq(0:ndose)-1, sep = ""),"Trial.Duration","# pts treated","# Stagger" ,"# stagger decision",
              paste("# pts in Dose ", seq(1:ndose), sep = ""), "% stagger",
              paste("% pts in Dose", 1:ndose,sep = ""), "actual.Nsim", "# pts DLT", "# pts Respond")

colnames(result_BOIN12) <- col_name

############# Combine results  ###############
Combined_sim_OC_results
Combined_DLT_res_results

merged_data <- merge(Combined_sim_OC_results, Combined_DLT_res_results, by = c("scenario", "Method"))

merged_data <- merged_data[,-dim(merged_data)[2]]
colnames(merged_data) <- col_name

Combined_results <- rbind(merged_data, result_BOIN12)

library(dplyr)
method_order <- c( "CUB_Weibull", "Strict_Stagger_Weibull", "BOIN12_Weibull","Adaptive_Stagger_Weibull")
Combined_results <- Combined_results %>%
  mutate(Method = factor(Method, levels = method_order)) %>%
  arrange(Scenario, Method)

write.csv(Combined_results,"TableS13/Sim_results_Weibull_TableS13.csv", row.names = FALSE)

