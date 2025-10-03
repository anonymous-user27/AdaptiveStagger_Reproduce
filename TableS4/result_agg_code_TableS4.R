###################################################################################################################
#                      Bayesian Phase I-II Designs with Adaptive Rules for Staggering Patient Entry
#             
#  This file R code to aggregate simulation results and get numbers presented in  Table S4 (low tox not predictive of DLT)
###################################################################################################################



source("Sim_codes_functions/Result_aggregate_code.R")
source("Sim_codes_functions/DLT_Response_count_aggregate_code.R")


Combined_sim_OC_results <- data.frame()
Combined_DLT_res_results <- data.frame()
############# Adaptive Stagger  ###############

Design = "Adaptive_Stagger_nonpredict"
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


############# Combine results  ###############
Combined_sim_OC_results
Combined_DLT_res_results

merged_data <- merge(Combined_sim_OC_results, Combined_DLT_res_results, by = c("scenario", "Method"))

merged_data <- merged_data[,-dim(merged_data)[2]]
colnames(merged_data) <- col_name


write.csv(Combined_results,"TableS4/Sim_results_lowtox_nonpredict_TableS4.csv", row.names = FALSE)

