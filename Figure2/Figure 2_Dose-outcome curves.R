library(ggplot2)
library(latex2exp)
library(gridExtra)
library(cowplot)


source("Sim_codes/Sim.scenarios.R")
library(locfit)

dose_level <- c(1,2,3,4,5)


Scenarios

eff_list <- list(
  eff1 <- Scenarios[[1]]$prob_E,
  eff2 <- Scenarios[[2]]$prob_E,
  eff3 <- Scenarios[[3]]$prob_E,
  eff4 <- Scenarios[[4]]$prob_E,
  eff5 <- Scenarios[[5]]$prob_E,
  eff6 <- Scenarios[[6]]$prob_E
)

tox_list <- list(
  tox1 <- Scenarios[[1]]$prob_T,
  tox2 <- Scenarios[[2]]$prob_T,
  tox3 <- Scenarios[[3]]$prob_T,
  tox4 <- Scenarios[[4]]$prob_T,
  tox5 <- Scenarios[[5]]$prob_T,
  tox6 <- Scenarios[[6]]$prob_T
)

lowtox_list <- list(
  lowtoxtox1 <- Scenarios[[1]]$lambda_d * 30,
  lowtoxtox2 <- Scenarios[[2]]$lambda_d * 30,
  lowtoxtox3 <- Scenarios[[3]]$lambda_d * 30,
  lowtoxtox4 <- Scenarios[[4]]$lambda_d * 30,
  lowtoxtox5 <- Scenarios[[5]]$lambda_d * 30,
  lowtoxtox6 <- Scenarios[[6]]$lambda_d * 30
)



data.list <- list()
for(i in 1:length(eff_list)){
  
  data.list[[((i))]] <- data.frame(
    'dose' = dose_level, 
    'efficacy' = unlist(eff_list[[i]]), 
    'toxicity' = unlist(tox_list[[i]]),
    'Lowgradetoxicity' = unlist(lowtox_list[[i]])
  )
  
  
}
ggplot(data.list[[1]]) + 
  geom_line(aes(x = dose, y = efficacy), color = "darkred", size = 1.2) + 
  geom_line(aes(x = dose, y = toxicity), color = "darkblue", size = 1.2) + 
  geom_point(aes(x = dose, y = efficacy), size = 4) +
  geom_point(aes(x = dose, y = toxicity), size = 4, shape = 2) + 
  geom_hline(yintercept = 0.25, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = 0.35, color = "darkblue", linetype = "dashed") +
  theme_bw() + 
  labs(x = TeX(r'(log(Dose))'),
       title = paste0("Scenario", 1)) + 
  scale_x_continuous(breaks = c(1,2,3,4,5)) + 
  geom_line(aes(x = dose, y = Lowgradetoxicity/80), color = "darkgreen", size = 1.2) +
  geom_point(aes(x = dose, y = Lowgradetoxicity/80), shape = 5, color = "darkgreen", size = 4) + 
  scale_y_continuous(
    name = "Toxicity or Efficacy \n Probability", limits = c(0,1),
    sec.axis = sec_axis(~ . * 80, name = "Low-grade Intensity")
  )

legend_plot <- ggplot(data.list[[1]]) + 
  geom_hline(yintercept = 0.25, color = "darkred", linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = 0.35, color = "darkblue", linetype = "dashed", size = 0.3) +
  
  geom_line(aes(x = dose, y = efficacy, color = "Efficacy"), size = 0.3) + 
  geom_line(aes(x = dose, y = toxicity, color = "DLT"), size = 0.3) + 
  geom_line(aes(x = dose, y = Lowgradetoxicity / 80, color = "Low-grade Toxicities"), size = 0.3) +
  
  geom_point(aes(x = dose, y = efficacy, color = "Efficacy", shape = "Efficacy"), size = 0.8, fill = "white") +
  geom_point(aes(x = dose, y = toxicity, color = "DLT", shape = "DLT"), size = 0.8, fill = "white") + 
  geom_point(aes(x = dose, y = Lowgradetoxicity / 80, color = "Low-grade Toxicities", shape = "Low-grade Toxicities"), size = 0.8, fill = "white") + 
  
  theme_bw() + 
  labs(x = TeX(r'(log(Dose))'),
       title = paste0("Scenario ", 1),
       color = "Outcome",
       shape = "Outcome") + 
  scale_x_continuous(breaks = c(1,2,3,4,5)) + 
  
  scale_color_manual(
    values = c("Efficacy" = "darkred", "DLT" = "darkblue", "Low-grade Toxicities" = "darkgreen")
  ) +
  scale_shape_manual(
    values = c("Efficacy" = 21, "DLT" = 24, "Low-grade Toxicities" = 23)
  ) +
  scale_y_continuous(
    name = "DLT or Efficacy \n Probability", limits = c(0,1),
    sec.axis = sec_axis(~ . * 80, name = "Number of Low-grade Toxicities")
  ) +
  theme(
    axis.text.x = element_text(size = 6), 
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    plot.title = element_text(size = 8, face = "bold",hjust = 0.5),
    legend.text = element_text(size = 7.5),
    legend.title = element_text(size = 8),
    legend.position = "bottom" # ?????????????????????
  )

plot.list <- lapply(seq_len(6), function(i) {
  ggplot(data.list[[i]]) + 
    geom_hline(yintercept = 0.25, color = "darkred", linetype = "dashed", size = 0.3) +
    geom_hline(yintercept = 0.35, color = "darkblue", linetype = "twodash", size = 0.3) +
    
    geom_line(aes(x = dose, y = efficacy, color = "Efficacy"), size = 0.3) + 
    geom_line(aes(x = dose, y = toxicity, color = "DLT"), size = 0.3) + 
    geom_line(aes(x = dose, y = Lowgradetoxicity / 80, color = "Low-grade Toxicities"), size = 0.3) +
    
    geom_point(aes(x = dose, y = efficacy, color = "Efficacy", shape = "Efficacy"), size = 0.8, fill = "white") +
    geom_point(aes(x = dose, y = toxicity, color = "DLT", shape = "DLT"), size = 0.8, fill = "white") + 
    geom_point(aes(x = dose, y = Lowgradetoxicity / 80, color = "Low-grade Toxicities", shape = "Low-grade Toxicities"), size = 0.8, fill = "white") + 
    
    theme_bw() + 
    labs(x = 'Dose Level',
         title = paste0("Scenario ", i),
         color = "Outcome",
         shape = "Outcome") + 
    scale_x_continuous(breaks = c(1,2,3,4,5)) + 
    
    scale_color_manual(
      values = c("Efficacy" = "darkred", "DLT" = "darkblue", "Low-grade Toxicities" = "darkgreen")
    ) +
    scale_shape_manual(
      values = c("Efficacy" = 21, "DLT" = 24, "Low-grade Toxicities" = 23)
    ) +
    scale_y_continuous(
      name = "DLT or Efficacy \n Probability", limits = c(0,1),
      sec.axis = sec_axis(~ . * 80, name = "Number of \n Low-grade Toxicities")
    ) +
    theme(
      axis.text.x = element_text(size = 6), 
      axis.text.y = element_text(size = 6),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      plot.title = element_text(size = 8, face = "bold",hjust = 0.5),
      legend.position = "none" # ?????????????????????
    )
})

# 
 legend <- get_legend(legend_plot + theme(legend.position = "bottom"))

# combine grid of plots + legend
combined_plot <- plot_grid(
  plot_grid(plotlist = plot.list, ncol = 2, nrow = 3),
  legend,
  ncol = 1,
  rel_heights = c(1, 0.08)  # tweak as you like
)

ggsave("Figure2/Figure2_dose_outcome_curve.jpg", combined_plot, limitsize = FALSE, width = 5, height = 6, dpi = 1000)


