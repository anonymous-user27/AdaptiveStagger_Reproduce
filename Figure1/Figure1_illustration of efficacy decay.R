
library(ggplot2)

t <- c(0, 7,14,21,28)
decay <- 0.05

y1 <- round(c(0.4, 0.95*0.4,0.9*0.4,0.82*0.4,0.7*0.4),2)

df <- data.frame( time =c(0, 7,14,21,28),
                    y =  y1)

p <- ggplot(data=df, aes(x=time, y=y, group=1, label=y)) +
  geom_path()+
  geom_point() + 
  geom_text(nudge_y = 0.02, size=2.3)+
  ylim(0,0.5) +   
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 7,14,21,28)) +
  xlab("Treatment Delay Time (days)") + ylab("Efficacy Probability") +
  theme(axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.title.y=element_text(size=7))

p

ggsave("Figure1/Figure1_example_efficacy_decay.jpg", width = 3.2, height=3, limitsize = FALSE, dpi=1000)

