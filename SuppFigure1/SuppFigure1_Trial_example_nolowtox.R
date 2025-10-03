# Load required packages
library(ggplot2)
library(dplyr)
library(ggforce)  # for rounded rectangles if desired

# Read the data
df <- read.csv("SuppFigure1/Trial_illustration_nolowtox_FigureS1.csv")  


# Convert Dose and Initial dose to factors for plotting
df$Dose <- as.factor(df$Dose)
df$Initial.dose <- as.factor(df$Initial.dose)

# Color for YE
df$fill_color <- ifelse(df$YE == 1, "green", "orange")


# apply rule: if VT < VE, then set YE = 0 and VE = VT
df$YE <- ifelse(df$VT < df$VE, 0, df$YE)
df$VE <- ifelse(df$VT < df$VE, df$VT, df$VE)


## if the patient die before response, then we need record YE=0 at VE=VT



df$OutcomeCombo <- with(df, factor(
  ifelse(YT == 0 & YE == 1, "No DLT + Response",
         ifelse(YT == 0 & YE == 0, "No DLT + No Response",
                ifelse(YT == 1 & YE == 1, "DLT + Response", "DLT + No Response"))),
  levels = c("No DLT + Response", "No DLT + No Response", "DLT + Response", "DLT + No Response")
))

# Add jittered y position
# Add jitter to y-position
set.seed(42)
df$y_jittered <- as.numeric(df$Dose) + runif(nrow(df), -0.25, 0.25)

# Base size for the square
square_size <- 1.0

# Plot
df <- df %>%
  arrange(t.enrollment) %>%
  mutate(
    cohort_index = rep(1:(n()/3), each = 3)[1:n()],
    within_cohort = rep(c(-0.2, 0, 0.2), length.out = n()),
    y_jittered = as.numeric(Dose) + within_cohort
  )


df$StaggerStatus <- ifelse(
  df$idx_stagger == 1 & df$t.enrollment != df$t.treated, "Stagger",
  ifelse(df$idx_stagger == 1 & df$t.enrollment == df$t.treated, "No stagger", NA)
)
# Square dimensions
half_width <- 5.5
half_height <- 0.08

# Plot
p <- ggplot() +
  # --- Dashed squares for pending treatment ---
  # Dashed rectangle for pending treatment
  geom_rect(data = df %>% filter(t.enrollment != t.treated),
            aes(xmin = t.enrollment - half_width, xmax = t.enrollment + half_width,
                ymin = as.numeric(Initial.dose) + within_cohort - half_height,
                ymax = as.numeric(Initial.dose) + within_cohort + half_height),
            color = "darkred", fill = NA, linetype = "dashed") +
  
  # ID inside dashed rectangle
  geom_text(data = df %>% filter(t.enrollment != t.treated),
            aes(x = t.enrollment, y = as.numeric(Initial.dose) + within_cohort, label = ID),
            size = 4, color = "black") +
  
  
  # --- Solid squares for treated patients ---
  geom_rect(data = df,
            aes(xmin = t.treated - half_width, xmax = t.treated + half_width,
                ymin = y_jittered - half_height, ymax = y_jittered + half_height,
                fill = OutcomeCombo),
            color = "black")  +
  scale_fill_manual(
    name = NULL,
    values = c( "No DLT + Response" = "#90ee90",
                "No DLT + No Response" = "#6baed6",
                "DLT + Response" = "#f4a582",
                "DLT + No Response" = "#e34a33")
  )+
  # ID inside solid squares
  geom_text(data = df,
            aes(x = t.treated, y = y_jittered, label = ID),
            size = 4, color = "black") +
  
  # --- Lines to outcome markers ---
  # 1. Line: from right edge of square to event time
  # Toxicity line (above center)
  geom_segment(data = df %>% filter(VT > 0),
               aes(x = t.treated + half_width,
                   y = y_jittered + 0.05,         # shift up
                   xend = t.treated + VT,
                   yend = y_jittered + 0.05),
               color = "black") +
  
  # Toxicity marker: "x" or "0"
  geom_text(data = df %>% filter(VT > 0),
            aes(x = t.treated + VT + 0.8,
                y = y_jittered + 0.05,
                label = ifelse(YT == 1, "x", "0")),
            size = 4, color = "black", vjust = 0.3)+
  # Response line (below center)
  geom_segment(data = df %>% filter(VE > 0),
               aes(x = t.treated + half_width,
                   y = y_jittered - 0.05,         # shift down
                   xend = t.treated + VE,
                   yend = y_jittered - 0.05),
               color = "darkblue", linetype = "longdash") +
  # Response marker: "r" for responders, open circle for non-responders
  geom_text(data = df %>% filter(VE > 0 & YE == 1),
            aes(x = t.treated + VE + 0.6,    # smaller horizontal offset
                y = y_jittered - 0.05),      # match response line y-shift
            label = "R",
            size = 3,
            color = "darkblue",
            vjust = 0.4,                     # vertical tweak (optional)
            fontface = "bold") +
  geom_point(data = df %>% filter(VE > 0 & YE == 0),
             aes(x = t.treated + VE + 0.8,
                 y = y_jittered - 0.05),
             shape = 21,       # circle with border
             size = 0.6,         # bold appearance
             fill = "darkblue",    # fill color
             color = "darkblue",  # border color
             stroke = 1)     +  # border thickness
  # --- Styling ---
  # scale_fill_identity() +
  scale_y_continuous(breaks = 1:5, name = "Dose Level") +
  labs(x = "Days Since Study Start", y = "Dose Level")+
  theme_minimal() +
  theme(panel.grid.minor = element_blank())+
  scale_x_continuous(
    name = "Days Since Study Start",
    breaks = seq(0, 450, 50),
    minor_breaks = seq(0, 480, 25)
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "plain", size = 15),
    axis.title.y = element_text(face = "plain", size = 15),
    axis.text.x = element_text(face = "plain", size = 15),
    axis.text.y = element_text(face = "plain", size = 15),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),
    legend.position = "None",
  ) + guides(
    fill = guide_legend(override.aes = list(
      shape = NA,        # No shape
      colour = NA,       # No stroke
      alpha = 1,         # Fully visible fill
      size = 6           # Size of color tile
    )),
    color = guide_legend(override.aes = list(
      shape = 17,
      size = 4
    ))
  )

# Display plot
print(p)

ggsave("SuppFigure1/SuppFigure1_Trial_example_nolowtox.jpg", p, limitsize = FALSE, width = 10, height = 8, dpi = 1000)
getwd()

