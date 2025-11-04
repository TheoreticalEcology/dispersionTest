# Dispersion test paper
# summarising figure
# 
# set 25

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())


# scores 1 to 3
df <- tibble(
  test = c("parametric Pearson Residuals", "nonparametric Pearson Residuals",
           "simulation-based response variance"),
  Speed = c(3,2,1),
  GLM = c(3,3,3),
  `GLM (small-data)` = c(1,2,1),
  `GLMM (few RE groups)` = c(2,3,3),
  `GLMM (many RE groups)` = c(1,3,2)
)

# plot plot

df %>% pivot_longer(2:6, names_to = "variable", values_to = "value") %>%
  mutate(test = fct_relevel(test,"simulation-based response variance",
                            "nonparametric Pearson Residuals",
                            "parametric Pearson Residuals"),
         variable = fct_relevel(variable,"GLM","GLM (small-data)",
                                "GLMM (few RE groups)", "GLMM (many RE groups)", 
                                "Speed" )) %>%
  mutate(test.n = as.numeric(test),
         variable.n = as.numeric(variable)) -> df2

df2$value.t <- "++"
df2$value.t[df2$value == 1] <- "-"
df2$value.t[df2$value == 2] <- "+"

ggplot(df2, aes(x=variable.n, y=test.n, fill=value, col=value))+
  geom_text(aes(label= value.t), size=15)+
  #geom_point(shape=21,size=20) +
  #scale_radius(range = c(3, 25),guide = "none" ) +
  scale_fill_gradient(low="coral1", high="aquamarine3", name = "Score") +
  scale_color_gradient(low="coral1", high="aquamarine3", name = "Score") +
  scale_y_continuous(name = "", breaks=1:3, limits=c(0.6,3.4),
                     labels= c("Parametric \n Pearson  residuals",
                               "Nonparametric \n Pearson residuals",
                               "Simulation-based \n response variance")) +
  scale_x_continuous(name="",  breaks=1:5, limits=c(0.7,5.3),
                     labels= c("GLM \n", "GLM \n (\"small-data\")",
                               "GLMM \n (few RE groups)", "GLMM \n (many RE groups)",
                               "Speed \n"),
                     position = "top") +
  geom_vline(xintercept = c(1.5,2.5,3.5,3.5,4.5), color="white", size=1.5)+
  geom_hline(yintercept = c(1.5,2.5), color="white", size=1.5)+
  theme(panel.background = element_rect(colour="white", fill="gray95"),
        axis.line=element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")




ggsave(filename = "figures/9_testsEvaluation_sign.jpeg", height = 4.5, width = 9)












#### radarChart 


dft <- df %>% pivot_longer(2:6,
                           names_to = "variables",
                           values_to = "score") %>%
  pivot_wider(names_from = variables, values_from = score)

dfbasic <- data.frame(dft[1:2,-1])
dfbasic[1,] <- 5
dfbasic[2,] <- 0

dfbasic <- bind_rows(dfbasic, dft[,-1])

radarchart(dfbasic) 

areas <- c(rgb(1, 0, 0, 0.25),
           rgb(0, 1, 0, 0.25),
           rgb(0, 0, 1, 0.25))

radarchart(dfbasic,
           cglty = 1,       # Grid line type
           cglcol = "gray", # Grid line color
           pcol = 2:4,      # Color for each line
           plwd = 2,        # Width for each line
           plty = 1,        # Line type for each line
           pfcol = areas)   # Color of the areas   

legend("topright",
       legend = df$test,
       bty = "n", pch = 20, col = areas,
       text.col = "grey25", pt.cex = 2)

par(mfrow=c(1,3))
radarchart(dfbasic[1:3,],
           cglty = 1,       # Grid line type
           cglcol = "gray", # Grid line color
           pcol = 2,      # Color for each line
           plwd = 2,        # Width for each line
           plty = 1,        # Line type for each line
           pfcol = areas[1],
           title="Parametric Pearson Residuals test")
radarchart(dfbasic[c(1,2,4),],
           cglty = 1,       # Grid line type
           cglcol = "gray", # Grid line color
           pcol = 3,      # Color for each line
           plwd = 2,        # Width for each line
           plty = 1,        # Line type for each line
           pfcol = areas[2],
           title="Nonparametric Pearson Residuals test")
radarchart(dfbasic[c(1,2,5),],
           cglty = 1,       # Grid line type
           cglcol = "gray", # Grid line color
           pcol = 4,      # Color for each line
           plwd = 2,        # Width for each line
           plty = 1,        # Line type for each line
           pfcol = areas[3],
           title="Simulation-based response variance test")

