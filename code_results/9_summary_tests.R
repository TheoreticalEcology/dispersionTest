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
  Speed = c(3,1,2),
  GLM = c(3,3,3),
  `GLM (small-data)` = c(1,2,1),
  `GLMM (few RE groups)` = c(2,3,3),
  `GLMM (many RE groups)` = c(1,3,2)
)

# plot plot

df %>% pivot_longer(2:6, names_to = "variable", values_to = "value") %>%
  mutate(test = fct_relevel(test,
                            "parametric Pearson Residuals",
                            "nonparametric Pearson Residuals",
                            "simulation-based response variance"
                            ),
         variable = fct_relevel(variable,"GLM","GLM (small-data)",
                                "GLMM (few RE groups)", "GLMM (many RE groups)", 
                                "Speed" )) %>%
  mutate(test.n = as.numeric(test),
         variable.n = as.numeric(variable)) -> df2

df2$value.t <- "++"
df2$value.t[df2$value == 1] <- "-"
df2$value.t[df2$value == 2] <- "+"

ggplot(df2, aes(x=variable.n, y=test.n, fill=value, col=value))+
  geom_text(aes(label= value.t), size=15) +
  #geom_point(shape=21,size=20) +
  #scale_radius(range = c(3, 25),guide = "none" ) +
  scale_fill_gradient(low="coral1", high="aquamarine3", name = "Score") +
  scale_color_gradient(low="coral1", high="aquamarine3", name = "Score") +
  scale_y_continuous(name = "", breaks=1:3, limits=c(0.6,3.4),
                     labels= c("Parametric \n Pearson  residuals",
                               "Nonparametric \n Pearson residuals",
                               "Simulation-based \n response variance"
                               )) +
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
        axis.text = element_text(size=12.5),
        legend.position = "none")



ggsave(filename = "figures/9_testsEvaluation_sign.jpeg", height = 4.5, width = 9)
