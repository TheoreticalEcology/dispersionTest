# Dispersion test paper
# summarising figure
# Melina 
# set 25

library(tidyverse)
#install.packages("fmsb")
library(fmsb)

# scores 1 to 3
# df <- data.frame(
#   test = c("parametric Pearson Residuals", "nonparametric Pearson Residuals",
#            "simulation-based response variance"),
#   speed = c(3,1,2),
#   GLMs = c(3,3,3),
#   GLMs_special_cases = c(1,2,2),
#   GLMMs_few_REgroups = c(2,3,3),
#   GLMMs_many_REgroups = c(1,3,3)
# )

# scores 1 to 5
df <- data.frame(
  test = c("parametric Pearson Residuals", "nonparametric Pearson Residuals",
           "simulation-based response variance"),
  speed = c(5,2,3),
  GLMs = c(5,5,5),
  GLMs_special_cases = c(1,4,3),
  GLMMs_few_REgroups = c(2,5,4),
  GLMMs_many_REgroups = c(1,4,4)
)

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
