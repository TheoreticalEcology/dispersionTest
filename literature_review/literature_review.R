###### Script to perform literature reveiw in PUBMED ###### 
# 11 Dec 2025

library(europepmc)
library(tidyverse)
library(readxl)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot() + 
  theme(panel.background = element_rect(color="black")))


# saved objects
load("literature_review/literature_review.Rdata")

# QUERY WORDS ----

eco_query <- '("ecology" OR "ecolog*")'

glm_terms <- 'AND ("count data" OR "poisson" OR "negative binomial" OR "generalized poisson" OR "generalised poisson" OR "conway-maxwell poisson" OR "binomial" OR "beta-binomial" OR "binomial proportion")'

count_terms <- 'AND ("count data" OR "poisson" OR "negative binomial" OR "generalized poisson" OR "generalised poisson" OR "conway-maxwell poisson")'

prop_terms <- 'AND ("binomial" OR "beta-binomial" OR "binomial proportion")'

disp_terms <- 'AND ("overdispersion" OR "over dispersion" OR "over-dispersion" OR "underdispersion" OR "under dispersion" OR "under-dispersion" OR "dispersion")'

# TRENDS ----
resultEco <- epmc_hits_trend(query = eco_query, period=2005:2025,
                             synonym = FALSE)

resultGLM <- epmc_hits_trend(query = paste(eco_query, glm_terms), 
                             period=2005:2025, synonym = FALSE)

resultCount <- epmc_hits_trend(query = paste(eco_query, count_terms), 
                             period=2005:2025, synonym = FALSE)

resultProp <- epmc_hits_trend(query = paste(eco_query, prop_terms), 
                               period=2005:2025, synonym = FALSE)

resultGLMDisp <- epmc_hits_trend(query = paste(eco_query, glm_terms, 
                                            disp_terms), 
                              period=2005:2025, synonym = FALSE)

resultCountDisp <- epmc_hits_trend(query = paste(eco_query, count_terms,
                                               disp_terms), 
                                 period=2005:2025, synonym = FALSE)

resultPropDisp <- epmc_hits_trend(query = paste(eco_query, prop_terms,
                                                 disp_terms), 
                                   period=2005:2025, synonym = FALSE)

results <- resultEco %>% rename(eco_hits = query_hits) %>%
  mutate(glm_hits = resultGLM$query_hits,
         count_hits = resultCount$query_hits,
         prop_hits = resultProp$query_hits,
         glm_disp_hits = resultGLMDisp$query_hits,
         count_disp_hits = resultCountDisp$query_hits,
         prop_disp_hits = resultPropDisp$query_hits) %>%
  mutate(glm_within_eco = glm_hits*100/eco_hits,
         disp_within_glm = glm_disp_hits*100/glm_hits,
         disp_within_count = count_disp_hits*100/count_hits,
         disp_within_prop = prop_disp_hits*100/prop_hits
         ) 



## Figure trends ----

# only ecological studies using GlM for count data
figA <- ggplot(results, aes(x=year, y=glm_within_eco)) + 
  geom_point(col="coral")+
  geom_line(col="coral") +
  scale_y_continuous(limits=c(0,9),breaks=c(0,2,4,6,8),
                     name="Proportion of studies (%)") + 
  labs(tag="A)", title = "GLMs/GLMMs within Ecology")


# only dispersion within GlM studies
figB <- ggplot(results, aes(x=year, y=disp_within_glm)) + 
  geom_point(col="aquamarine3")+
  geom_line(col="aquamarine3") +
  scale_y_continuous(breaks=c(0,5,10,15,20,25),labels=c(0,5,10,15,20,25),
                     limits=c(0,25),
    name="Proportion of studies (%)") + 
  labs(tag="B)", title="Dispersion within GLMs/GLMMs")

figA + figB 
ggsave("figures/box1_trend_S1.1.jpeg", height = 4, width=8)


# count + prop dispersion with GLM studies
results %>% select(year, disp_within_count, disp_within_prop) %>%
  pivot_longer(2:3, names_to = "model", values_to = "Proportion") %>%
ggplot( aes(x=year, y=Proportion, col=model)) + 
  scale_color_manual(values=c("aquamarine3","coral"), labels=c("Count data", "Discrete proportions"),
                     name="")+
  geom_line() +
  geom_point()+
  scale_y_continuous(breaks=c(5,10,15,20,25,30),labels=c(5,10,15,20,25,30),
                     limits=c(5,30),
                     name="Proportion of studies (%)") +
  theme(legend.position = "inside",
        legend.text = element_text(size=13),
        legend.position.inside = c(0.5,0.2))
# legend: Annual trends for the proportion of ecological studies using GLMs for count AND/OR discrete proportion data that mention dispersion and related words in the text.  
ggsave("figures/box1_trend.jpeg", height = 4, width=5)




# PAPERS QUERY ----

papers <- epmc_search(query = paste(eco_query, glm_terms, disp_terms), 
                      limit=10000, synonym = FALSE)
# 7634 records


# journals
journals <- sort(unique(papers$journalTitle))

# ONLY eco / bot / zoo journals (identified manually)
ecoJ <- c("Wetlands (Wilmington)","Weed Res","Virus Evol", "Water Air Soil Pollut","Urban For Urban Green","Theor Popul Biol","Theor Ecol", "Sci Rep", "Sci Adv","Restor Ecol", "Proc Natl Acad Sci U S A", "Primates", "PLoS One", "Plants (Basel)", "Plant Environ Interact", "Plant Divers", "Philos Trans R Soc Lond B Biol Sci", "PeerJ",  "Ornithology Research", "Oecologia", "New Phytol", "Nature", "Nat Plants","Nat Ecol Evol" ,  "Mov Ecol", "Mol Ecol", "Mol Ecol Resour", "Microb Ecol", "Methods Ecol Evol", "Mar Biol", "Mamm Res", "Limnologica","Landsc Urban Plan", "Landsc Ecol","J Plant Res", "J Ornithol" , "J Mammal", "J Insect Sci", "J Insect Conserv",  "J Fish Biol" ,"J Evol Biol", "J Ecol", "J Appl Ecol", "J Anim Ecol", "Int J Evol Biol" ,   "Insects" ,"Glob Ecol Biogeogr", "Glob Chang Biol",  "Funct Ecol" , "For Ecol Manage", "Freshw Biol", "Evolution" ,"Evol Ecol Res", "Evol Ecol" ,  "Ecosphere"   , "Ecology.",  "Ecography", "Ecohealth","Ecol Appl", "Ecol Entomol","Ecol Evol","Ecol Lett", "Ecol Monogr","Ecol Res","Ecology", "Divers Distrib","Conserv Biol", "BMC Ecol", "BMC Ecol Evol", "Biodivers Conserv", "Behav Ecol", "Am Nat" )


# subseting papers:
# - only selected journals and 
# - from 2025
# - openAccess
papers25 <- papers %>% filter(journalTitle %in% ecoJ,
                              pubYear ==2025,
                              isOpenAccess == "Y")

# taking 200 randomly to read
set.seed(1)
papers200 <- papers25 %>% slice_sample(n=200)


# writing table ----
write_csv2(papers200, file="literature_review/literature_review2025.csv")

# saving objects ----
save(results,papers25,papers200, file="literature_review/literature_review.Rdata")


# REVIEWED LITERATURE ----

# with the randomly selected results from the literature searche above, we read the papers and manually retreived the following information:

# in_scope: y or n - if it is an ecological study
# count_prop_data: y or n  - if there is data analysis using count and/or discrete proportion data
# model: the distributions mentioned in the analysis
# dispersion_mentioned: y or n if dispersion (over or underdispersion) was mentioned in the text
# which_dispersion: which dispersion issue was mentioned: over for overdispersion, under for underdispersion, both for both types
# check_dispersion: y or n if they explicitly or not checked for dispersion issues  
# test_dispersion: y or no if they explicitly tested for dispersion issues 
# test_used: which test they say they used to test for dispersion issues 
# comments: comments on what they did
# OBS: extra observations from the papers

revlit <- read_excel("literature_review/literature_review2025.xlsx") %>%
  dplyr::select(id, title, in_scope,count_prop_data, model, dispersion_mentioned, 
         which_dispersion, check_dispersion, test_dispersion, test_used, 
         comments, OBS) %>%
  filter(!is.na(in_scope))
#155 reviewed papers

revlit100 <- revlit %>% filter(in_scope=="y", count_prop_data =="y")

# dispersion mentioned
revlit100 %>% count(dispersion_mentioned)

# which dispersion
dispy <- revlit100 %>% filter(dispersion_mentioned == "y")
dispy %>% count(which_dispersion)

# did they check overdispersion?
dispy %>% count(check_dispersion)

# test dispersion
dispy %>% count(test_dispersion)

# which tests
tests <- dispy %>% filter(test_dispersion == "y") %>% count(test_used)

(dharma <- sum(tests$n[grepl("DHARMa", tests$test_used)]))
dharma  / sum(tests$n)
(performance.test <- sum(tests$n[grepl("performance", tests$test_used)]))
performance.test/sum(tests$n)
(pearson <- sum(tests$n[grepl("pearson", tests$test_used)]))
pearson/sum(tests$n)
(model.comp <- sum(tests$n[grepl("model", tests$test_used)]))
model.comp/sum(tests$n)


# Distributions
distr <- revlit100 %>% dplyr::select(id, model) %>% 
  separate(2,c("m1","m2", "m3", "m4","m5"), sep=", ") %>%
  pivot_longer(2:6, names_to = "name", values_to="model",values_drop_na = TRUE) 
sort(table(distr$model), decreasing=T)

# 

