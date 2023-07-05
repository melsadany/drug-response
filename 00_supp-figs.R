################################################################################
#                               supplementary figs                            #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response"
setwd(project.dir)
################################################################################
spark.mph <- read_tsv("data/derivatives/spark-rm-mph-effect.tsv")
abcd.mph <- read_tsv("data/derivatives/abcd-sst-data-pcs.tsv")

samples.meta <- rbind(spark.mph %>% select(IID, age = Dependent_age, sex)%>%mutate(cohort = "SPARK"), 
                      abcd.mph %>% select(IID, age=interview_age, sex)%>%mutate(cohort = "ABCD")) %>%
  mutate(cohort = as.factor(cohort), age = as.numeric(age), sex = as.factor(sex))
p1<- samples.meta %>%
  ggplot(aes(x=age, fill=cohort)) +
  geom_histogram() +
  scale_fill_manual(values = redblu.col)
p2 <- samples.meta %>%
  ggplot(aes(x=sex, fill=cohort)) +
  geom_bar() +
  scale_fill_manual(values = redblu.col)
patchwork::wrap_plots(p1,p2)
################################################################################