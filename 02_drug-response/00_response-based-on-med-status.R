################################################################################
#                          approach 3 for drug response                        #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(ggpubr)
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response"
setwd(project.dir)
################################################################################
recode_events <- function(event) {
  t <- data.frame(raw = event) %>%
    mutate(recoded = ifelse(event == "baseline_year_1_arm_1", 0,
                            ifelse(event == "1_year_follow_up_y_arm_1", 1,
                                   ifelse(event == "2_year_follow_up_y_arm_1", 2,
                                          ifelse(event == "3_year_follow_up_y_arm_1", 3,4)))))
  return(t$recoded)
}

################################################################################
################################################################################
abcd.raw.dir <- "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core"
abcd.deriv.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/abcd5"
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv")
colnames(abcd.pgs) <- sub("corrected_", "", colnames(abcd.pgs))
demo <- read_csv(paste0(abcd.deriv.dir, "/age-sex-by-eventname.csv"))
abcd.meds.r <- read_rds(paste0(abcd.deriv.dir, "/../meds/abcd5/abcd5-meds-matrix.rds"))
drop.sub <- read_csv(paste0(abcd.deriv.dir, "/../meds/abcd5/subjects-missing-med-name.csv"))
abcd.meds <- anti_join(abcd.meds.r, drop.sub)
abcd.pred <- read_rds("../data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(predicted = m) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
################################################################################
################################################################################
# identify responsiveness
mph.inc <- abcd.meds %>%
  select(IID, eventname, methylphenidate) %>%
  filter(methylphenidate == 1) %>%
  group_by(IID) %>%
  mutate(responsive = ifelse(sum(methylphenidate)>1, 1, 0)) %>%
  mutate(coded_eventname = recode_events(eventname))
################################################################################
# check correlation between responsiveness and PGS
mph.pgs <- inner_join(mph.inc %>% 
                        distinct(IID, .keep_all = T) %>%
                        select(IID, responsive, coded_eventname),
                      abcd.pgs) %>%
  pivot_longer(cols = colnames(abcd.pgs)[-1], names_to = "PGS", values_to = "pgs_value") 
mph.pgs %>%
  mutate(responsive = as.factor(responsive)) %>%
  ggplot(aes(x=responsive, y= pgs_value)) +
  geom_violin() +
  facet_wrap(~PGS, ncol = 12) +
  stat_compare_means(size = 2.5)
################################################################################
# check correlation between responsiveness and predicted response
mph.pred <- inner_join(mph.inc %>% 
                        distinct(IID, .keep_all = T) %>%
                        select(IID, responsive, coded_eventname),
                      abcd.pred)
mph.pred %>%
  mutate(responsive = as.factor(responsive)) %>%
  ggplot(aes(x=responsive, y= predicted)) +
  geom_violin() +
  stat_compare_means(size = 3)
