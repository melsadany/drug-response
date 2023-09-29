################################################################################
#                               cbcl dates infer                               #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response"
setwd(project.dir)
################################################################################
# I'm tying to infer the date that cbcl was reported in in
# I'm taking the age as a reference from cbcl and sleep data from rm
sleep.rm <- read_tsv("data/rawdata/spark-rm/raw/ChildSurvey20181031.csv") %>%
  mutate(IID = ParticipantID, sleep_rm_age_mos = Dependent_age_at_eval_start_months)
tr.cbcl <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/cbcl/static_data/CBCL_item-level_imputed.tsv")
inferred.cbcl <- inner_join(sleep.rm %>% select(IID, sleep_rm_age_mos, date = StartTime) , 
                           tr.cbcl %>% select(IID, cbcl_age = age_months)) %>%
  mutate(age_diff = cbcl_age - sleep_rm_age_mos) %>% 
  mutate(sleep_date = mdy(sub(" .*", "", date))) %>%
  mutate(cbcl_date = sleep_date%m+% months(age_diff)) %>%
  mutate(birthdate = sleep_date%m-% months(sleep_rm_age_mos))
write_tsv(inferred.cbcl, "data/derivatives/inferred-cbcl-dates-from-age-by-sleep-rm.tsv")
################################################################################