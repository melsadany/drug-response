################################################################################
#                   predcit drug response for ABCD participants                #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response"
setwd(project.dir)
################################################################################
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meph-meta_by-eventname.rds") %>%
  filter(base_meph == 1 | year1_meph == 1 | year2_meph == 1 | year3_meph == 1)
################################################################################
################################################################################
load("data/derivatives/m-outputs/abcd/combined-output.RData")
m.output <- to.save[[1]]
rm(to.save)
gc()
################################################################################
################################################################################
################################################################################
######################## correlate response to ABCL att ########################
################################################################################
################################################################################
################################################################################
# get correlation between taking MPH and att_problems total in ABCD


# this is weird. they only have ABCl collected once per subject!!
# stop the analysis here

abcd.abcl <- read_tsv("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_4_0/abcd_abcls01.txt")[-1,] %>%
  dplyr::select(IID = subjectkey, starts_with("interview"), sex, eventname, ends_with("_r")) %>%
  drop_na() %>%
  as.data.frame() %>%
  filter(IID %in% m.output$IID) %>%
  select(1:5,abcl_scr_prob_attention_r)
abcd.abcl <- cbind(abcd.abcl[c(1,2,4,5)], abcd.abcl[c(3,6)]%>%mutate_all(function(x) as.numeric(x)))

meds.abcl.meta <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meds-matrix_CMAP-meds.rds") %>%
  as.data.frame() %>%
  filter(ID %in% abcd.abcl$IID) %>%
  select(1:6, methylphenidate)

meta <- left_join(left_join(abcd.abcl, meds.abcl.meta%>%rename(IID=ID, age_med = interview_age)), 
                  m.output)
meta %>%
  filter(methylphenidate==1) %>%
  pivot_longer(cols = starts_with("model_"), names_to = "model", values_to = "predicted_response") %>%
  filter(model %in% paste0("model_", c(1,3,5,7,9,11,13,15))) %>%
  ggplot(aes(x=predicted_response, y=abcl_scr_prob_attention_r)) +
  geom_point()+
  geom_smooth(method = "lm") +
  facet_wrap("model", scales = "free") +
  stat_cor(color = "red")

################################################################################



################################################################################


################################################################################

################################################################################