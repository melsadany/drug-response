################################################################################
#                 change in other neurocog tasks based on drug status          #
################################################################################
# packages setup ----------------------------------------------------------
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(ggpubr);library(ggExtra);library(ggh4x);library(ggridges)
####
# project dir setup -------------------------------------------------------
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response"
setwd(project.dir)
####
# ABCD demo and directory -------------------------------------------------
abcd.raw.dir <- "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core"
age <- read_csv(paste0(abcd.raw.dir, "/abcd-general/abcd_y_lt.csv")) %>% 
  select(IID = src_subject_id, eventname, interview_age)
sex <- read_csv(paste0(abcd.raw.dir, "/gender-identity-sexual-health/gish_p_gi.csv")) %>%
  mutate(sex = ifelse(demo_sex_v2 == 2, "Female", ifelse(demo_sex_v2 == 1, "Male", ifelse(demo_sex_v2 == 3, "intersex_M", NA)))) %>%
  select(IID = src_subject_id, sex) %>%
  distinct(IID, .keep_all = T)
# combine age and sex of all participants by eventname in one dataframe
demo <- full_join(age, sex)
####
# ABCD meds ---------------------------------------------------------------
# list of adhd meds of interest
all.adhd.meds <- data.frame(drug = c("methylphenidate", "adderall", "concerta", "vyvanse", 
                                     "dextroamphetamine",  "ritalin", "intuniv", "strattera",
                                     "tenex", "amphetamine", "dexmethylphenidate", "lisdexamfetamine",
                                     "atomoxetine", "clonidine", "guanfacine", "albuterol"
))
adhd.meds <- data.frame(drug = c("methylphenidate",
                                 "concerta",
                                 "stim", "non_stim", "neither",
                                 # "dextroamphetamine", "dexmethylphenidate", 
                                 "amphetamine", "lisdexamfetamine",
                                 "intuniv",
                                 "atomoxetine", "clonidine", "guanfacine", "albuterol"
))
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/abcd5-meds-matrix.rds") %>%
  as.data.frame() %>%
  select(c(1:2), any_of(all.adhd.meds$drug)) %>%
  # mutate(methylphenidate = ifelse((methylphenidate+ritalin+concerta)>=1,1,0),
  mutate(methylphenidate = ifelse((methylphenidate+ritalin)>=1 & concerta == 0,1,0),
         amphetamine = ifelse((adderall+amphetamine)>=1,1,0),
         lisdexamfetamine = ifelse((vyvanse+lisdexamfetamine)>=1,1,0),
         stim = ifelse((methylphenidate+adderall+amphetamine+concerta+ritalin+vyvanse+lisdexamfetamine)>=1,1,0),
         non_stim = ifelse((intuniv+strattera+tenex+atomoxetine+clonidine+guanfacine)>=1,1,0),
         # guanfacine = ifelse((guanfacine+tenex+intuniv)>=1,1,0),
         guanfacine = ifelse((guanfacine+tenex)>=1 & intuniv == 0,1,0),
         atomoxetine = ifelse((atomoxetine+strattera)>=1,1,0),
         neither = ifelse(stim >= 1 | non_stim >= 1, 0, 1)) %>% 
  select(-c(ritalin, adderall, tenex, strattera, vyvanse))
####
# ABCD PGS file -----------------------------------------------------------
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))
####
# ABCD predicted MPH response ---------------------------------------------
abcd.pred <- read_rds("../data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(predicted = m) %>%
  # filter(IID %in% sst.r1.deltas$IID) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
####
# nihtbx task ---------------------------------------------------------------
abcd.nihtbx <- read_csv(paste0(abcd.raw.dir, "/neurocognition/nc_y_nihtb.csv"))
# keep variables of interest 
# only keeping the anxiety/frustration/irritability/happiness statuses pre and post task
abcd.nihtbx.filt <- abcd.nihtbx %>%
  select(IID = src_subject_id, eventname, 
         contains("agecorrected")) %>%
  filter(rowSums(is.na(abcd.nihtbx%>%select(contains("agecorrected"))))<=5)
  # drop_na()
# combine nihtbx, age, sex, meds data 
abcd.nihtbx.filt <- inner_join(inner_join(demo, abcd.nihtbx.filt), abcd.meds) %>% 
  # keep this order for OCD pref
  select(IID, eventname, interview_age,sex, 
         starts_with("nih"), 
         colnames(abcd.meds)) %>%
  drop_na()
####
# nihtbx performance based on PGS and categorized by drug -----------------
pgs.nihtbx.drug <- inner_join(abcd.nihtbx.filt, abcd.pgs)
# # Create a custom function for the linear fit
# fit_lm <- function(data) {
#   lm(task_score ~ pg_score, data = data)
# }

pgs.nihtbx.drug %>% 
  # mutate(adq = as.factor(ntile(`ADHD-Demontis`, 4))) %>%
  # mutate(cogq = as.factor(ntile(cog_gFactor, 4))) %>%
  mutate(adq = cut(`ADHD-Demontis`, breaks = quantile(`ADHD-Demontis`, probs = seq(0, 1, length.out = 5)))) %>%
  mutate(cogq = cut(cog_gFactor, breaks = quantile(cog_gFactor, probs = seq(0, 1, length.out = 5)))) %>%
  pivot_longer(cols = starts_with("nih"), names_to = "task", values_to = "task_score") %>%
  pivot_longer(cols = colnames(abcd.pgs)[-1], names_to = "pgs", values_to = "pg_score") %>%
  pivot_longer(cols = adhd.meds$drug, names_to = "drug", values_to = "drug_status") %>%
  pivot_longer(cols = c(adq, cogq), names_to = "quantile", values_to = "q") %>%
  mutate(drug = ifelse(drug_status == 0, "noth", drug)) %>%
  filter(drug %in% c("stim", "non_stim", "neither", "methylphenidate", "guanfacine")) %>%
  filter(grepl("ADHD", pgs) | grepl("gFa", pgs)) %>%
  mutate(task =sub("_agecorrected", "", task)) %>% mutate(task =sub("nihtbx_", "", task)) %>%
  # ggplot(aes(x=pg_score, y=task_score, color = drug, group = interaction(q,drug))) +
  ggplot(aes(x=pg_score, y=task_score, color = drug)) +
  # geom_point(size=0.1)+
  geom_smooth(method = "gam", se = F, na.rm = T) +
  # geom_smooth(method = segmented(fit_lm, psi=0), se = FALSE, aes(group = q)) +
  # geom_step(aes(xend = pg_score), alpha = 0.5) +
  stat_cor(show.legend = F, size = 3, na.rm = T)+
  # facet_wrap("adq", scales = "free") +
  facet_grid2(rows = vars(task), cols = vars(pgs)) +
  labs(x="PGS", caption = paste0("n(samples):\n",
                                 "\tstim: ", sum(pgs.nihtbx.drug$stim), "\n",
                                 "\tnon_stim: ", sum(pgs.nihtbx.drug$non_stim), "\n",
                                 "\tneither: ", sum(pgs.nihtbx.drug$neither)))

pgs.nihtbx.drug %>% 
  mutate(adq = as.factor(ntile(`ADHD-Demontis`, 4))) %>%
  mutate(cogq = as.factor(ntile(cog_gFactor, 4))) %>%
  # mutate(adq = cut(`ADHD-Demontis`, breaks = quantile(`ADHD-Demontis`, probs = seq(0, 1, length.out = 5)))) %>%
  # mutate(cogq = cut(cog_gFactor, breaks = quantile(cog_gFactor, probs = seq(0, 1, length.out = 5)))) %>%
  pivot_longer(cols = starts_with("nih"), names_to = "task", values_to = "task_score") %>%
  pivot_longer(cols = colnames(abcd.pgs)[-1], names_to = "pgs", values_to = "pg_score") %>%
  pivot_longer(cols = adhd.meds$drug, names_to = "drug", values_to = "drug_status") %>%
  pivot_longer(cols = c(adq, cogq), names_to = "quantile", values_to = "q") %>%
  mutate(drug = ifelse(drug_status == 0, "noth", drug)) %>%
  filter(drug %in% c("stim", "non_stim", "neither", "methylphenidate", "guanfacine")) %>%
  filter(grepl("ADHD", pgs) | grepl("gFa", pgs)) %>%
  mutate(task =sub("_agecorrected", "", task)) %>% mutate(task =sub("nihtbx_", "", task)) %>%
  ggplot(aes(x=as.factor(q), y=task_score, color = drug, group = interaction(drug,q))) +
  geom_boxplot(outlier.size = 0.3)+
  # stat_compare_means(label = "p.format", show.legend = F, aes(group=drug)) +
  # stat_compare_means(label = "p.signif") +
  facet_grid2(rows = vars(task), cols = vars(pgs))
  
  
  
  



pgs.nihtbx.drug %>% 
  mutate(adq = as.factor(ntile(`ADHD-Demontis`, 4))) %>%
  mutate(cogq = as.factor(ntile(cog_gFactor, 4))) %>%
  pivot_longer(cols = starts_with("nih"), names_to = "task", values_to = "task_score") %>%
  # pivot_longer(cols = colnames(abcd.pgs)[-1], names_to = "pgs", values_to = "pg_score") %>%
  pivot_longer(cols = adhd.meds$drug, names_to = "drug", values_to = "drug_status") %>%
  # pivot_longer(cols = c(adq, cogq), names_to = "quantile", values_to = "q") %>%
  mutate(drug = ifelse(drug_status == 0, "noth", drug)) %>%
  filter(drug %in% c("stim", "non_stim", "neither")) %>%
  # filter(grepl("ADHD", pgs) | grepl("gFa", pgs)) %>%
  mutate(task =sub("_agecorrected", "", task)) %>% mutate(task =sub("nihtbx_", "", task)) %>%
  ggplot(aes(x=`ADHD-Demontis`, y=cog_gFactor, color = task_score)) +
  geom_point(size=0.1)+
  # geom_smooth(method = "glm", se = F, na.rm = T) +
  # stat_cor(show.legend = F, size = 3, na.rm = T)+
  facet_grid2(cols = vars(drug), rows = vars(task))+
  scale_color_gradient(low = redblu.col[2], high = redblu.col[1])

####

####

####



####
