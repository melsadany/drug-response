################################################################################
#             distinguishing different drugs based on beh and performance      #
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
                                 "stim", "non_stim",
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
         atomoxetine = ifelse((atomoxetine+strattera)>=1,1,0)) %>% 
  select(-c(ritalin, adderall, tenex, strattera, vyvanse))
####
# ABCD PGS file -----------------------------------------------------------
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))
####
# correlation between PGS and drug type -----------------------------------
pgs.med <- right_join(abcd.pgs, abcd.meds)
pgs.med %>%
  pivot_longer(cols = adhd.meds$drug, names_to = "drug", values_to = "status") %>%
  filter(status >=1) %>%
  pivot_longer(cols = colnames(abcd.pgs)[-1], names_to = "PGS", values_to = "score") %>%
  ggplot(aes(x=score, y=drug)) +
  geom_boxplot(outlier.size = 0.3)+
  facet_wrap("PGS") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  labs(title = "distribution of PGS by drug/drug-category",
       x="", y="")
####
# correlation between age and drug type -----------------------------------
demo.med <- right_join(demo, abcd.meds)
demo.med %>%
  pivot_longer(cols = adhd.meds$drug, names_to = "drug", values_to = "status") %>%
  filter(status >=1) %>%
  pivot_longer(cols = colnames(demo)[3], names_to = "demo", values_to = "score") %>%
  ggplot(aes(x=score, y=drug)) +
  # geom_boxplot(outlier.size = 0.3)+
  geom_density_ridges() +
  labs(title = "distribution of age by drug/drug-category",
       x="", y="")
####
# ABCD CBCL data ----------------------------------------------------------
abcd.cbcl <- read_csv(paste0(abcd.raw.dir, "/mental-health/mh_p_cbcl.csv"))
# keep variables of interest (syndrome subscales and dsm5)
# only keeping the raw scores "_r"
abcd.cbcl.filt <- abcd.cbcl %>%
  select(IID = src_subject_id, eventname, 
         syn_raw_attention = cbcl_scr_syn_attention_r,
         syn_raw_anxdep = cbcl_scr_syn_anxdep_r,
         syn_raw_withdep = cbcl_scr_syn_withdep_r,
         syn_raw_somatic = cbcl_scr_syn_somatic_r,
         syn_raw_social = cbcl_scr_syn_social_r,
         syn_raw_thought = cbcl_scr_syn_thought_r,
         syn_raw_rulebreak = cbcl_scr_syn_rulebreak_r,
         syn_raw_aggressive = cbcl_scr_syn_aggressive_r,
         syn_raw_internal = cbcl_scr_syn_internal_r,
         syn_raw_external = cbcl_scr_syn_external_r,
         syn_raw_totprob = cbcl_scr_syn_totprob_r,
         dsm5_raw_depress = cbcl_scr_dsm5_depress_r,
         dsm5_raw_anxdisord = cbcl_scr_dsm5_anxdisord_r,
         dsm5_raw_somaticpr = cbcl_scr_dsm5_somaticpr_r,
         dsm5_raw_adhd = cbcl_scr_dsm5_adhd_r,
         dsm5_raw_opposit = cbcl_scr_dsm5_opposit_r,
         dsm5_raw_conduct = cbcl_scr_dsm5_conduct_r) %>%
  drop_na()
# combine cbcl, age, sex, meds data 
abcd.cbcl.filt <- inner_join(inner_join(demo, abcd.cbcl.filt), abcd.meds) %>% 
  # keep this order for OCD pref
  select(IID, eventname, interview_age,sex, 
         starts_with("syn"), starts_with("dsm5"), 
         colnames(abcd.meds)) %>%
  drop_na()
####
# CBCL correction for age and sex -----------------------------------------
cbcl.as.corrected <- cbind(abcd.cbcl.filt %>% 
                             select(IID, interview_age, sex, eventname), 
                           apply(abcd.cbcl.filt %>% 
                                   select(starts_with("syn"), 
                                          starts_with("dsm5")), 
                                 MARGIN = 2, FUN = function(x) {
                                   residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                 data = abcd.cbcl.filt %>% mutate(y = x), 
                                                 family = poisson()))
                                 }))
# combine CBCL raw with CBCL corrected
cbcl.all <- inner_join(abcd.cbcl.filt,
                       cbcl.as.corrected %>% 
                         rename_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), 
                                   .funs = function(x) sub("raw", "as", x)))
####
# correlation between CBCL and drug type -----------------------------------
cbcl.all %>%
  pivot_longer(cols = adhd.meds$drug, names_to = "drug", values_to = "status") %>%
  filter(status >=1) %>%
  pivot_longer(cols = c(contains("as"),-contains("totprob")), names_to = "CBCL", values_to = "score") %>%
  ggplot(aes(x=score, y=drug)) +
  geom_boxplot(outlier.size = 0.3)+
  facet_wrap("CBCL", scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  labs(title = "distribution of CBCL by drug/drug-category",
       x="", y="")
####
# ABCD SST data -----------------------------------------------------------
sst.raw <- read_csv(paste0(abcd.raw.dir, "/imaging/mri_y_tfmr_sst_beh.csv"))
question.dict <- data.frame(q0 = colnames(sst.raw), description = t(sst.raw)[,1], row.names = NULL)
# only keeping 4 questions of interest for the first run of the task only
sst.r1 <- left_join(sst.raw[-1,] %>%
                      select(IID = src_subject_id, eventname,
                             e_raw_correct_go = tfmri_sst_r1_beh_crgo_nt,
                             e_raw_correct_stop = tfmri_sst_r1_beh_crs_nt,
                             e_raw_stop_doesnot_stop = tfmri_sst_r1_beh_ssds_nt,
                             e_raw_no_response_on_go = tfmri_sst_r1_beh_nrgo_nt,
                             tfmri_sst_beh_switchflag) %>%
                      filter(is.na(tfmri_sst_beh_switchflag) == F),
                    left_join(abcd.meds, demo)) %>% drop_na() %>%
  mutate_at(.vars = vars(3:7, interview_age), .funs = function(x) as.numeric(x))
# combine sst, age, sex, meds data 
abcd.sst.filt <- inner_join(inner_join(demo, sst.r1), abcd.meds) %>% 
  # keep this order for OCD pref
  select(IID, eventname, interview_age,sex, 
         starts_with("e_"), 
         colnames(abcd.meds)) %>%
  drop_na()
####
# SST correction for age and sex -----------------------------------------
sst.as.corrected <- cbind(abcd.sst.filt %>% 
                            select(IID, interview_age, sex, eventname), 
                          apply(abcd.sst.filt %>% 
                                  select(starts_with("e_")), 
                                MARGIN = 2, FUN = function(x) {
                                  residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                data = abcd.sst.filt %>% mutate(y = x), 
                                                family = poisson()))
                                }))
# combine sst raw with sst corrected
sst.all <- inner_join(abcd.sst.filt,
                      sst.as.corrected %>% 
                        rename_at(.vars = vars(starts_with("e_")), 
                                  .funs = function(x) sub("raw", "as", x)))
####
# correlation between SST and drug type -----------------------------------
sst.all %>%
  pivot_longer(cols = adhd.meds$drug, names_to = "drug", values_to = "status") %>%
  filter(status >=1) %>%
  pivot_longer(cols = contains("as"), names_to = "SST", values_to = "score") %>%
  ggplot(aes(x=score, y=drug)) +
  geom_boxplot(outlier.size = 0.3)+
  facet_wrap("SST", scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  labs(title = "distribution of SST by drug/drug-category",
       x="", y="")
####

####
