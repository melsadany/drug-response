################################################################################
#                  check CBCL change by status of taking ADHD drugs            #
################################################################################
rm(list = ls());gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
################################################################################
## project dir setup
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/2026")
setwd(project.dir)
################################################################################
################################################################################
## ABCD demo and directory
abcd.raw.dir <- correct_path("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core")
age <- read_csv(paste0(abcd.raw.dir, "/abcd-general/abcd_y_lt.csv")) %>% 
  select(IID = src_subject_id, eventname, interview_age)
sex <- read_csv(paste0(abcd.raw.dir, "/gender-identity-sexual-health/gish_p_gi.csv")) %>%
  mutate(sex = ifelse(demo_sex_v2 == 2, "Female", ifelse(demo_sex_v2 == 1, "Male", ifelse(demo_sex_v2 == 3, "intersex_M", NA)))) %>%
  select(IID = src_subject_id, sex) %>%
  distinct(IID, .keep_all = T)
demo <- full_join(age, sex)
################################################################################
################################################################################
## ABCD meds
mph.equv <- c("methylphenidate","ritalin","concerta", "methylin", "aptensio", 
              "daytrana","quillivant", "quillichew", "jornay", "adhansia", 
              "metadate","relexxii") # https://www.drugs.com/medical-answers/brands-methylphenidate-3510739/
other.adhd.meds <- c("adderall", "vyvanse", "dextroamphetamine","intuniv", 
                     "strattera","tenex", "amphetamine", "dexmethylphenidate", 
                     "lisdexamfetamine","atomoxetine", "clonidine", "guanfacine")
abcd.other.adhd <- read_rds(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/abcd5-meds-matrix.rds")) %>%
  as.data.frame() %>% select(c(1:2), any_of(other.adhd.meds))
abcd.mph <- read_rds(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/abcd5-meds-matrix.rds")) %>%
  as.data.frame() %>% select(c(1:2), any_of(mph.equv)) %>% ungroup() %>%
  mutate(MPH_only = methylphenidate) %>%
  pivot_longer(cols = -c(1:2, MPH_only)) %>% group_by(IID, eventname, MPH_only) %>% 
  summarise(MPH = sum(value))

length(unique(abcd.mph$IID[abcd.mph$MPH==0])) #11693
length(unique(abcd.mph$IID[abcd.mph$MPH>0])) # 823
length(unique(abcd.mph$IID[abcd.mph$MPH==1])) #789
length(unique(abcd.mph$IID[abcd.mph$MPH==2])) #160
length(unique(abcd.mph$IID[abcd.mph$MPH==3])) #7
length(unique(abcd.mph$IID[abcd.mph$MPH==4])) #0

################################################################################
################################################################################
## ABCD PGS file
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB"), contains("EXTRA")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))
################################################################################
################################################################################
## ABCD predicted MPH response
abcd.pred <- read_rds("../data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(predicted = m) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
################################################################################
################################################################################
## ABCD SST data
sst.raw <- read_csv(paste0(abcd.raw.dir, "/imaging/mri_y_tfmr_sst_beh.csv"))
question.dict <- data.frame(q0 = colnames(sst.raw), description = t(sst.raw)[,1], row.names = NULL)
sst.r1 <- left_join(sst.raw[-1,] %>%
                      select(IID = src_subject_id, eventname,
                             e_raw_correct_go = tfmri_sst_r1_beh_crgo_nt,
                             e_raw_correct_stop = tfmri_sst_r1_beh_crs_nt,
                             e_raw_stop_doesnot_stop = tfmri_sst_r1_beh_ssds_nt,
                             e_raw_no_response_on_go = tfmri_sst_r1_beh_nrgo_nt,
                             tfmri_sst_beh_switchflag) %>%
                      filter(is.na(tfmri_sst_beh_switchflag) == F),
                    left_join(abcd.mph, demo)) %>% drop_na() %>%
  mutate_at(.vars = vars(3:7, interview_age), .funs = function(x) as.numeric(x))
colnames(sst.r1)
summary(sst.r1)
################################################################################
################################################################################
## SST correction for age and sex
sst.as.corrected <- cbind(sst.r1 %>% 
                             select(IID, interview_age, sex, eventname), 
                           apply(sst.r1 %>%  select(starts_with("e_raw")), 
                                 MARGIN = 2, FUN = function(x) {
                                   residuals(glm(y ~ interview_age + sex, 
                                                 data = sst.r1 %>% mutate(y = x), 
                                                 family = poisson()))
                                   }))
# combine CBCL raw with CBCL corrected
sst.all <- inner_join(sst.r1,
                       sst.as.corrected %>% 
                         rename_at(.vars = vars(starts_with("e_raw")), 
                                   .funs = function(x) sub("raw", "as", x)))

## correction for other meds
sst.r12 <- inner_join(sst.r1, abcd.other.adhd)
sst.asm.corrected <- cbind(sst.r12 %>% select(IID, interview_age, sex, eventname), 
                            apply(sst.r12 %>% select(starts_with("e_raw")), 
                                  MARGIN = 2, FUN = function(x) {
                                    residuals(glm(y ~ ., family = poisson(),
                                                  data = sst.r12 %>% 
                                                    select(interview_age,sex,other.adhd.meds) %>%
                                                    mutate(y = x)))
                                  }))
sst.all2 <- inner_join(sst.r12,
                        sst.asm.corrected %>% 
                          rename_at(.vars = vars(starts_with("e_raw")), 
                                    .funs = function(x) sub("raw", "asm", x)))
################################################################################
################################################################################
## SST deltas by MPH status for participants
colnames(sst.all)

sst.mph.all.deltas <- sst.all %>% mutate(MPH = as.numeric(MPH>0)) %>%
  group_by(IID) %>% filter(length(unique(MPH))==2) %>% ungroup() %>%
  group_by(IID, MPH) %>% 
  summarise_at(.vars = vars(starts_with(c("e_"))), .funs = function(x) mean(x)) %>%
  ungroup() %>% distinct(IID, MPH, .keep_all = T) %>%
  pivot_longer(cols = starts_with(c("e_"))) %>%
  arrange(IID, name, MPH) %>% mutate(MPH = as.factor(MPH)) %>%
  pivot_wider(names_from = MPH, values_from = value, id_cols = c(IID, name)) %>%
  mutate(delta = `1` - `0`)

sst.mph.only.deltas <- sst.all %>% mutate(MPH = as.numeric(MPH_only>0)) %>%
  group_by(IID) %>% filter(length(unique(MPH))==2) %>% ungroup() %>%
  group_by(IID, MPH) %>% 
  summarise_at(.vars = vars(starts_with(c("e_"))), .funs = function(x) mean(x)) %>%
  ungroup() %>% distinct(IID, MPH, .keep_all = T) %>%
  pivot_longer(cols = starts_with(c("e_"))) %>%
  arrange(IID, name, MPH) %>% mutate(MPH = as.factor(MPH)) %>%
  pivot_wider(names_from = MPH, values_from = value, id_cols = c(IID, name)) %>%
  mutate(delta = `1` - `0`)

sst.mph.all.deltas2 <- sst.all2 %>% mutate(MPH = as.numeric(MPH>0)) %>%
  group_by(IID) %>% filter(length(unique(MPH))==2) %>% ungroup() %>%
  group_by(IID, MPH) %>% 
  summarise_at(.vars = vars(starts_with(c("e_"))), .funs = function(x) mean(x)) %>%
  ungroup() %>% distinct(IID, MPH, .keep_all = T) %>%
  pivot_longer(cols = starts_with(c("e_"))) %>%
  arrange(IID, name, MPH) %>% mutate(MPH = as.factor(MPH)) %>%
  pivot_wider(names_from = MPH, values_from = value, id_cols = c(IID, name)) %>%
  mutate(delta = `1` - `0`)
sst.mph.only.deltas2 <- sst.all2 %>% mutate(MPH = as.numeric(MPH_only>0)) %>%
  group_by(IID) %>% filter(length(unique(MPH))==2) %>% ungroup() %>%
  group_by(IID, MPH) %>% 
  summarise_at(.vars = vars(starts_with(c("e_"))), .funs = function(x) mean(x)) %>%
  ungroup() %>% distinct(IID, MPH, .keep_all = T) %>%
  pivot_longer(cols = starts_with(c("e_"))) %>%
  arrange(IID, name, MPH) %>% mutate(MPH = as.factor(MPH)) %>%
  pivot_wider(names_from = MPH, values_from = value, id_cols = c(IID, name)) %>%
  mutate(delta = `1` - `0`)

################################################################################
################################################################################
## deltas by PGS
abcd.pgs %>% pivot_longer(-1, names_to = "pgs_name",values_to = "PGS") %>%
  inner_join(sst.mph.all.deltas %>% select(1:2,delta)) %>%
  filter(grepl("reaction|gFac|ADHD",pgs_name))%>%
  mutate(name = case_when(grepl("_as_",name)~str_replace_all(paste0(sub("e_as_"," ",name)," \n(corrected for age and sex)"),"_"," "),
                          grepl("_raw_",name)~str_replace_all(paste0(sub("e_raw_"," ",name)," \n(RAW)"),"_"," ")),
         pgs_name = str_replace_all(sub("-"," ",sub("cog_","",pgs_name)),"_"," "),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age and sex"))%>%
  ggplot(aes(PGS,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),cols=vars(pgs_name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH/alternatives - off MPH/alternatives)')
ggsave2("../2026/figs/abcd-mph-and-alt-effectiveness-by-SST-on-off-average-delta-cog-and-ADHD-PGS.png",7,14)
abcd.pgs %>% pivot_longer(-1, names_to = "pgs_name",values_to = "PGS") %>%
  inner_join(sst.mph.only.deltas %>% select(1:2,delta)) %>%
  filter(grepl("reaction|gFac|ADHD",pgs_name))%>%
  mutate(name = case_when(grepl("_as_",name)~str_replace_all(paste0(sub("e_as_"," ",name)," \n(corrected for age and sex)"),"_"," "),
                          grepl("_raw_",name)~str_replace_all(paste0(sub("e_raw_"," ",name)," \n(RAW)"),"_"," ")),
         pgs_name = str_replace_all(sub("-"," ",sub("cog_","",pgs_name)),"_"," "),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age and sex"))%>%
  ggplot(aes(PGS,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),cols=vars(pgs_name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH - off MPH)')
ggsave2("../2026/figs/abcd-mph-only-effectiveness-by-SST-on-off-average-delta-cog-and-ADHD-PGS.png",7,14)


abcd.pgs %>% pivot_longer(-1, names_to = "pgs_name",values_to = "PGS") %>%
  inner_join(sst.mph.all.deltas2 %>% select(1:2,delta)) %>%
  filter(grepl("reaction|gFac|ADHD",pgs_name))%>%
  mutate(name = case_when(grepl("_asm_",name)~str_replace_all(paste0(sub("e_asm_"," ",name)),"_"," "),
                          grepl("_raw_",name)~str_replace_all(paste0(sub("e_raw_"," ",name)," \n(RAW)"),"_"," ")),
         pgs_name = str_replace_all(sub("-"," ",sub("cog_","",pgs_name)),"_"," "),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age, sex, and other ADHD meds"))%>%
  ggplot(aes(PGS,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),cols=vars(pgs_name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH/alternatives - off MPH/alternatives)')
ggsave2("../2026/figs/abcd-mph-and-alt-effectiveness-by-SST-on-off-average-delta-cog-and-ADHD-PGS-ADHD-meds-corrected.png",7,14)
abcd.pgs %>% pivot_longer(-1, names_to = "pgs_name",values_to = "PGS") %>%
  inner_join(sst.mph.only.deltas2 %>% select(1:2,delta)) %>%
  filter(grepl("reaction|gFac|ADHD",pgs_name))%>%
  mutate(name = case_when(grepl("_asm_",name)~str_replace_all(paste0(sub("e_asm_"," ",name)),"_"," "),
                          grepl("_raw_",name)~str_replace_all(paste0(sub("e_raw_"," ",name)," \n(RAW)"),"_"," ")),
         pgs_name = str_replace_all(sub("-"," ",sub("cog_","",pgs_name)),"_"," "),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age, sex, and other ADHD meds"))%>%
  ggplot(aes(PGS,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),cols=vars(pgs_name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH - off MPH)')
ggsave2("../2026/figs/abcd-mph-only-effectiveness-by-SST-on-off-average-delta-cog-and-ADHD-PGS-ADHD-meds-corrected.png",7,14)

################################################################################
################################################################################
## deltas to predicted
abcd.pred %>% 
  inner_join(sst.mph.all.deltas %>% select(1:2,delta)) %>%
  mutate(name = case_when(grepl("_as_",name)~str_replace_all(paste0(sub("e_as_"," ",name)," \n(corrected for age and sex)"),"_"," "),
                          grepl("_raw_",name)~str_replace_all(paste0(sub("e_raw_"," ",name)," \n(RAW)"),"_"," ")),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age and sex"))%>%
  ggplot(aes(predicted,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH/alternatives - off MPH/alternatives)',
                 x="Predicted MPH response")
ggsave2("../2026/figs/abcd-mph-and-alt-SST-effectiveness-by-model-pred.png",4,14)
abcd.pred %>% 
  inner_join(sst.mph.only.deltas %>% select(1:2,delta)) %>%
  mutate(name = case_when(grepl("_as_",name)~str_replace_all(paste0(sub("e_as_"," ",name)," \n(corrected for age and sex)"),"_"," "),
                          grepl("_raw_",name)~str_replace_all(paste0(sub("e_raw_"," ",name)," \n(RAW)"),"_"," ")),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age and sex"))%>%
  ggplot(aes(predicted,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH - off MPH)',
                 x="Predicted MPH response")
ggsave2("../2026/figs/abcd-mph-only-SST-effectiveness-by-model-pred.png",4,16)


abcd.pred %>% 
  inner_join(sst.mph.all.deltas2 %>% select(1:2,delta)) %>%
  mutate(name = case_when(grepl("_asm_",name)~str_replace_all(paste0(sub("e_asm_"," ",name)),"_"," "),
                          grepl("_raw_",name)~str_replace_all(paste0(sub("e_raw_"," ",name)," \n(RAW)"),"_"," ")),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age, sex, and other ADHD meds"))%>%
  ggplot(aes(predicted,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH/alternatives - off MPH/alternatives)',
                 x="Predicted MPH response")
ggsave2("../2026/figs/abcd-mph-and-alt-SST-effectiveness-by-model-pred-ADHD-meds-corrected.png",4,14)
abcd.pred %>% 
  inner_join(sst.mph.only.deltas2 %>% select(1:2,delta)) %>%
  mutate(name = case_when(grepl("_asm_",name)~str_replace_all(paste0(sub("e_asm_"," ",name)),"_"," "),
                          grepl("_raw_",name)~str_replace_all(paste0(sub("e_raw_"," ",name)," \n(RAW)"),"_"," ")),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age, sex, and other ADHD meds"))%>%
  ggplot(aes(predicted,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH - off MPH)',
                 x="Predicted MPH response")
ggsave2("../2026/figs/abcd-mph-only-SST-effectiveness-by-model-pred-ADHD-meds-corrected.png",4,16)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

