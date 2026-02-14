################################################################################
#                  check CBCL change by status of taking ADHD drugs            #
################################################################################
# packages setup ----------------------------------------------------------
rm(list = ls());gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(ggpubr);library(ggExtra);library(ggh4x);library(pwr);library(patchwork);library(ggpmisc)
################################################################################
################################################################################
## project dir setup
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response"
setwd(project.dir)
################################################################################
################################################################################
## ABCD demo and directory
abcd.raw.dir <- "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core"
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
## ABCD CBCL data
abcd.cbcl <- read_csv(paste0(abcd.raw.dir, "/mental-health/mh_p_cbcl.csv"))
# keep variables of interest (syndrome subscales and dsm5)
# only keeping the raw scores "_r"
abcd.cbcl.filt <- abcd.cbcl %>%
  rename_at(.vars = vars(starts_with("cbcl_scr_syn_")), .funs = function(x){
    sub("cbcl_scr_syn_", "syn_raw_", sub("_r$","",x))}) %>%
  rename_at(.vars = vars(starts_with("cbcl_scr_dsm5_")), .funs = function(x){
    sub("cbcl_scr_dsm5_", "dsm5_raw_", sub("_r$","",x))}) %>%
  select(IID = src_subject_id, eventname, 
         starts_with(c("syn_raw","dsm5_raw"))&!ends_with(c("_nm","_m","_t"))) %>%
  drop_na()
colnames(abcd.cbcl.filt)

## combine cbcl, age, sex, meds data 
abcd.cbcl.filt <- inner_join(inner_join(demo, abcd.cbcl.filt), abcd.mph) %>% 
  # keep this order for OCD pref
  select(IID, eventname, interview_age,sex, 
         starts_with("syn"), starts_with("dsm5"), 
         colnames(abcd.mph)) %>%
  drop_na()
colnames(abcd.cbcl.filt)

################################################################################
################################################################################
## CBCL correction for age and sex
cbcl.as.corrected <- cbind(abcd.cbcl.filt %>% 
                             select(IID, interview_age, sex, eventname), 
                           apply(abcd.cbcl.filt %>%  select(starts_with(c("syn","dsm5"))), 
                                 MARGIN = 2, FUN = function(x) {
                                   residuals(glm(y ~ interview_age + sex, 
                                                 data = abcd.cbcl.filt %>% mutate(y = x), 
                                                 family = poisson()))
                                   }))
# combine CBCL raw with CBCL corrected
cbcl.all <- inner_join(abcd.cbcl.filt,
                       cbcl.as.corrected %>% 
                         rename_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), 
                                   .funs = function(x) sub("raw", "as", x)))

## correction for other meds
abcd.cbcl.filt2 <- inner_join(abcd.cbcl.filt, abcd.other.adhd)
cbcl.asm.corrected <- cbind(abcd.cbcl.filt2 %>% select(IID, interview_age, sex, eventname), 
                            apply(abcd.cbcl.filt2 %>% select(starts_with(c("syn","dsm5"))), 
                                  MARGIN = 2, FUN = function(x) {
                                    residuals(glm(y ~ ., family = poisson(),
                                                  data = abcd.cbcl.filt2 %>% 
                                                    select(interview_age,sex,other.adhd.meds) %>%
                                                    mutate(y = x)))
                                    }))
cbcl.all2 <- inner_join(abcd.cbcl.filt2,
                        cbcl.asm.corrected %>% 
                          rename_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), 
                                    .funs = function(x) sub("raw", "asm", x)))


################################################################################
################################################################################
## CBCL deltas by MPH status for participants
colnames(cbcl.all)

cbcl.mph.all.deltas <- cbcl.all %>% mutate(MPH = as.numeric(MPH>0)) %>%
  group_by(IID) %>% filter(length(unique(MPH))==2) %>% ungroup() %>%
  group_by(IID, MPH) %>% 
  summarise_at(.vars = vars(starts_with(c("syn","dsm5"))), .funs = function(x) mean(x)) %>%
  ungroup() %>% distinct(IID, MPH, .keep_all = T) %>%
  pivot_longer(cols = starts_with(c("syn","dsm5"))) %>%
  arrange(IID, name, MPH) %>% mutate(MPH = as.factor(MPH)) %>%
  pivot_wider(names_from = MPH, values_from = value, id_cols = c(IID, name)) %>%
  mutate(delta = `1` - `0`)
cbcl.mph.all.deltas2 <- cbcl.all2 %>% mutate(MPH = as.numeric(MPH>0)) %>%
  group_by(IID) %>% filter(length(unique(MPH))==2) %>% ungroup() %>%
  group_by(IID, MPH) %>% 
  summarise_at(.vars = vars(starts_with(c("syn","dsm5"))), .funs = function(x) mean(x)) %>%
  ungroup() %>% distinct(IID, MPH, .keep_all = T) %>%
  pivot_longer(cols = starts_with(c("syn","dsm5"))) %>%
  arrange(IID, name, MPH) %>% mutate(MPH = as.factor(MPH)) %>%
  pivot_wider(names_from = MPH, values_from = value, id_cols = c(IID, name)) %>%
  mutate(delta = `1` - `0`)

cbcl.mph.only.deltas <- cbcl.all %>% mutate(MPH = as.numeric(MPH_only>0)) %>%
  group_by(IID) %>% filter(length(unique(MPH))==2) %>% ungroup() %>%
  group_by(IID, MPH) %>% 
  summarise_at(.vars = vars(starts_with(c("syn","dsm5"))), .funs = function(x) mean(x)) %>%
  ungroup() %>% distinct(IID, MPH, .keep_all = T) %>%
  pivot_longer(cols = starts_with(c("syn","dsm5"))) %>%
  arrange(IID, name, MPH) %>% mutate(MPH = as.factor(MPH)) %>%
  pivot_wider(names_from = MPH, values_from = value, id_cols = c(IID, name)) %>%
  mutate(delta = `1` - `0`)
cbcl.mph.only.deltas2 <- cbcl.all2 %>% mutate(MPH = as.numeric(MPH_only>0)) %>%
  group_by(IID) %>% filter(length(unique(MPH))==2) %>% ungroup() %>%
  group_by(IID, MPH) %>% 
  summarise_at(.vars = vars(starts_with(c("syn","dsm5"))), .funs = function(x) mean(x)) %>%
  ungroup() %>% distinct(IID, MPH, .keep_all = T) %>%
  pivot_longer(cols = starts_with(c("syn","dsm5"))) %>%
  arrange(IID, name, MPH) %>% mutate(MPH = as.factor(MPH)) %>%
  pivot_wider(names_from = MPH, values_from = value, id_cols = c(IID, name)) %>%
  mutate(delta = `1` - `0`)

################################################################################
################################################################################
## deltas by PGS
abcd.pgs %>% pivot_longer(-1, names_to = "pgs_name",values_to = "PGS") %>%
  inner_join(cbcl.mph.all.deltas %>% select(1:2,delta)%>%filter(grepl("adhd|attention",name))) %>%
  filter(grepl("reaction|gFac|ADHD",pgs_name))%>%
  mutate(name = case_when(grepl("_as_",name)~paste0(sub("_as_"," ",name)," \n(corrected for age and sex)"),
                          grepl("_raw_",name)~paste0(sub("_raw_"," ",name)," \n(RAW)")),
         pgs_name = str_replace_all(sub("-"," ",sub("cog_","",pgs_name)),"_"," "),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age and sex"))%>%
  ggplot(aes(PGS,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),cols=vars(pgs_name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH/alternatives - off MPH/alternatives)')
ggsave2("../2026/figs/abcd-mph-and-alt-effectiveness-by-cbcl-on-off-average-delta-cog-and-ADHD-PGS.png",7,9)
abcd.pgs %>% pivot_longer(-1, names_to = "pgs_name",values_to = "PGS") %>%
  inner_join(cbcl.mph.only.deltas %>% select(1:2,delta)%>%filter(grepl("adhd|attention",name))) %>%
  filter(grepl("reaction|gFac|ADHD",pgs_name))%>%
  mutate(name = case_when(grepl("_as_",name)~paste0(sub("_as_"," ",name)," \n(corrected for age and sex)"),
                          grepl("_raw_",name)~paste0(sub("_raw_"," ",name)," \n(RAW)")),
         pgs_name = str_replace_all(sub("-"," ",sub("cog_","",pgs_name)),"_"," "),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age and sex"))%>%
  ggplot(aes(PGS,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),cols=vars(pgs_name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH - off MPH)')
ggsave2("../2026/figs/abcd-mph-only-effectiveness-by-cbcl-on-off-average-delta-cog-and-ADHD-PGS.png",7,9)


abcd.pgs %>% pivot_longer(-1, names_to = "pgs_name",values_to = "PGS") %>%
  inner_join(cbcl.mph.all.deltas2 %>% select(1:2,delta)%>%filter(grepl("adhd|attention",name))) %>%
  filter(grepl("reaction|gFac|ADHD",pgs_name))%>%
  mutate(name = case_when(grepl("_asm_",name)~paste0(sub("_asm_"," ",name)),
                          grepl("_raw_",name)~paste0(sub("_raw_"," ",name)," \n(RAW)")),
         pgs_name = str_replace_all(sub("-"," ",sub("cog_","",pgs_name)),"_"," "),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age, sex, and other ADHD meds"))%>%
  ggplot(aes(PGS,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),cols=vars(pgs_name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH/alternatives - off MPH/alternatives)')
ggsave2("../2026/figs/abcd-mph-and-alt-effectiveness-by-cbcl-on-off-average-delta-cog-and-ADHD-PGS-ADHD-meds-corrected.png",7,9)
abcd.pgs %>% pivot_longer(-1, names_to = "pgs_name",values_to = "PGS") %>%
  inner_join(cbcl.mph.only.deltas2 %>% select(1:2,delta)%>%filter(grepl("adhd|attention",name))) %>%
  filter(grepl("reaction|gFac|ADHD",pgs_name))%>%
  mutate(name = case_when(grepl("_asm_",name)~paste0(sub("_asm_"," ",name)),
                          grepl("_raw_",name)~paste0(sub("_raw_"," ",name)," \n(RAW)")),
         pgs_name = str_replace_all(sub("-"," ",sub("cog_","",pgs_name)),"_"," "),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age, sex, and other ADHD meds"))%>%
  ggplot(aes(PGS,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),cols=vars(pgs_name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH - off MPH)')
ggsave2("../2026/figs/abcd-mph-only-effectiveness-by-cbcl-on-off-average-delta-cog-and-ADHD-PGS-ADHD-meds-corrected.png",7,9)


################################################################################
################################################################################
## deltas to predicted
abcd.pred %>% 
  inner_join(cbcl.mph.all.deltas %>% select(1:2,delta)%>%filter(grepl("adhd|attention",name))) %>%
  mutate(name = case_when(grepl("_as_",name)~paste0(sub("_as_"," ",name)," \n(corrected for age and sex)"),
                          grepl("_raw_",name)~paste0(sub("_raw_"," ",name)," \n(RAW)")),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age and sex"))%>%
  ggplot(aes(predicted,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH/alternatives - off MPH/alternatives)',
                 x="Predicted MPH response")
ggsave2("../2026/figs/abcd-mph-and-alt-effectiveness-by-model-pred.png",4,9)
abcd.pred %>% 
  inner_join(cbcl.mph.only.deltas %>% select(1:2,delta)%>%filter(grepl("adhd|attention",name))) %>%
  mutate(name = case_when(grepl("_as_",name)~paste0(sub("_as_"," ",name)," \n(corrected for age and sex)"),
                          grepl("_raw_",name)~paste0(sub("_raw_"," ",name)," \n(RAW)")),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age and sex"))%>%
  ggplot(aes(predicted,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH - off MPH)',
                 x="Predicted MPH response")
ggsave2("../2026/figs/abcd-mph-only-effectiveness-by-model-pred.png",4,9)


abcd.pred %>% 
  inner_join(cbcl.mph.all.deltas2 %>% select(1:2,delta)%>%filter(grepl("adhd|attention",name))) %>%
  mutate(name = case_when(grepl("_asm_",name)~paste0(sub("_asm_"," ",name)),
                          grepl("_raw_",name)~paste0(sub("_raw_"," ",name)," \n(RAW)")),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age, sex, and other ADHD meds"))%>%
  ggplot(aes(predicted,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH/alternatives - off MPH/alternatives)',
                 x="Predicted MPH response")
ggsave2("../2026/figs/abcd-mph-and-alt-effectiveness-by-model-pred-ADHD-meds-corrected.png",4,9)
abcd.pred %>% 
  inner_join(cbcl.mph.only.deltas2 %>% select(1:2,delta)%>%filter(grepl("adhd|attention",name))) %>%
  mutate(name = case_when(grepl("_asm_",name)~paste0(sub("_asm_"," ",name)),
                          grepl("_raw_",name)~paste0(sub("_raw_"," ",name)," \n(RAW)")),
         tt=ifelse(grepl("RAW",name),"RAW", "corrected for age, sex, and other ADHD meds"))%>%
  ggplot(aes(predicted,delta))+
  geom_point(shape=1,alpha=0.3)+geom_smooth(method="lm",se=F,color=redblu.col.2[2])+ci_ribbon1+
  ggh4x::facet_nested(rows = vars(tt,name),scales="free")+
  ggpubr::stat_cor(label.y.npc = 0.95,color="red")+
  bw.theme+ labs(y='Delta (on MPH - off MPH)',
                 x="Predicted MPH response")
ggsave2("../2026/figs/abcd-mph-only-effectiveness-by-model-pred-ADHD-meds-corrected.png",4,9)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

