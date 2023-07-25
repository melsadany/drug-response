################################################################################
#                   predcit drug response for ABCD participants                #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(ggpubr);library(ggExtra);library(ggh4x)
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response"
setwd(project.dir)
################################################################################
det.cap <- paste0("delta per question = score_on_MPH - score_off_MPH", "\n", 
                  "corrected for:", "\n",
                  "\tinterview_age + sex + interview_age:sex","\n",
                  "\tclonidine + adderall + concerta + ", "\n",
                  "\tvyvanse + dextroamphetamine +", "\n",
                  "\tritalin + intuniv +tenex + ", "\n", 
                  "\tamphetamine+ guanfacine +","\n",
                  "\tlisdexamfetamine + dexmethylphenidate")
################################################################################
abcd.pred <- read_rds("../data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(predicted = m) %>%
  # filter(IID %in% sst.r1.deltas$IID) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))
#########################
abcd.raw.dir <- "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core"
age <- read_csv(paste0(abcd.raw.dir, "/abcd-general/abcd_y_lt.csv")) %>% 
  select(IID = src_subject_id, eventname, interview_age)
sex <- read_csv(paste0(abcd.raw.dir, "/gender-identity-sexual-health/gish_p_gi.csv")) %>%
  mutate(sex = ifelse(demo_sex_v2 == 2, "Female", ifelse(demo_sex_v2 == 1, "Male", ifelse(demo_sex_v2 == 3, "intersex_M", NA)))) %>%
  select(IID = src_subject_id, sex) %>%
  distinct(IID, .keep_all = T)
demo <- full_join(age, sex)
#########################
# meds
adhd.meds <- data.frame(drug = c("methylphenidate", "adderall", "concerta", "vyvanse", 
                                 "dextroamphetamine",  "ritalin", "intuniv", "strattera",
                                 "tenex", "amphetamine", "dexmethylphenidate", "lisdexamfetamine",
                                 "atomoxetine", "clonidine", "guanfacine"
                                 ))
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/abcd5-meds-matrix.rds") %>%
  as.data.frame() %>%
  select(c(1:2), any_of(adhd.meds$drug))
#########################
# decided to keep count of answers in trial run 1 (r1) only
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
                    left_join(abcd.meds, demo)) %>% drop_na() %>%
  mutate_at(.vars = vars(3:7, interview_age), .funs = function(x) as.numeric(x))
hist(sst.r1)
#########################
# get participants with 2 or more measurements
sst.r1.tom <- sst.r1 %>%
  group_by(IID) %>%
  filter(n_distinct(eventname) >= 2, n_distinct(methylphenidate) >= 2) %>%
  ungroup()
# correct for age, sex, their interaction, other meds
# age-sex only
sst.r1.tom.as.corrected <- cbind(sst.r1.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                              apply(sst.r1.tom %>% 
                                select(starts_with("e_")), 
                              MARGIN = 2, FUN = function(x) {
                                residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                                           data = sst.r1.tom %>% mutate(y = x), 
                                                                           family = poisson()))
                                }))
# age, sex, interaction, and other ADHD meds
sst.r1.tom.asmeds.corrected <- cbind(sst.r1.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                                 apply(sst.r1.tom %>% 
                                         select(starts_with("e_")), 
                                       MARGIN = 2, FUN = function(x) {
                                         residuals(glm(y ~ interview_age + sex + interview_age:sex
                                                       + clonidine + adderall + concerta + vyvanse + ritalin +
                                                         intuniv + tenex + guanfacine + dexmethylphenidate
                                                       , 
                                                       data = sst.r1.tom %>% mutate(y = x), 
                                                       family = poisson()))
                                       }))
# combine raw and corrected sst data
sst.r1.tom.all <- inner_join(sst.r1.tom,
                             inner_join(sst.r1.tom.as.corrected %>% rename_at(.vars = vars(starts_with("e_")), 
                                                                              .funs = function(x) sub("e_raw_", "e_as_", x)),
                                        sst.r1.tom.asmeds.corrected %>% rename_at(.vars = vars(starts_with("e_")), 
                                                                                  .funs = function(x) sub("e_raw_", "e_asmeds_", x)))) %>% 
  group_by(IID) %>%
  filter(n_distinct(methylphenidate)==2) %>%
  group_by(IID, methylphenidate) %>%
  mutate_at(.vars = vars(starts_with("e_")), .funs = function(x) mean(x)) %>%
  ungroup() %>%
  distinct(IID, methylphenidate, .keep_all = T) %>%
  pivot_longer(cols = starts_with("e_"), names_to = "question", values_to = "val") %>%
  arrange(IID, question, methylphenidate) %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  pivot_wider(names_from = methylphenidate, values_from = val, id_cols = c(IID, sex, question)) %>%
  mutate(value_type = factor(ifelse(grepl("e_raw_", question), "raw data", ifelse(grepl("e_as_", question), 
                                                                           "corrected for age, sex, and interaction", 
                                                                           "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  mutate(delta = `1` - `0`) # being on MPH - not being on MPH
###
rm(sst.r1.tom.asmeds.corrected)
rm(sst.r1.tom.as.corrected)
gc()
###
################################################################################
################################## plots #######################################
# scatter plot for before and after score for raw, age&sex corrected, and age&sex&meds corrected
p1.s <- inner_join(sst.r1.tom.all, abcd.pgs) %>%
  mutate(ADHD_decile = as.factor(ntile(`ADHD-Demontis`, 10))) %>%
  mutate(question = as.factor(sub("e_.*?_","",question))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  ggplot(aes(x=`1`, y=`0`, alpha = ADHD_decile))+
  geom_point(position = "jitter", size=1) +
  scale_color_manual(values = redblu.col)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="on - MPH", y = "off - MPH")
# density plot for on/off MPH
p1.d <- sst.r1.tom.all %>%
  mutate(question = as.factor(sub("e_.*?_","",question))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  pivot_longer(cols = c(`0`,`1`), names_to = "MPH_status", values_to = "score") %>%
  ggplot(aes(x=score, color = MPH_status))+
  geom_density(position = "identity") +
  scale_color_manual(values = redblu.col)+
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="")
# paired boxplot
p1.b <- sst.r1.tom.all %>% 
  # filter(value_type == "raw data") %>%
  mutate(question = as.factor(sub("e_.*?_","",question))) %>%
  arrange(IID, question) %>%
  rename(off_MPH = `0`, on_MPH = `1`) %>%
  ggpaired(cond1 = "off_MPH", cond2 = "on_MPH", y = "val", 
           line.color = "gray", line.size = 0.1, palette = "jco", 
           ylab = "", xlab = "", point.size = 0.3) +
  facet_grid2(rows = vars(value_type), cols = vars(question), scales = "free_y", independent = "y")+
  theme_minimal() +
  stat_compare_means(paired = T, size = 3, method = "t.test")
# correlation between predicted response to MPH and participants performance in SST whether data is corrected or not
abcd.c.1 <- inner_join(abcd.pred, sst.r1.tom.all %>% select(-c(`0`,`1`))%>% pivot_wider(names_from = "question", values_from = "delta"))
p1.p <- corr.table(abcd.c.1%>%select(predicted), abcd.c.1 %>% select(starts_with("e_")), method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  mutate(question = as.factor(sub("e_.*?_","",V2))) %>%
  mutate(value_type = factor(ifelse(grepl("e_raw_", V2), "raw data", ifelse(grepl("e_as_", V2), 
                                                                                  "corrected for age, sex, and interaction", 
                                                                                  "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(pval < 0.1, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.1%>%distinct(IID)), "\n",
                                              det.cap))
################################################################################
# cbcl
abcd.cbcl <- read_csv(paste0(abcd.raw.dir, "/mental-health/mh_p_cbcl.csv"))
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

abcd.cbcl.filt <- inner_join(inner_join(abcd.cbcl.filt, 
                            abcd.meds), demo) %>% select(IID, eventname, interview_age,sex, 
                                          starts_with("syn"), starts_with("dsm5"), 
                                          colnames(abcd.meds)) %>%
  drop_na()
#########################
# get participants with 2 or more measurements
cbcl.tom <- abcd.cbcl.filt %>%
  group_by(IID) %>%
  filter(n_distinct(eventname) >= 2, n_distinct(methylphenidate) >= 2) %>%
  ungroup()
# correct for age, sex, their interaction, other meds
# age-sex only
cbcl.tom.as.corrected <- cbind(cbcl.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                                 apply(cbcl.tom %>% 
                                         select(starts_with("syn"), starts_with("dsm5")), 
                                       MARGIN = 2, FUN = function(x) {
                                         residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                       data = cbcl.tom %>% mutate(y = x), 
                                                       family = poisson()))
                                       }))
# age, sex, interaction, and other ADHD meds
cbcl.tom.asmeds.corrected <- cbind(cbcl.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                                     apply(cbcl.tom %>% 
                                             select(starts_with("syn"), starts_with("dsm5")), 
                                           MARGIN = 2, FUN = function(x) {
                                             residuals(glm(y ~ interview_age + sex + interview_age:sex
                                                           + clonidine + adderall + concerta + vyvanse + ritalin +
                                                             intuniv + tenex + guanfacine + dexmethylphenidate
                                                           , 
                                                           data = cbcl.tom %>% mutate(y = x), 
                                                           family = poisson()))
                                           }))
# combine raw and corrected sst data
cbcl.tom.all <- inner_join(cbcl.tom,
                             inner_join(cbcl.tom.as.corrected %>% rename_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), 
                                                                              .funs = function(x) ifelse(grepl("syn", x), sub("syn_raw_", "syn_as_", x),
                                                                                                         sub("dsm5_raw_", "dsm5_as_", x))),
                                        cbcl.tom.asmeds.corrected %>% rename_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), 
                                                                                .funs = function(x) ifelse(grepl("syn", x), sub("syn_raw_", "syn_asmeds_", x),
                                                                                                           sub("dsm5_raw_", "dsm5_asmeds_", x))))) %>% 
  group_by(IID) %>%
  filter(n_distinct(methylphenidate)==2) %>%
  group_by(IID, methylphenidate) %>%
  mutate_at(.vars = vars(starts_with("e_")), .funs = function(x) mean(x)) %>%
  ungroup() %>%
  distinct(IID, methylphenidate, .keep_all = T) %>%
  pivot_longer(cols = c(starts_with("syn"), starts_with("dsm5")), names_to = "question", values_to = "val") %>%
  arrange(IID, question, methylphenidate) %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  pivot_wider(names_from = methylphenidate, values_from = val, id_cols = c(IID, sex, question)) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", question), "raw data", ifelse(grepl("as_", question), 
                                                                                  "corrected for age, sex, and interaction", 
                                                                                  "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  mutate(subscale = ifelse(grepl("dsm5", question), "dsm5", "syn")) %>%
  mutate(delta = `1` - `0`) # being on MPH - not being on MPH
###
rm(cbcl.tom.asmeds.corrected)
rm(cbcl.tom.as.corrected)
gc()
###
################################################################################
################################## plots #######################################
# scatter plot for before and after score for raw, age&sex corrected, and age&sex&meds corrected
p2.s.syn <- inner_join(cbcl.tom.all, abcd.pgs) %>%
  mutate(ADHD_decile = as.factor(ntile(`ADHD-Demontis`, 10))) %>%
  filter(subscale=="syn") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  ggplot(aes(x=`1`, y=`0`, alpha = ADHD_decile))+
  geom_point(position = "jitter", size=1) +
  scale_color_manual(values = redblu.col)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="on - MPH", y = "off - MPH")
p2.s.dsm5 <- inner_join(cbcl.tom.all, abcd.pgs) %>%
  mutate(ADHD_decile = as.factor(ntile(`ADHD-Demontis`, 10))) %>%
  filter(subscale=="dsm5") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  ggplot(aes(x=`1`, y=`0`, alpha = ADHD_decile))+
  geom_point(position = "jitter", size=1) +
  scale_color_manual(values = redblu.col)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="on - MPH", y = "off - MPH")
# density plot for on/off MPH
p2.d.syn <- cbcl.tom.all %>%
  filter(subscale=="syn") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  pivot_longer(cols = c(`0`,`1`), names_to = "MPH_status", values_to = "score") %>%
  ggplot(aes(x=score, color = MPH_status))+
  geom_density(position = "identity") +
  scale_color_manual(values = redblu.col)+
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="")
p2.d.dsm5 <- cbcl.tom.all %>%
  filter(subscale=="dsm5") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  pivot_longer(cols = c(`0`,`1`), names_to = "MPH_status", values_to = "score") %>%
  ggplot(aes(x=score, color = MPH_status))+
  geom_density(position = "identity") +
  scale_color_manual(values = redblu.col)+
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="")
# paired boxplot
p2.b.syn <- cbcl.tom.all %>% 
  filter(subscale=="syn") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  arrange(IID, question) %>%
  rename(off_MPH = `0`, on_MPH = `1`) %>%
  ggpaired(cond1 = "off_MPH", cond2 = "on_MPH", y = "val", 
           line.color = "gray", line.size = 0.1, palette = "jco", 
           ylab = "", xlab = "", point.size = 0.3) +
  facet_grid2(rows = vars(value_type), cols = vars(question), scales = "free_y", independent = "y")+
  theme_minimal() +
  stat_compare_means(paired = T, size = 3, method = "t.test")
p2.b.dsm5 <- cbcl.tom.all %>% 
  filter(subscale=="dsm5") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  arrange(IID, question) %>%
  rename(off_MPH = `0`, on_MPH = `1`) %>%
  ggpaired(cond1 = "off_MPH", cond2 = "on_MPH", y = "val", 
           line.color = "gray", line.size = 0.1, palette = "jco", 
           ylab = "", xlab = "", point.size = 0.3) +
  facet_grid2(rows = vars(value_type), cols = vars(question), scales = "free_y", independent = "y")+
  theme_minimal() +
  stat_compare_means(paired = T, size = 3, method = "t.test")
# correlation between predicted response to MPH and participants performance in SST whether data is corrected or not
abcd.c.2 <- inner_join(abcd.pred, cbcl.tom.all %>% select(-c(`0`,`1`))%>% pivot_wider(names_from = "question", values_from = "delta"))
p2.p <- corr.table(abcd.c.2%>%select(predicted), abcd.c.2 %>% select(starts_with("syn"), starts_with("dsm5")), method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  mutate(question = as.factor(ifelse(grepl("syn", V2), sub("syn_.*?_", "syn_", V2),
                                     sub("dsm5_.*?_", "dsm5_", V2)))) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                            "corrected for age, sex, and interaction", 
                                                                            "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(pval < 0.1, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.2%>%distinct(IID))))
# 
# inner_join(abcd.pred, cbcl.tom.all) %>%
#   filter(subscale=="dsm5") %>%
#   mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
#                                      sub("dsm5_.*?_", "dsm5_", question)))) %>%
#   mutate(value_type = as.factor(value_type)) %>%
#   ggplot(aes(x=`1`, y=`0`, alpha = predicted))+
#   geom_point(position = "jitter") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
#               independent = T) +
#   labs(x="on - MPH", y = "off - MPH")
# 
# 
################################################################################
# brief problems monitor bpm
abcd.bpm <- read_csv(paste0(abcd.raw.dir, "/mental-health/mh_y_bpm.csv")) %>%
  select(IID = src_subject_id, eventname, 
         bpm_raw_att = bpm_y_scr_attention_r, 
         bpm_raw_int = bpm_y_scr_internal_r,
         bpm_raw_ext = bpm_y_scr_external_r, 
         bpm_raw_tot = bpm_y_scr_totalprob_r)
abcd.bpm <- inner_join(inner_join(abcd.bpm, abcd.meds), demo) %>%
  drop_na()
#########################
# get participants with 2 or more measurements
bpm.tom <- abcd.bpm %>%
  group_by(IID) %>%
  filter(n_distinct(eventname) >= 2, n_distinct(methylphenidate) >= 2) %>%
  ungroup()
# correct for age, sex, their interaction, other meds
# age-sex only
bpm.tom.as.corrected <- cbind(bpm.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                                 apply(bpm.tom %>% 
                                         select(starts_with("bpm_")), 
                                       MARGIN = 2, FUN = function(x) {
                                         residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                       data = bpm.tom %>% mutate(y = x), 
                                                       family = poisson()))
                                       }))
# age, sex, interaction, and other ADHD meds
bpm.tom.asmeds.corrected <- cbind(bpm.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                                     apply(bpm.tom %>% 
                                             select(starts_with("bpm_")), 
                                           MARGIN = 2, FUN = function(x) {
                                             residuals(glm(y ~ interview_age + sex + interview_age:sex
                                                           + clonidine + adderall + concerta + vyvanse + ritalin +
                                                             intuniv + tenex + guanfacine + dexmethylphenidate
                                                           , 
                                                           data = bpm.tom %>% mutate(y = x), 
                                                           family = poisson()))
                                           }))
# combine raw and corrected sst data
bpm.tom.all <- inner_join(bpm.tom,
                             inner_join(bpm.tom.as.corrected %>% rename_at(.vars = vars(starts_with("bpm_")), 
                                                                              .funs = function(x) sub("bpm_raw_", "bpm_as_", x)),
                                        bpm.tom.asmeds.corrected %>% rename_at(.vars = vars(starts_with("bpm_")), 
                                                                                  .funs = function(x) sub("bpm_raw_", "bpm_asmeds_", x)))) %>% 
  group_by(IID) %>%
  filter(n_distinct(methylphenidate)==2) %>%
  group_by(IID, methylphenidate) %>%
  mutate_at(.vars = vars(starts_with("bpm_")), .funs = function(x) mean(x)) %>%
  ungroup() %>%
  distinct(IID, methylphenidate, .keep_all = T) %>%
  pivot_longer(cols = starts_with("bpm_"), names_to = "question", values_to = "val") %>%
  arrange(IID, question, methylphenidate) %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  pivot_wider(names_from = methylphenidate, values_from = val, id_cols = c(IID, sex, question)) %>%
  mutate(value_type = factor(ifelse(grepl("bpm_raw_", question), "raw data", ifelse(grepl("bpm_as_", question), 
                                                                                  "corrected for age, sex, and interaction", 
                                                                                  "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  mutate(delta = `1` - `0`) # being on MPH - not being on MPH
###
rm(bpm.tom.asmeds.corrected)
rm(bpm.tom.as.corrected)
gc()
###
################################################################################
################################## plots #######################################
# scatter plot for before and after score for raw, age&sex corrected, and age&sex&meds corrected
p3.s <- inner_join(bpm.tom.all, abcd.pgs) %>%
  mutate(ADHD_decile = as.factor(ntile(`ADHD-Demontis`, 10))) %>%
  mutate(question = as.factor(sub("bpm_.*?_","",question))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  ggplot(aes(x=`1`, y=`0`, alpha = ADHD_decile))+
  geom_point(position = "jitter", size=1) +
  scale_color_manual(values = redblu.col)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="on - MPH", y = "off - MPH")
# density plot for on/off MPH
p3.d <- bpm.tom.all %>%
  mutate(question = as.factor(sub("bpm_.*?_","",question))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  pivot_longer(cols = c(`0`,`1`), names_to = "MPH_status", values_to = "score") %>%
  ggplot(aes(x=score, color = MPH_status))+
  geom_density(position = "identity") +
  scale_color_manual(values = redblu.col)+
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="")
# paired boxplot
p3.b <- bpm.tom.all %>% 
  # filter(value_type == "raw data") %>%
  mutate(question = as.factor(sub("bpm_.*?_","",question))) %>%
  arrange(IID, question) %>%
  rename(off_MPH = `0`, on_MPH = `1`) %>%
  ggpaired(cond1 = "off_MPH", cond2 = "on_MPH", y = "val", 
           line.color = "gray", line.size = 0.1, palette = "jco", 
           ylab = "", xlab = "", point.size = 0.3) +
  facet_grid2(rows = vars(value_type), cols = vars(question), scales = "free_y", independent = "y")+
  theme_minimal() +
  stat_compare_means(paired = T, size = 3, method = "t.test")
# correlation between predicted response to MPH and participants performance in SST whether data is corrected or not
abcd.c.3 <- inner_join(abcd.pred, bpm.tom.all %>% select(-c(`0`,`1`))%>% pivot_wider(names_from = "question", values_from = "delta"))
p3.p <- corr.table(abcd.c.3%>%select(predicted), abcd.c.3 %>% select(starts_with("bpm_")), method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  mutate(question = as.factor(sub("bpm_.*?_","",V2))) %>%
  mutate(value_type = factor(ifelse(grepl("bpm_raw_", V2), "raw data", ifelse(grepl("bpm_as_", V2), 
                                                                            "corrected for age, sex, and interaction", 
                                                                            "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(pval < 0.1, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.3%>%distinct(IID))))
################################################################################
# ASR
abcd.asr <- read_csv(paste0(abcd.raw.dir, "/mental-health/mh_p_asr.csv"))
abcd.asr.filt <- abcd.asr %>%
  select(IID = src_subject_id, eventname, 
         dsm5_raw_perstr = asr_scr_perstr_r,
         syn_raw_anxdep = asr_scr_anxdep_r,
         syn_raw_withdrawn = asr_scr_withdrawn_r,
         syn_raw_somatic = asr_scr_somatic_r,
         syn_raw_thought = asr_scr_thought_r,
         syn_raw_attention = asr_scr_attention_r,
         syn_raw_aggressive = asr_scr_aggressive_r,
         syn_raw_rulebreak = asr_scr_rulebreak_r,
         syn_raw_intrusive = asr_scr_intrusive_r,
         syn_raw_internal = asr_scr_internal_r,
         syn_raw_external = asr_scr_external_r,
         syn_raw_totprob = asr_scr_totprob_r,
         dsm5_raw_depress = asr_scr_depress_r,
         dsm5_raw_anxdisord = asr_scr_anxdisord_r,
         dsm5_raw_somaticpr = asr_scr_somaticpr_r,
         dsm5_raw_avoidant = asr_scr_avoidant_r,
         dsm5_raw_adhd = asr_scr_adhd_r,
         dsm5_raw_antisocial = asr_scr_antisocial_r,
         dsm5_raw_inattention = asr_scr_inattention_r,
         dsm5_raw_hyperactive = asr_scr_hyperactive_r) %>%
  drop_na()
abcd.asr.filt <- inner_join(inner_join(abcd.asr.filt, 
                                        abcd.meds), demo) %>% select(IID, eventname, interview_age,sex, 
                                                                     starts_with("syn"), starts_with("dsm5"), 
                                                                     colnames(abcd.meds)) %>%
  drop_na()
#########################
# get participants with 2 or more measurements
asr.tom <- abcd.asr.filt %>%
  group_by(IID) %>%
  filter(n_distinct(eventname) >= 2, n_distinct(methylphenidate) >= 2) %>%
  ungroup()
# correct for age, sex, their interaction, other meds
# age-sex only
asr.tom.as.corrected <- cbind(asr.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                               apply(asr.tom %>% 
                                       select(starts_with("syn"), starts_with("dsm5")), 
                                     MARGIN = 2, FUN = function(x) {
                                       residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                     data = asr.tom %>% mutate(y = x), 
                                                     family = poisson()))
                                     }))
# age, sex, interaction, and other ADHD meds
asr.tom.asmeds.corrected <- cbind(asr.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                                   apply(asr.tom %>% 
                                           select(starts_with("syn"), starts_with("dsm5")), 
                                         MARGIN = 2, FUN = function(x) {
                                           residuals(glm(y ~ interview_age + sex + interview_age:sex
                                                         + clonidine + adderall + concerta + vyvanse + ritalin +
                                                           intuniv + tenex + guanfacine + dexmethylphenidate
                                                         , 
                                                         data = asr.tom %>% mutate(y = x), 
                                                         family = poisson()))
                                         }))
# combine raw and corrected sst data
asr.tom.all <- inner_join(asr.tom,
                           inner_join(asr.tom.as.corrected %>% rename_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), 
                                                                          .funs = function(x) ifelse(grepl("syn", x), sub("syn_raw_", "syn_as_", x),
                                                                                                     sub("dsm5_raw_", "dsm5_as_", x))),
                                      asr.tom.asmeds.corrected %>% rename_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), 
                                                                              .funs = function(x) ifelse(grepl("syn", x), sub("syn_raw_", "syn_asmeds_", x),
                                                                                                         sub("dsm5_raw_", "dsm5_asmeds_", x))))) %>% 
  group_by(IID) %>%
  filter(n_distinct(methylphenidate)==2) %>%
  group_by(IID, methylphenidate) %>%
  mutate_at(.vars = vars(starts_with("e_")), .funs = function(x) mean(x)) %>%
  ungroup() %>%
  distinct(IID, methylphenidate, .keep_all = T) %>%
  pivot_longer(cols = c(starts_with("syn"), starts_with("dsm5")), names_to = "question", values_to = "val") %>%
  arrange(IID, question, methylphenidate) %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  pivot_wider(names_from = methylphenidate, values_from = val, id_cols = c(IID, sex, question)) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", question), "raw data", ifelse(grepl("as_", question), 
                                                                                "corrected for age, sex, and interaction", 
                                                                                "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  mutate(subscale = ifelse(grepl("dsm5", question), "dsm5", "syn")) %>%
  mutate(delta = `1` - `0`) # being on MPH - not being on MPH
###
rm(asr.tom.asmeds.corrected)
rm(asr.tom.as.corrected)
gc()
###
################################################################################
################################## plots #######################################
# scatter plot for before and after score for raw, age&sex corrected, and age&sex&meds corrected
p4.s.syn <- inner_join(asr.tom.all, abcd.pgs) %>%
  mutate(ADHD_decile = as.factor(ntile(`ADHD-Demontis`, 10))) %>%
  filter(subscale=="syn") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  ggplot(aes(x=`1`, y=`0`, alpha = ADHD_decile))+
  geom_point(position = "jitter", size=1) +
  scale_color_manual(values = redblu.col)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="on - MPH", y = "off - MPH")
p4.s.dsm5 <- inner_join(asr.tom.all, abcd.pgs) %>%
  mutate(ADHD_decile = as.factor(ntile(`ADHD-Demontis`, 10))) %>%
  filter(subscale=="dsm5") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  ggplot(aes(x=`1`, y=`0`, alpha = ADHD_decile))+
  geom_point(position = "jitter", size=1) +
  scale_color_manual(values = redblu.col)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="on - MPH", y = "off - MPH")
# density plot for on/off MPH
p4.d.syn <- asr.tom.all %>%
  filter(subscale=="syn") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  pivot_longer(cols = c(`0`,`1`), names_to = "MPH_status", values_to = "score") %>%
  ggplot(aes(x=score, color = MPH_status))+
  geom_density(position = "identity") +
  scale_color_manual(values = redblu.col)+
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="")
p4.d.dsm5 <- asr.tom.all %>%
  filter(subscale=="dsm5") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  pivot_longer(cols = c(`0`,`1`), names_to = "MPH_status", values_to = "score") %>%
  ggplot(aes(x=score, color = MPH_status))+
  geom_density(position = "identity") +
  scale_color_manual(values = redblu.col)+
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="")
# paired boxplot
p4.b.syn <- asr.tom.all %>% 
  filter(subscale=="syn") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  arrange(IID, question) %>%
  rename(off_MPH = `0`, on_MPH = `1`) %>%
  ggpaired(cond1 = "off_MPH", cond2 = "on_MPH", y = "val", 
           line.color = "gray", line.size = 0.1, palette = "jco", 
           ylab = "", xlab = "", point.size = 0.3) +
  facet_grid2(rows = vars(value_type), cols = vars(question), scales = "free_y", independent = "y")+
  theme_minimal() +
  stat_compare_means(paired = T, size = 3, method = "t.test")
p4.b.dsm5 <- asr.tom.all %>% 
  filter(subscale=="dsm5") %>%
  mutate(question = as.factor(ifelse(grepl("syn", question), sub("syn_.*?_", "syn_", question),
                                     sub("dsm5_.*?_", "dsm5_", question)))) %>%
  arrange(IID, question) %>%
  rename(off_MPH = `0`, on_MPH = `1`) %>%
  ggpaired(cond1 = "off_MPH", cond2 = "on_MPH", y = "val", 
           line.color = "gray", line.size = 0.1, palette = "jco", 
           ylab = "", xlab = "", point.size = 0.3) +
  facet_grid2(rows = vars(value_type), cols = vars(question), scales = "free_y", independent = "y")+
  theme_minimal() +
  stat_compare_means(paired = T, size = 3, method = "t.test")
# correlation between predicted response to MPH and participants performance in SST whether data is corrected or not
abcd.c.4 <- inner_join(abcd.pred, asr.tom.all %>% select(-c(`0`,`1`))%>% pivot_wider(names_from = "question", values_from = "delta"))
p4.p <- corr.table(abcd.c.4%>%select(predicted), abcd.c.4 %>% select(starts_with("syn"), starts_with("dsm5")), method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  mutate(question = as.factor(ifelse(grepl("syn", V2), sub("syn_.*?_", "syn_", V2),
                                     sub("dsm5_.*?_", "dsm5_", V2)))) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                          "corrected for age, sex, and interaction", 
                                                                          "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(pval < 0.1, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.4%>%distinct(IID))))
################################################################################
# culture and environment variables 
abcd.ce <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/abcd5/culture-environment_c.rds") %>%
  mutate(psb_y_ss_mean = psb_y_ss_mean*3, psb_p_ss_mean = psb_p_ss_mean*3) %>%
  select(IID = src_subject_id, eventname, 
         ce_raw_wills_problem_solving = wps_ss_sum,
         # ce_raw_grades = sag_grades_last_yr,
         ce_raw_prosocial_beh_youth = psb_y_ss_mean,
         ce_raw_prosocial_beh_parent = psb_p_ss_mean)
abcd.ce <- inner_join(inner_join(abcd.ce, abcd.meds), demo) %>%
  drop_na()
#########################
#########################
# get participants with 2 or more measurements
ce.tom <- abcd.ce %>%
  group_by(IID) %>%
  filter(n_distinct(eventname) >= 2, n_distinct(methylphenidate) >= 2) %>%
  ungroup()
# correct for age, sex, their interaction, other meds
# age-sex only
ce.tom.as.corrected <- cbind(ce.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                              apply(ce.tom %>% 
                                      select(starts_with("ce_")), 
                                    MARGIN = 2, FUN = function(x) {
                                      residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                    data = ce.tom %>% mutate(y = x), 
                                                    family = poisson()))
                                    }))
# age, sex, interaction, and other ADHD meds
ce.tom.asmeds.corrected <- cbind(ce.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                                  apply(ce.tom %>% 
                                          select(starts_with("ce_")), 
                                        MARGIN = 2, FUN = function(x) {
                                          residuals(glm(y ~ interview_age + sex + interview_age:sex
                                                        + clonidine + adderall + concerta + vyvanse + ritalin +
                                                          intuniv + tenex + guanfacine + dexmethylphenidate
                                                        , 
                                                        data = ce.tom %>% mutate(y = x), 
                                                        family = poisson()))
                                        }))
# combine raw and corrected sst data
ce.tom.all <- inner_join(ce.tom,
                          inner_join(ce.tom.as.corrected %>% rename_at(.vars = vars(starts_with("ce_")), 
                                                                        .funs = function(x) sub("ce_raw_", "ce_as_", x)),
                                     ce.tom.asmeds.corrected %>% rename_at(.vars = vars(starts_with("ce_")), 
                                                                            .funs = function(x) sub("ce_raw_", "ce_asmeds_", x)))) %>% 
  group_by(IID) %>%
  filter(n_distinct(methylphenidate)==2) %>%
  group_by(IID, methylphenidate) %>%
  mutate_at(.vars = vars(starts_with("ce_")), .funs = function(x) mean(x)) %>%
  ungroup() %>%
  distinct(IID, methylphenidate, .keep_all = T) %>%
  pivot_longer(cols = starts_with("ce_"), names_to = "question", values_to = "val") %>%
  arrange(IID, question, methylphenidate) %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  pivot_wider(names_from = methylphenidate, values_from = val, id_cols = c(IID, sex, question)) %>%
  mutate(value_type = factor(ifelse(grepl("ce_raw_", question), "raw data", ifelse(grepl("ce_as_", question), 
                                                                                    "corrected for age, sex, and interaction", 
                                                                                    "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  mutate(delta = `1` - `0`) # being on MPH - not being on MPH
###
rm(ce.tom.asmeds.corrected)
rm(ce.tom.as.corrected)
gc()
###
################################################################################
################################## plots #######################################
# scatter plot for before and after score for raw, age&sex corrected, and age&sex&meds corrected
p5.s <- inner_join(ce.tom.all, abcd.pgs) %>%
  mutate(ADHD_decile = as.factor(ntile(`ADHD-Demontis`, 10))) %>%
  mutate(question = as.factor(sub("ce_.*?_","",question))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  ggplot(aes(x=`1`, y=`0`, alpha = ADHD_decile))+
  geom_point(position = "jitter", size=1) +
  scale_color_manual(values = redblu.col)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="on - MPH", y = "off - MPH")
# density plot for on/off MPH
p5.d <- ce.tom.all %>%
  mutate(question = as.factor(sub("ce_.*?_","",question))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  pivot_longer(cols = c(`0`,`1`), names_to = "MPH_status", values_to = "score") %>%
  ggplot(aes(x=score, color = MPH_status))+
  geom_density(position = "identity") +
  scale_color_manual(values = redblu.col)+
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="")

# paired boxplot
p5.b <- ce.tom.all %>% 
  # filter(value_type == "raw data") %>%
  mutate(question = as.factor(sub("ce_.*?_","",question))) %>%
  arrange(IID, question) %>%
  rename(off_MPH = `0`, on_MPH = `1`) %>%
  ggpaired(cond1 = "off_MPH", cond2 = "on_MPH", y = "val", 
           line.color = "gray", line.size = 0.1, palette = "jco", 
           ylab = "", xlab = "", point.size = 0.3) +
  facet_grid2(rows = vars(value_type), cols = vars(question), scales = "free_y", independent = "y")+
  theme_minimal() +
  stat_compare_means(paired = T, size = 3, method = "t.test")
# correlation between predicted response to MPH and participants performance in SST whether data is corrected or not
abcd.c.5 <- inner_join(abcd.pred, ce.tom.all %>% select(-c(`0`,`1`))%>% pivot_wider(names_from = "question", values_from = "delta"))
p5.p <- corr.table(abcd.c.5%>%select(predicted), abcd.c.5 %>% select(starts_with("ce_")), method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  mutate(question = as.factor(sub("ce_.*?_","",V2))) %>%
  mutate(value_type = factor(ifelse(grepl("ce_raw_", V2), "raw data", ifelse(grepl("ce_as_", V2), 
                                                                              "corrected for age, sex, and interaction", 
                                                                              "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(pval < 0.1, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.5%>%distinct(IID))))
################################################################################
pdf("figs/abcd5/sst-cbcl-asr-ce-mph-on-off-s-alpha-plots.pdf", width = 15)
print(p1.s+labs(title = "SST"));print(p2.s.dsm5+labs(title = "CBCL DSM5 by parent"))
print(p2.s.syn+labs(title = "CBCL syndromes by parent"));print(p3.s+labs(title = "BPM by youth"))
print(p4.s.dsm5+labs(title = "ASR DSM5 by parent"));print(p4.s.syn+labs(title = "ASR syndromes by parent"))
print(p5.s+labs(title = "Culture & Environment"))
dev.off()
pdf("figs/abcd5/sst-cbcl-asr-ce-mph-on-off-d-plots.pdf", width = 15)
print(p1.d+labs(title = "SST"));print(p2.d.dsm5+labs(title = "CBCL DSM5 by parent"))
print(p2.d.syn+labs(title = "CBCL syndromes by parent"));print(p3.d+labs(title = "BPM by youth"))
print(p4.d.dsm5+labs(title = "ASR DSM5 by parent"));print(p4.d.syn+labs(title = "ASR syndromes by parent"))
print(p5.d+labs(title = "Culture & Environment"))
dev.off()
pdf("figs/abcd5/sst-cbcl-asr-ce-mph-on-off-b-plots.pdf", width = 15)
print(p1.b+labs(title = "SST"));print(p2.b.dsm5+labs(title = "CBCL DSM5 by parent"))
print(p2.b.syn+labs(title = "CBCL syndromes by parent"));print(p3.b+labs(title = "BPM by youth"))
print(p4.b.dsm5+labs(title = "ASR DSM5 by parent"));print(p4.b.syn+labs(title = "ASR syndromes by parent"))
print(p5.b+labs(title = "Culture & Environment"))
dev.off()
patchwork::wrap_plots(p1.p+labs(title = "SST"), 
                      p2.p+labs(title = "CBCL by parent"), 
                      p3.p+labs(title = "BPM by youth"), 
                      p4.p+labs(title = "ASR by parent"), 
                      p5.p+labs(title = "Culture & Environment"), 
                      nrow = 1)
################################################################################
# pgs
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))

# correlation between pgs and participants performance in SST whether data is corrected or not
abcd.p.1 <- inner_join(abcd.pgs, sst.r1.tom.all %>% select(-c(`0`,`1`))%>% pivot_wider(names_from = "question", values_from = "delta"))
p1.p.p <- corr.table(abcd.p.1%>%select(colnames(abcd.pgs)[-1]), 
                     abcd.p.1 %>% select(starts_with("e_")), method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
  mutate(question = as.factor(sub("e_.*?_","",V2))) %>%
  mutate(value_type = factor(ifelse(grepl("e_raw_", V2), "raw data", ifelse(grepl("e_as_", V2), 
                                                                            "corrected for age, sex, and interaction", 
                                                                            "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.p.1%>%distinct(IID)), "\n",
                                              "* pval < 0.05", "\n", 
                                              "** pval < 0.01", "\n", 
                                              "*** FDR < 0.05", "\n", 
                                              det.cap))
# correlation between pgs and participants performance in SST whether data is corrected or not
abcd.p.2 <- inner_join(abcd.pgs, cbcl.tom.all %>% select(-c(`0`,`1`))%>% pivot_wider(names_from = "question", values_from = "delta"))
p2.p.p <- corr.table(abcd.p.2%>%select(colnames(abcd.pgs)[-1]), 
                     abcd.p.2 %>% select(starts_with("syn"), starts_with("dsm5")), method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
  mutate(question = as.factor(ifelse(grepl("syn", V2), sub("syn_.*?_", "syn_", V2),
                                     sub("dsm5_.*?_", "dsm5_", V2)))) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                          "corrected for age, sex, and interaction", 
                                                                          "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.p.2%>%distinct(IID))))
# correlation between pgs and participants performance in SST whether data is corrected or not
abcd.p.3 <- inner_join(abcd.pgs, bpm.tom.all %>% select(-c(`0`,`1`))%>% pivot_wider(names_from = "question", values_from = "delta"))
p3.p.p <- corr.table(abcd.p.3%>%select(colnames(abcd.pgs)[-1]), 
                     abcd.p.3 %>% select(starts_with("bpm_")), method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
  mutate(question = as.factor(sub("bpm_.*?_","",V2))) %>%
  mutate(value_type = factor(ifelse(grepl("bpm_raw_", V2), "raw data", ifelse(grepl("bpm_as_", V2), 
                                                                              "corrected for age, sex, and interaction", 
                                                                              "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.p.3%>%distinct(IID))))
# correlation between pgs and participants performance in SST whether data is corrected or not
abcd.p.4 <- inner_join(abcd.pgs, asr.tom.all %>% select(-c(`0`,`1`))%>% pivot_wider(names_from = "question", values_from = "delta"))
p4.p.p <- corr.table(abcd.p.4%>%select(colnames(abcd.pgs)[-1]), 
                   abcd.p.4 %>% select(starts_with("syn"), starts_with("dsm5")), method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
  mutate(question = as.factor(ifelse(grepl("syn", V2), sub("syn_.*?_", "syn_", V2),
                                     sub("dsm5_.*?_", "dsm5_", V2)))) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                          "corrected for age, sex, and interaction", 
                                                                          "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.p.4%>%distinct(IID))))
# correlation between pgs and participants performance in SST whether data is corrected or not
abcd.p.5 <- inner_join(abcd.pgs, ce.tom.all %>% select(-c(`0`,`1`))%>% pivot_wider(names_from = "question", values_from = "delta"))
p5.p.p <- corr.table(abcd.p.5%>%select(colnames(abcd.pgs)[-1]), 
                     abcd.p.5 %>% select(starts_with("ce_")), method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
  mutate(question = as.factor(sub("ce_.*?_","",V2))) %>%
  mutate(value_type = factor(ifelse(grepl("ce_raw_", V2), "raw data", ifelse(grepl("ce_as_", V2), 
                                                                             "corrected for age, sex, and interaction", 
                                                                             "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.p.5%>%distinct(IID))))
##############################
# combine pgs plots
patchwork::wrap_plots(p1.p.p+labs(title = "SST"), 
                      p2.p.p+labs(title = "CBCL by parent"), 
                      p3.p.p+labs(title = "BPM by youth"), 
                      p4.p.p+labs(title = "ASR by parent"), 
                      p5.p.p+labs(title = "Culture & Environment"), 
                      nrow = 1)
################################################################################
# same scatterplots as before, but change the colors to be ADHD pgs deciles
library(ggtern)
inner_join(sst.r1.tom.all, abcd.pgs) %>%
  mutate(ADHD_decile = as.factor(ntile(`ADHD-Demontis`, 10))) %>%
  mutate(question = as.factor(sub("ce_.*?_","",question))) %>%
  mutate(value_type = as.factor(value_type)) %>%
  mutate(ADHD_decile = as.factor(ntile(`ADHD-Demontis`, 10))) %>%
  ggtern(aes(x=`1`, y=`0`, z = `ADHD-Demontis`))+
  # ggplot(aes(x=`ADHD-Demontis`, y=ADHD_decile, alpha = ADHD_decile))+
  geom_point(size=1) +
  # scale_color_manual(values = redblu.col)+
  # geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid2(cols = vars(question), rows = vars(value_type), scales = "free", 
              independent = T) +
  labs(x="on - MPH", y = "off - MPH", z = "ADHD PGS") +
  theme_rgbw() +
  theme_hidetitles()
  theme_hidelabels()
  tlr2xy(coord = coord_tern(expand = F))
  coord_tern(expand = F)
  geom_Tline(Tintercept = 0)

pgs.pred <- inner_join(abcd.pgs,abcd.pred)
corr.table(pgs.pred %>% select(predicted),
           pgs.pred %>%select(colnames(abcd.pgs)[-1]), method = "spearman") %>%
  filter(V1 == "predicted", V2 != V1) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval < 0.05, "*", "")))+
  geom_tile()+
  geom_text(size=3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="")
################################################################################
################################################################################
################################################################################
# check correlation between cbcl and bpm
cbcl.bpm <- inner_join(abcd.cbcl.filt, abcd.bpm)
corr.table(cbcl.bpm %>% select(starts_with("bpm")), 
           cbcl.bpm %>% select(syn_raw_attention, syn_raw_external, syn_raw_internal, syn_raw_totprob),
           method = "spearman") %>%
  filter(V2 %in% c("syn_raw_attention", "syn_raw_external", "syn_raw_internal", "syn_raw_totprob"), 
         grepl("bpm", V1)) %>%
  ggplot(aes(x=V1, y=V2, fill=r, label = ifelse(pval < 1e-25, paste0(rho, ": ", round(r, 3), "\n",
                                                                     "p: *"), ""))) +
  geom_tile()+
  geom_text(size=3)+
  redblu.col.gradient + null_labs + my.guides +
  labs(caption = paste0("n(samples): ", nrow(cbcl.bpm), "\n", 
                        "* pval < 1e-25"))

################################################################################
################################################################################
# get correlation between cbcl deltas by item and pgs
abcd.cbcl <- read_csv(paste0(abcd.raw.dir, "/mental-health/mh_p_cbcl.csv"))
cbcl.scales <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/cbcl-scales.rds")
abcd.cbcl.filt <- abcd.cbcl %>%
  select(IID = src_subject_id, eventname, 
         starts_with("cbcl_q")) %>%
  drop_na()
abcd.cbcl.filt <- inner_join(inner_join(abcd.cbcl.filt, 
                                        abcd.meds), demo) %>% select(IID, eventname, interview_age,sex, 
                                                                     starts_with("cbcl_q"), 
                                                                     colnames(abcd.meds)) %>%
  drop_na()
colnames(abcd.cbcl.filt)[5:123] <- paste0("cbcl_", cbcl.scales$cbcl_item)
#########################
# get participants with 2 or more measurements
cbcl.tom <- abcd.cbcl.filt %>%
  group_by(IID) %>%
  filter(n_distinct(eventname) >= 2, n_distinct(methylphenidate) >= 2) %>%
  ungroup()
# correct for age, sex, their interaction, other meds
# age-sex only
cbcl.tom.as.corrected <- cbind(cbcl.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                               apply(cbcl.tom %>% 
                                       select(starts_with("cbcl_q")), 
                                     MARGIN = 2, FUN = function(x) {
                                       residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                     data = cbcl.tom %>% mutate(y = x), 
                                                     family = poisson()))
                                     }))
# age, sex, interaction, and other ADHD meds
cbcl.tom.asmeds.corrected <- cbind(cbcl.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                                   apply(cbcl.tom %>% 
                                           select(starts_with("cbcl_q")), 
                                         MARGIN = 2, FUN = function(x) {
                                           residuals(glm(y ~ interview_age + sex + interview_age:sex
                                                         + clonidine + adderall + concerta + vyvanse + ritalin +
                                                           intuniv + tenex + guanfacine + dexmethylphenidate
                                                         , 
                                                         data = cbcl.tom %>% mutate(y = x), 
                                                         family = poisson()))
                                         }))
# combine raw and corrected sst data
cbcl.tom.all <- inner_join(cbcl.tom,
                           inner_join(cbcl.tom.as.corrected %>% rename_at(.vars = vars(starts_with("cbcl_q")), 
                                                                          .funs = function(x) sub("cbcl_", "cbcl_as_", x)),
                                      cbcl.tom.asmeds.corrected %>% rename_at(.vars = vars(starts_with("cbcl_q")), 
                                                                              .funs = function(x) sub("cbcl_", "cbcl_asmeds_", x)))) %>% 
  group_by(IID) %>%
  filter(n_distinct(methylphenidate)==2) %>%
  group_by(IID, methylphenidate) %>%
  mutate_at(.vars = vars(starts_with("cbcl_")), .funs = function(x) mean(x)) %>%
  ungroup() %>%
  distinct(IID, methylphenidate, .keep_all = T) %>%
  pivot_longer(cols = c(starts_with("cbcl_")), names_to = "question", values_to = "val") %>%
  arrange(IID, question, methylphenidate) %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  pivot_wider(names_from = methylphenidate, values_from = val, id_cols = c(IID, sex, question)) %>%
  mutate(value_type = factor(ifelse(grepl("cbcl_q", question), "raw data", ifelse(grepl("as_", question), 
                                                                                "corrected for age, sex, and interaction", 
                                                                                "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  mutate(delta = `1` - `0`) # being on MPH - not being on MPH
###
rm(cbcl.tom.asmeds.corrected)
rm(cbcl.tom.as.corrected)
gc()
###
################################################################################
################################## plots
# correlation between predicted response to MPH and participants cbcl whether data is corrected or not
abcd.c.2 <- inner_join(abcd.pred, cbcl.tom.all %>% select(-c(`0`,`1`))%>% 
                         pivot_wider(names_from = "question", values_from = "delta"))
corr.table(abcd.c.2%>%select(predicted), abcd.c.2 %>% select(starts_with("cbcl_")), method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  mutate(question = V2) %>%
  mutate(value_type = factor(ifelse(grepl("cbcl_q", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                          "corrected for age, sex, and interaction", 
                                                                          "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  filter(grepl("cbcl_q", V2)) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(pval < 0.05, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.2%>%distinct(IID))))

################################################################################
# correlation between cbcl items delta and pgs
abcd.c.10 <- inner_join(abcd.pgs, cbcl.tom.all %>% select(-c(`0`,`1`))%>% 
                         pivot_wider(names_from = "question", values_from = "delta"))
corr.table(abcd.c.10%>%select(colnames(abcd.pgs)[-1]), abcd.c.10 %>% select(starts_with("cbcl_")), method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs)[-1], ! V2 %in% colnames(abcd.pgs)[-1]) %>%
  mutate(question = V2) %>%
  mutate(value_type = factor(ifelse(grepl("cbcl_q", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                            "corrected for age, sex, and interaction", 
                                                                            "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  filter(grepl("cbcl_q", V2)) %>%
  ggplot(aes(x=V1, y=question, fill = r, label = ifelse(pval < 0.05, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.10%>%distinct(IID))))+
  theme(axis.text.y.left = element_text(size=7), 
        axis.text.x.bottom = element_text(angle = 0, hjust = 0.5))
################################################################################
# correlation between taking the drug and cbcl
abcd.cbcl <- read_csv(paste0(abcd.raw.dir, "/mental-health/mh_p_cbcl.csv"))
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
abcd.cbcl.filt <- inner_join(inner_join(abcd.cbcl.filt, 
                                        abcd.meds), demo) %>% select(IID, eventname, interview_age,sex, 
                                                                     starts_with("syn"), starts_with("dsm5"), 
                                                                     colnames(abcd.meds)) %>%
  drop_na()
#########################
# correct for age, sex, their interaction
# age-sex only
cbcl.as.corrected <- cbind(abcd.cbcl.filt %>% select(IID, interview_age, sex, eventname), 
                           apply(abcd.cbcl.filt %>% 
                                   select(starts_with("syn"), starts_with("dsm5")), 
                                 MARGIN = 2, FUN = function(x) {
                                   residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                 data = abcd.cbcl.filt %>% mutate(y = x), 
                                                 family = poisson()))
                                 }))
# combine raw and corrected cbcl data
cbcl.all <- inner_join(abcd.cbcl.filt,
                       cbcl.as.corrected %>% rename_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), 
                                                       .funs = function(x) ifelse(grepl("syn", x), sub("syn_raw_", "syn_as_", x),
                                                                                  sub("dsm5_raw_", "dsm5_as_", x))))
############ plot
corr.table(cbcl.all %>% select(colnames(abcd.meds)[-c(1:2)]),
           cbcl.all %>% select(starts_with("syn"), starts_with("dsm5")),
           method = "spearman") %>%
  filter(V1 %in% colnames(abcd.meds)[-c(1:2)], ! V2 %in% colnames(abcd.meds)[-c(1:2)]) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                          "corrected for age, sex, and interaction", 
                                                                          "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", nrow(cbcl.all), "\n",
                        "* pval < 0.05", "\n", 
                        "** pval < 0.01", "\n", 
                        "*** FDR < 0.05"), 
       title = "correlation of taking a medication and cbcl")
################################################################################
# get cbcl deltas by drug status for participants
adhd.meds.deltas <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d <- adhd.meds$drug[i]
  # print(i)
  n <- nrow(cbcl.all %>%
                rename(drug = which(colnames(cbcl.all)==d)) %>%
                group_by(IID) %>%
                filter(n_distinct(eventname) >= 2, n_distinct(drug) >= 2) %>%
                ungroup() %>%
                distinct(IID))
  t <- cbcl.all %>%
    rename(drug = which(colnames(cbcl.all)==d)) %>%
    group_by(IID) %>%
    filter(n_distinct(eventname) >= 2, n_distinct(drug) >= 2) %>%
    ungroup() %>%
    group_by(IID, drug) %>%
    mutate_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), .funs = function(x) mean(x)) %>%
    ungroup() %>%
    distinct(IID, drug, .keep_all = T) %>%
    pivot_longer(cols = c(starts_with("syn"), starts_with("dsm5")), names_to = "question", values_to = "val") %>%
    arrange(IID, question, drug) %>%
    mutate(drug = as.factor(drug)) %>%
    pivot_wider(names_from = drug, values_from = val, id_cols = c(IID, sex, question)) %>%
    mutate(delta = `1` - `0`) %>% 
    mutate(drug = d) %>%
    mutate(n_samples = n)
  return(t)
}
# plot average delta per drug per cbcl
adhd.meds.deltas %>%
  group_by(question, drug) %>%
  summarise_at(.vars = vars(delta, n_samples), .funs = function(x) mean(x)) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", question), "raw data", ifelse(grepl("as_", question), 
                                                                          "corrected for age, sex, and interaction", 
                                                                          "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  ggplot(aes(x=drug, y=question, fill = delta)) +
  geom_tile()+
  geom_point(aes(size = n_samples), alpha = 0.1) +
  my.guides+null_labs+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "average delta (on_drug - off_drug)")+
  labs(title = "average deltas for cbcl problems per drug")
adhd.meds.deltas %>%
  group_by(question, drug) %>%
  summarise_at(.vars = vars(delta, n_samples), .funs = function(x) median(x)) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", question), "raw data", ifelse(grepl("as_", question), 
                                                                                "corrected for age, sex, and interaction", 
                                                                                "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  ggplot(aes(x=drug, y=question, fill = delta)) +
  geom_tile()+
  geom_point(aes(size = n_samples), alpha = 0.1) +
  my.guides+null_labs+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "average delta (on_drug - off_drug)")+
  labs(title = "median of deltas for cbcl problems per drug")
################################################################################
# correlation between predicted response and pgs
pgs.pred <- inner_join(abcd.pred, abcd.pgs)
corr.table(pgs.pred %>% select(predicted),
           pgs.pred %>% select(colnames(abcd.pgs)[-1]),
           method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, paste0(rho, ": ", round(r, 3), "\n",
                                                                    "p: ", round(pval, 5)), "")))+
  geom_tile() + geom_text(size=3)+
  my.guides+redblu.col.gradient+null_labs+
  labs(caption = paste0("n(samples): ", nrow(pgs.pred)),
       title = "correlation between predicted MPH response and PGS")+
  theme(plot.title = element_text(size=9, hjust = 1.3))
################################################################################
# average of SST deltas by med
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
                    left_join(abcd.meds, demo)) %>% drop_na() %>%
  mutate_at(.vars = vars(3:7, interview_age), .funs = function(x) as.numeric(x))
# correct for age, sex, their interaction
sst.r1.as.corrected <- cbind(sst.r1 %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                                 apply(sst.r1 %>% 
                                         select(starts_with("e_")), 
                                       MARGIN = 2, FUN = function(x) {
                                         residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                       data = sst.r1 %>% mutate(y = x), 
                                                       family = poisson()))
                                       }))
# combine raw and corrected sst data
sst.all <- inner_join(sst.r1,
                       sst.r1.as.corrected %>% rename_at(.vars = vars(starts_with("e_")), 
                                                       .funs = function(x) sub("e_raw_", "e_as_", x)))
############ plot
corr.table(sst.all %>% select(colnames(abcd.meds)[-c(1:2)]),
           sst.all %>% select(starts_with("e_")),
           method = "spearman") %>%
  filter(V1 %in% colnames(abcd.meds)[-c(1:2)], ! V2 %in% colnames(abcd.meds)[-c(1:2)]) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                          "corrected for age, sex, and interaction", 
                                                                          "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", nrow(sst.all), "\n",
                        "* pval < 0.05", "\n", 
                        "** pval < 0.01", "\n", 
                        "*** FDR < 0.05"), 
       title = "correlation of taking a medication and cbcl")
################################################################################
# get sst deltas by drug status for participants
sst.meds.deltas <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d <- adhd.meds$drug[i]
  # print(i)
  n <- nrow(sst.all %>%
              rename(drug = which(colnames(sst.all)==d)) %>%
              group_by(IID) %>%
              filter(n_distinct(eventname) >= 2, n_distinct(drug) >= 2) %>%
              ungroup() %>%
              distinct(IID))
  t <- sst.all %>%
    rename(drug = which(colnames(sst.all)==d)) %>%
    group_by(IID) %>%
    filter(n_distinct(eventname) >= 2, n_distinct(drug) >= 2) %>%
    ungroup() %>%
    group_by(IID, drug) %>%
    mutate_at(.vars = vars(starts_with("e_")), .funs = function(x) mean(x)) %>%
    ungroup() %>%
    distinct(IID, drug, .keep_all = T) %>%
    pivot_longer(cols = c(starts_with("e_")), names_to = "question", values_to = "val") %>%
    arrange(IID, question, drug) %>%
    mutate(drug = as.factor(drug)) %>%
    pivot_wider(names_from = drug, values_from = val, id_cols = c(IID, sex, question)) %>%
    mutate(delta = `1` - `0`) %>% 
    mutate(drug = d) %>%
    mutate(n_samples = n)
  return(t)
}
sst.meds.deltas %>%
  group_by(question, drug) %>%
  summarise_at(.vars = vars(delta, n_samples), .funs = function(x) mean(x)) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", question), "raw data", ifelse(grepl("as_", question), 
                                                                                "corrected for age, sex, and interaction", 
                                                                                "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  ggplot(aes(x=drug, y=question, fill = delta)) +
  geom_tile()+
  geom_point(aes(size = n_samples), alpha = 0.1) +
  my.guides+null_labs+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "average delta (on_drug - off_drug)")+
  labs(title = "average deltas for SST per drug")
################################################################################
# correlation between sst deltas per drug and pgs
sst.deltas.pgs <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d = adhd.meds$drug[i]
  t <- inner_join(sst.meds.deltas %>% filter(drug == d) %>% pivot_wider(names_from = question, 
                                                                   values_from = delta, 
                                                                   id_cols = c(IID, sex, drug)),
             abcd.pgs)
  ret <- corr.table(t %>% select(colnames(abcd.pgs)[-1]),
             t %>% select(starts_with("e_")),
             method = "spearman") %>%
    filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
    mutate(value_type = factor(ifelse(grepl("raw_", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                            "corrected for age, sex, and interaction", 
                                                                            "corrected for age, sex, interaction, and other ADHD meds")), 
                               levels = c("raw data", "corrected for age, sex, and interaction", 
                                          "corrected for age, sex, interaction, and other ADHD meds"))) %>%
    group_by(value_type) %>%
    mutate(FDR = p.adjust(pval, method = "fdr")) %>%
    mutate(drug = d)
  return(ret)
}

sst.deltas.pgs %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), cols = vars(drug), scales = "free_y") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", nrow(t), "\n",
                        "* pval < 0.05", "\n", 
                        "** pval < 0.01", "\n", 
                        "*** FDR < 0.05"), 
       title = "correlation of taking a medication and cbcl")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################### dropped ####################################
################################################################################
################################################################################
# little man task
abcd.lmt <- read_csv(paste0(abcd.raw.dir, "/neurocognition/nc_y_lmt.csv")) %>%
  select(IID = src_subject_id, eventname, 
         corr = lmt_scr_num_correct, wrong = lmt_scr_num_wrong) %>%
  drop_na()
abcd.lmt <- inner_join(inner_join(abcd.lmt, abcd.meds), demo) %>%
  drop_na()
lmt.corrected <- cbind(abcd.lmt %>% select(-c(corr, wrong)), 
                       lapply(abcd.lmt %>% select(corr, wrong), FUN = function(x) 
                         residuals(glm(y ~ interview_age + sex + interview_age:sex 
                                       + clonidine + adderall + concerta + vyvanse + dextroamphetamine +
                                         ritalin + intuniv + tenex + amphetamine+ guanfacine + lisdexamfetamine +
                                         dexmethylphenidate
                                       , data = abcd.lmt %>% mutate(y=x), family = poisson())))) %>%
  group_by(IID) %>%
  filter(n_distinct(methylphenidate)==2) %>%
  group_by(IID, methylphenidate) %>%
  mutate_at(.vars = vars(corr, wrong), .funs = function(x) mean(x)) %>%
  ungroup() %>%
  distinct(IID, methylphenidate, .keep_all = T)
####################################
ps7 <- lmt.corrected %>% 
  pivot_longer(cols = c("corr", "wrong"), names_to = "question", values_to = "score") %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  ggpaired(x = "methylphenidate", y = "score",  group = "methylphenidate",
           color = "methylphenidate", 
           line.color = "gray", line.size = 0.1, palette = "jco") +
  facet_wrap("question", scales = "free_y", ncol = 2) + theme_minimal() +
  stat_compare_means(paired = T)
#########################
# get deltas per participant
# being on MPH - not being on MPH
lmt.deltas <- do.call(cbind, apply(lmt.corrected %>% select(corr, wrong),
                                   MARGIN = 2, FUN = function(x) lmt.corrected %>% mutate(y = x) %>%
                                     group_by(IID) %>%                                        
                                     summarise(delta = (y[methylphenidate== 1] - y[methylphenidate== 0])))) %>%
  select(IID = 1, ends_with("delta"))
abcd.c.7 <- inner_join(lmt.deltas, abcd.pred)

p7 <-  corr.table(abcd.c.7%>%select(predicted), abcd.c.7 %>% select(ends_with("delta")), method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval < 0.1, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.7)))
################################################################################
################################################################################

