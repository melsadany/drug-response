################################################################################
#                          CBCL/SST change by MPH status                       #
################################################################################
rm(list = ls())
gc()
source("~/LSS/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
project.dir <- "~/LSS/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response"
setwd(project.dir)
################################################################################
################################################################################
abcd.raw.dir <- "~/LSS/jmichaelson-sdata/ABCD/abcd_release_5_0/core"

# demographics
age <- read_csv(paste0(abcd.raw.dir, "/abcd-general/abcd_y_lt.csv")) %>% 
  select(IID = src_subject_id, eventname, interview_age)
sex <- read_csv(paste0(abcd.raw.dir, "/gender-identity-sexual-health/gish_p_gi.csv")) %>%
  mutate(sex = ifelse(demo_sex_v2 == 2, "Female", ifelse(demo_sex_v2 == 1, "Male", ifelse(demo_sex_v2 == 3, "intersex_M", NA)))) %>%
  select(IID = src_subject_id, sex) %>%
  distinct(IID, .keep_all = T)
demo <- full_join(age, sex)
rm(age);rm(sex);gc()

# predicted MPH response
abcd.pred <- read_rds("../data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(predicted = m) %>% mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
# ABCD PGS
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))
# ABCD meds
adhd.meds <- data.frame(drug = c("methylphenidate", "adderall", "concerta", "vyvanse", 
                                 "dextroamphetamine",  "ritalin", "intuniv", "strattera",
                                 "tenex", "amphetamine", "dexmethylphenidate", "lisdexamfetamine",
                                 "atomoxetine", "clonidine", "guanfacine"))
abcd.meds <- read_rds("~/LSS/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/abcd5-meds-matrix.rds") %>%
  as.data.frame() %>%
  select(c(1:2), any_of(adhd.meds$drug))
################################################################################
################################################################################
# ABCD SST
# decided to keep count of answers in trial run 1 (r1) only
sst.raw <- read_csv(paste0(abcd.raw.dir, "/imaging/mri_y_tfmr_sst_beh.csv"))
sst.r1 <- left_join(sst.raw %>%
                      select(IID = src_subject_id, eventname,
                             # e_raw_correct_go = tfmri_sst_r1_beh_crgo_nt,
                             # e_raw_correct_stop = tfmri_sst_r1_beh_crs_nt,
                             e_raw_stop_doesnot_stop = tfmri_sst_r1_beh_ssds_nt,
                             e_raw_no_response_on_go = tfmri_sst_r1_beh_nrgo_nt,
                             tfmri_sst_beh_switchflag) %>%
                      filter(is.na(tfmri_sst_beh_switchflag) == F),
                    left_join(abcd.meds, demo)) %>% drop_na() %>%
  mutate_at(.vars = vars(3:5, interview_age), .funs = function(x) as.numeric(x))
rm(sst.raw);gc()
#################
# clean SST data and correct for age, sex, interaction, and taking other meds than MPH
# then get the delta of being on MPH - off MPH

# get participants with 2 or more measurements
sst.r1.tom <- sst.r1 %>%
  group_by(IID) %>%
  filter(n_distinct(eventname) >= 2& n_distinct(methylphenidate) >= 2) %>%
  ungroup()
# correct for age, sex, their interaction, other meds
# then, get the delta
sst.r1.tom.asmeds.corrected <- cbind(sst.r1.tom %>% 
                                       select(IID, interview_age, sex, eventname, methylphenidate), 
                                     apply(sst.r1.tom %>% 
                                             select(starts_with("e_")), 
                                           MARGIN = 2, FUN = function(x) {
                                             residuals(glm(y ~ interview_age + sex + interview_age:sex
                                                           + clonidine + adderall + concerta + vyvanse + ritalin +
                                                             intuniv + tenex + guanfacine + dexmethylphenidate, 
                                                           data = sst.r1.tom %>% mutate(y = x), 
                                                           family = poisson()))
                                           })) %>%
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
  mutate(delta = `1` - `0`)


################################################################################
################################################################################
# ABCD CBCL
abcd.cbcl <- read_csv(paste0(abcd.raw.dir, "/mental-health/mh_p_cbcl.csv")) %>%
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

abcd.cbcl.filt <- inner_join(inner_join(abcd.cbcl,abcd.meds),demo) %>% 
  select(IID, eventname, interview_age,sex, 
         starts_with("syn"), starts_with("dsm5"), 
         colnames(abcd.meds)) %>%
  drop_na()
#########################
# get participants with 2 or more measurements
cbcl.tom <- abcd.cbcl.filt %>%
  group_by(IID) %>%
  filter(n_distinct(eventname) >= 2& n_distinct(methylphenidate) >1) %>%
  ungroup()

# correct for age, sex, their interaction, other meds than MPH
# then get the delta of on MPH - off MPH
cbcl.tom.asmeds.corrected <- cbind(cbcl.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                                   apply(cbcl.tom %>% 
                                           select(starts_with("syn"), starts_with("dsm5")), 
                                         MARGIN = 2, FUN = function(x) {
                                           residuals(glm(y ~ interview_age + sex + interview_age:sex
                                                         + clonidine + adderall + concerta + vyvanse + ritalin +
                                                           intuniv + tenex + guanfacine + dexmethylphenidate, 
                                                         data = cbcl.tom %>% mutate(y = x), 
                                                         family = poisson()))
                                         })) %>%
  group_by(IID, methylphenidate) %>%
  mutate_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), .funs = function(x) mean(x)) %>%
  ungroup() %>%
  distinct(IID, methylphenidate, .keep_all = T) %>%
  pivot_longer(cols = c(starts_with("syn"), starts_with("dsm5")), names_to = "question", values_to = "val") %>%
  arrange(IID, question, methylphenidate) %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  pivot_wider(names_from = methylphenidate, values_from = val, id_cols = c(IID, sex, question)) %>%
  mutate(subscale = ifelse(grepl("dsm5", question), "dsm5", "syn")) %>%
  mutate(delta = `1` - `0`) # being on MPH - not being on MPH

################################################################################
################################################################################
# combine deltas
abcd.deltas <- full_join(sst.r1.tom.asmeds.corrected %>%
                           pivot_wider(names_from = "question", values_from = "delta", id_cols = "IID"),
                         cbcl.tom.asmeds.corrected %>%
                           pivot_wider(names_from = "question", values_from = "delta", id_cols = "IID")) %>%
  rename_all(.funs = function(x) sub("_raw", "", sub("e_", "SST_", x)))
################################################################################
# correlation between PGS and 
deltas.pgs <- inner_join(abcd.deltas, abcd.pgs) %>% left_join(abcd.pred)
corr.table(deltas.pgs %>% select(starts_with("SST"), starts_with("dsm5"), starts_with("syn")),
           deltas.pgs %>% select(colnames(abcd.pgs)[-1], "predicted")) %>%
  mutate(FDR = p.adjust(pval, method = "fdr"),
         cat = ifelse(grepl("SST", V1), "SST", ifelse(grepl("syn", V1), "CBCL subscales", "DSM5 subscales"))) %>%
  filter(V1 %in% c(sub("_raw", "", cbcl.tom.asmeds.corrected$question),
                   sub("_raw", "", sub("e_", "SST_", sst.r1.tom.asmeds.corrected$question))),
         V2 %in% c(colnames(abcd.pgs), "predicted")) %>% 
  mutate(V1 = factor(V1, levels= c(colnames(deltas.pgs %>% select(starts_with("SST"), 
                                                                  starts_with("dsm5"), 
                                                                  starts_with("syn"))))),
         V1 = sub("SST_", "", V1), V1 = sub("syn_", "", sub("dsm5_", "", V1))) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "**", ifelse(pval < 0.05, "*",ifelse(pval < 0.1, ".", ""))))) +
  geom_tile() + geom_text(color = "white") +
  ggh4x::facet_grid2(cols = vars(cat), scales = "free", space = "free") +
  redblack.col.gradient + my.guides
ggsave("figs/0724_report/ABCD-MPH-CBCL-and-SST-deltas-w-PGS-and-predicted.png", bg = "white",
       width = 8, height = 8, units = "in", dpi = 360)
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
################################################################################
################################################################################
################################################################################