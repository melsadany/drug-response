################################################################################
#                   predcit drug response for ABCD participants                #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response"
setwd(project.dir)
################################################################################
################################################################################
################################################################################
abcd.raw.dir <- "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core"
sst.raw <- read_csv(paste0(abcd.raw.dir, "/imaging/mri_y_tfmr_sst_beh.csv"))
question.dict <- data.frame(q0 = colnames(sst.raw), description = t(sst.raw)[,1], row.names = NULL)

age <- read_csv(paste0(abcd.raw.dir, "/abcd-general/abcd_y_lt.csv")) %>% 
  select(IID = src_subject_id, eventname, interview_age)
sex <- read_csv(paste0(abcd.raw.dir, "/gender-identity-sexual-health/gish_p_gi.csv")) %>%
  mutate(sex = ifelse(demo_sex_v2 == 2, "Female", ifelse(demo_sex_v2 == 1, "Male", ifelse(demo_sex_v2 == 3, "intersex_M", NA)))) %>%
  select(IID = src_subject_id, sex) %>%
  distinct(IID, .keep_all = T)
demo <- full_join(age, sex)
#########################
# meds
adhd.meds <- data.frame(group = c(rep("stim", 5), rep("nonstim", 4), rep("antidep", 4)), 
                        drug = c("methylphenidate", "dextroamphetamine", 
                                 "amphetamine", "dexmethylphenidate", "lisdexamfetamine",
                                 "atomexetine", "clonidine", "guanfacine", "viloxazine",
                                 "bupropion", "desipramine", "imipramine", "nortriptyline"))
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/abcd5-meds-matrix.rds") %>%
  as.data.frame() %>%
  select(c(1:2), any_of(adhd.meds$drug[adhd.meds$group %in% c("stim", "nonstim")]))


#########################
# decided to keep count of answers in trial run 1 (r1) only
sst.r1 <- left_join(sst.raw[-1,] %>%
                      select(IID = src_subject_id, eventname,
                             ends_with("_nt")&contains("r1"), 
                             tfmri_sst_beh_switchflag,tfmri_sst_beh_performflag,
                             -c(contains("_beh_go_nt"), contains("_beh_s_nt"))) %>%
                      filter(is.na(tfmri_sst_beh_switchflag) == F),
                    left_join(abcd.meds, demo)) %>% drop_na() %>%
  mutate_at(.vars = vars(ends_with("_nt"), ends_with("flag"), ends_with("age")), .funs = function(x) as.numeric(x))
hist(sst.r1)
#########################
# get participants with 2 or more measurements
sst.r1.tom <- sst.r1 %>%
  group_by(IID) %>%
  filter(n_distinct(eventname) >= 2, sum(methylphenidate) >= 1) %>%
  ungroup()

# correct for age, sex, their interaction, other meds
sst.r1.tom.corrected <- cbind(sst.r1.tom %>% select(IID, interview_age, sex, eventname, methylphenidate), 
                              apply(sst.r1.tom %>% 
                                select(starts_with("tfmri"), -contains("flag")), 
                              MARGIN = 2, FUN = function(x) {
                                residuals(glm(y ~ interview_age + sex + interview_age:sex
                                              + clonidine + guanfacine + dexmethylphenidate
                                              , 
                                                                           data = sst.r1.tom %>% mutate(y = x), 
                                                                           family = poisson()))
                                })) %>%
  group_by(IID) %>%
  filter(n_distinct(methylphenidate)==2) %>%
  group_by(IID, methylphenidate) %>%
  mutate_at(.vars = vars(starts_with("tfmri")), .funs = function(x) mean(x)) %>%
  ungroup() %>%
  distinct(IID, methylphenidate, .keep_all = T)
#########################
library(ggpubr)
ps1 <- sst.r1.tom.corrected %>% 
  select(-contains("flag")) %>%
  pivot_longer(cols = starts_with("tfmri"), names_to = "question", values_to = "val") %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  ggpaired(x = "methylphenidate", y = "val",  group = "methylphenidate",
           color = "methylphenidate", 
           line.color = "gray", line.size = 0.4, palette = "jco") +
  facet_wrap("question", scales = "free_y", ncol = 4) + theme_minimal() +
  stat_compare_means(paired = T)
#########################
# get deltas per participant
# being on MPH - not being on MPH
sst.r1.deltas <- do.call(cbind, apply(sst.r1.tom.corrected %>% select(starts_with("tfmri")), 
                                      MARGIN = 2, FUN = function(x) sst.r1.tom.corrected %>% mutate(y = x) %>%
                                        group_by(IID) %>%
                                        summarise(delta = (y[methylphenidate== 1] - y[methylphenidate== 0])))) %>%
  select(IID = 1, ends_with("delta"))

################################################################################
abcd.pred <- read_rds("../data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(predicted = m) %>%
  # filter(IID %in% sst.r1.deltas$IID) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
abcd.c <- inner_join(abcd.pred, sst.r1.deltas)

p1 <- corr.table(abcd.c%>%select(predicted), abcd.c %>% select(starts_with("tfmri")), method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval < 0.1, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c), "\n",
                                              "delta per question = score_on_MPH - score_off_MPH", "\n", 
                                              "scores were first derived from:", "\n",
                                              "residuals(glm(y ~ interview_age + sex + interview_age:sex", "\n", 
                                              "\t\t\t + clonidine + guanfacine + dexmethylphenidate, family = poisson()))"))
################################################################################
# cbcl
abcd.cbcl <- read_csv(paste0(abcd.raw.dir, "/mental-health/mh_p_cbcl.csv"))
abcd.cbcl.imp <- abcd.cbcl %>%
  select(src_subject_id, eventname, 
         syn_attention = cbcl_scr_syn_attention_r,
         syn_anxdep = cbcl_scr_syn_anxdep_r,
         syn_withdep = cbcl_scr_syn_withdep_r,
         syn_somatic = cbcl_scr_syn_somatic_r,
         syn_social = cbcl_scr_syn_social_r,
         syn_thought = cbcl_scr_syn_thought_r,
         syn_rulebreak = cbcl_scr_syn_rulebreak_r,
         syn_aggressive = cbcl_scr_syn_aggressive_r,
         syn_internal = cbcl_scr_syn_internal_r,
         syn_external = cbcl_scr_syn_external_r,
         syn_totprob = cbcl_scr_syn_totprob_r) %>%
  drop_na()

abcd.cbcl.imp <- inner_join(inner_join(abcd.cbcl.imp %>% rename(IID = src_subject_id), 
                            abcd.meds), demo)

abcd.cbcl.imp <- abcd.cbcl.imp %>% select(IID, eventname, interview_age,sex, 
                                          starts_with("syn"), colnames(abcd.meds)) %>%
  drop_na()
cbcl.corrected <- cbind(abcd.cbcl.imp %>% select(IID , sex, eventname, methylphenidate), 
                        apply(abcd.cbcl.imp%>%select(starts_with("syn")), MARGIN = 2, 
                              function(x) residuals(glm(y ~ interview_age + sex + interview_age:sex + 
                                                          clonidine + guanfacine, 
                                                        data = abcd.cbcl.imp %>% mutate(y = x), 
                                                        family = poisson()))))
ps2 <- cbcl.corrected %>% 
  group_by(IID) %>%
  filter(n_distinct(eventname) >= 2, sum(methylphenidate) >= 1) %>% 
  ungroup() %>% 
  pivot_longer(cols = starts_with("syn"), names_to = "question", values_to = "val") %>%
  mutate(methylphenidate == as.factor(methylphenidate)) %>%
  group_by(IID, methylphenidate, question) %>%
  summarise(avg = mean(val)) %>%
  ungroup() %>%
  group_by(IID) %>%
  filter(n() == 22) %>%
  ggpaired(x = "methylphenidate", y = "avg",  group = "methylphenidate",
           color = "methylphenidate", 
           line.color = "gray", line.size = 0.4, palette = "jco") +
  facet_wrap("question", scales = "free_y", ncol = 4) + theme_minimal() +
  stat_compare_means(paired = T)

cbcl.deltas <- do.call(cbind, lapply(cbcl.corrected %>% 
                                       group_by(IID) %>%
                                       filter(n_distinct(eventname) >= 2, sum(methylphenidate) >= 1) %>% 
                                       ungroup() %>%
                                       select(starts_with("syn")), 
                                     FUN = function(x) 
                                       cbcl.corrected %>% 
                                       group_by(IID) %>%
                                       filter(n_distinct(eventname) >= 2, sum(methylphenidate) >= 1) %>%
                                       ungroup() %>%
                                       mutate(y = x) %>%
                                       mutate(methylphenidate = as.factor(methylphenidate)) %>%
                                       group_by(IID, methylphenidate) %>%
                                       summarise(avg = mean(y)) %>%
                                       pivot_wider(names_from = methylphenidate, values_from = avg) %>%
                                       rename(off = `0`, on = `1`) %>%
                                       drop_na() %>%                                       
                                       mutate(delta = (on - off)) %>%
                                       select(IID, delta))) %>%
  select(IID = 1, starts_with("delta"))
colnames(cbcl.deltas)[2:ncol(cbcl.deltas)] <- paste0(colnames(cbcl.corrected%>%select(starts_with("syn"))), ".delta")

abcd.c.2 <- inner_join(cbcl.deltas, abcd.pred)
p2 <- corr.table(abcd.c.2%>%select(predicted), abcd.c.2 %>% select(ends_with("delta")), method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval < 0.1, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.2)))
################################################################################
# brief problems monitor bpm
abcd.bpm <- read_csv(paste0(abcd.raw.dir, "/mental-health/mh_y_bpm.csv")) %>%
  select(IID = src_subject_id, eventname, 
         att = bpm_y_scr_attention_r, int = bpm_y_scr_internal_r,
         ext = bpm_y_scr_external_r, tot = bpm_y_scr_totalprob_r)
abcd.bpm <- inner_join(inner_join(abcd.bpm, abcd.meds), demo) %>%
  drop_na()
bpm.corrected <- cbind(abcd.bpm %>% select(-ends_with("t")), 
                       lapply(abcd.bpm %>% select(ends_with("t")), FUN = function(x) 
                         residuals(glm(y ~ interview_age + sex + interview_age:sex 
                                       + clonidine + guanfacine + dexmethylphenidate
                                       , data = abcd.bpm %>% mutate(y=x), family = poisson())))) %>%
  group_by(IID) %>%
  filter(n_distinct(methylphenidate)==2) %>%
  group_by(IID, methylphenidate) %>%
  mutate_at(.vars = vars(ends_with("t")), .funs = function(x) mean(x)) %>%
  ungroup() %>%
  distinct(IID, methylphenidate, .keep_all = T)
####################################
ps3 <- bpm.corrected %>% 
  pivot_longer(cols = ends_with("t"), names_to = "question", values_to = "score") %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  ggpaired(x = "methylphenidate", y = "score",  group = "methylphenidate",
           color = "methylphenidate", 
           line.color = "gray", line.size = 0.4, palette = "jco") +
  facet_wrap("question", scales = "free_y", ncol = 2) + theme_minimal() +
  stat_compare_means(paired = T)
#########################
# get deltas per participant
# being on MPH - not being on MPH
bpm.deltas <- do.call(cbind, apply(bpm.corrected %>% select(ends_with("t")),
                                   MARGIN = 2, FUN = function(x) bpm.corrected %>% mutate(y = x) %>%
                                     group_by(IID) %>%                                        
                                     summarise(delta = (y[methylphenidate== 1] - y[methylphenidate== 0])))) %>%
  select(IID = 1, ends_with("delta"))
abcd.c.5 <- inner_join(bpm.deltas, abcd.pred)

p5 <-  corr.table(abcd.c.5%>%select(predicted), abcd.c.5 %>% select(ends_with("delta")), method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval < 0.1, paste0("ρ: ", round(r, 3), ",  p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.5)))
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
                                       + clonidine + guanfacine + dexmethylphenidate
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
           line.color = "gray", line.size = 0.4, palette = "jco") +
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
# pgs
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", " "))

abcd.c.3 <- inner_join(abcd.pgs, sst.r1.deltas) %>%
  mutate_at(.vars = vars(colnames(abcd.pgs)[-1]), .funs = function(x) scale(x, scale = T, center = T)[,1])

p3 <- corr.table(abcd.c.3%>%select(colnames(abcd.pgs)[-1]), abcd.c.3 %>% select(ends_with("delta")), method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval < 0.05, paste0("ρ: ", round(r, 3), "\n","p: ", round(pval, 4)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.3)))
abcd.c.4 <- inner_join(abcd.pgs, cbcl.deltas) %>%
  mutate_at(.vars = vars(colnames(abcd.pgs)[-1]), .funs = function(x) scale(x, scale = T, center = T)[,1])
p4 <- corr.table(abcd.c.4%>%select(colnames(abcd.pgs)[-1]), abcd.c.4 %>% select(ends_with("delta")), method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval < 0.05, paste0("ρ: ", round(r, 3), "\n","p: ", round(pval, 4)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.4)))

abcd.c.6 <- inner_join(abcd.pgs, bpm.deltas) %>%
  mutate_at(.vars = vars(colnames(abcd.pgs)[-1]), .funs = function(x) scale(x, scale = T, center = T)[,1])

p6 <- corr.table(abcd.c.6%>%select(colnames(abcd.pgs)[-1]), abcd.c.6 %>% select(ends_with("delta")), method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval < 0.05, paste0("ρ: ", round(r, 3), "\n","p: ", round(pval, 4)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.6)))

abcd.c.8 <- inner_join(abcd.pgs, lmt.deltas) %>%
  mutate_at(.vars = vars(colnames(abcd.pgs)[-1]), .funs = function(x) scale(x, scale = T, center = T)[,1])

p8 <- corr.table(abcd.c.8%>%select(colnames(abcd.pgs)[-1]), abcd.c.8 %>% select(ends_with("delta")), method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval < 0.05, paste0("ρ: ", round(r, 3), "\n","p: ", round(pval, 4)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.8), "\n",
                                              "delta per question = score_on_MPH - score_off_MPH", "\n", 
                                              "scores were first derived from:", "\n",
                                              "residuals(glm(y ~ interview_age + sex + interview_age:sex + clonidine + guanfacine + dexmethylphenidate, family = poisson()))"))

ps4 <- corr.table(abcd.c.3%>%select(starts_with("ADHD")),
           abcd.c.3%>%select(starts_with("cog")), method = "spearman") %>%
  filter(V1 != V2) %>%
  ggplot(aes(x=V1, y=V2, fill=r, label=ifelse(pval<0.05, paste0("ρ: ", round(r, 3), ", ",
                                                                "p: ", round(pval, 5)), "")))+
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.3)))
ps5 <- corr.table(abcd.c.4%>%select(starts_with("ADHD")),
                  abcd.c.4%>%select(starts_with("cog")), method = "spearman") %>%
  filter(V1 != V2) %>%
  ggplot(aes(x=V1, y=V2, fill=r, label=ifelse(pval<0.05, paste0("ρ: ", round(r, 3), ", ",
                                                                "p: ", round(pval, 5)), "")))+
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.4)))
ps6 <- corr.table(abcd.c.6%>%select(starts_with("ADHD")),
                  abcd.c.6%>%select(starts_with("cog")), method = "spearman") %>%
  filter(V1 != V2) %>%
  ggplot(aes(x=V1, y=V2, fill=r, label=ifelse(pval<0.05, paste0("ρ: ", round(r, 3), ", ",
                                                                "p: ", round(pval, 5)), "")))+
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.6)))
ps8 <- corr.table(abcd.c.8%>%select(starts_with("ADHD")),
                  abcd.c.8%>%select(starts_with("cog")), method = "spearman") %>%
  filter(V1 != V2) %>%
  ggplot(aes(x=V1, y=V2, fill=r, label=ifelse(pval<0.05, paste0("ρ: ", round(r, 3), ", ",
                                                                "p: ", round(pval, 5)), "")))+
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c.8)))

################################################################################


################################################################################
patchwork::wrap_plots(ps1,ps3,ps2,ps7,nrow = 2, widths = c(4,2)) #sst-cbcl-bpm-lmt-mph-status
patchwork::wrap_plots(ps4+theme(axis.text.x = element_blank()),
                      ps5+theme(axis.text.x = element_blank()),
                      ps6+theme(axis.text.x = element_blank()),
                      ps8,nrow = 4) #sst-cbcl-bpm-lmt-pgs-corr
patchwork::wrap_plots(p1,p2,p5,p7) # sst-cbcl-bpm-lmt-delta-predicted-mph
patchwork::wrap_plots(p3+theme(axis.text.x = element_blank()),
                      p4+theme(axis.text.x = element_blank()),
                      p6+theme(axis.text.x = element_blank()),
                      p8, 
                      ncol = 1, heights = c(8,11,4,2)) #sst-cbcl-bpm-lmt-delta-pgs


