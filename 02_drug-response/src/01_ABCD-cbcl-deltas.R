################################################################################
#                  check CBCL change by status of taking ADHD drugs            #
################################################################################
# packages setup ----------------------------------------------------------
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(ggpubr);library(ggExtra);library(ggh4x)
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
                                 "ritalin", "intuniv", "strattera",
                                 "tenex", "amphetamine", "lisdexamfetamine",
                                 "atomoxetine", "clonidine", "guanfacine"
))
adhd.meds <- data.frame(drug = c("methylphenidate",
                                 "concerta",
                                 "stim", "non_stim",
                                 "amphetamine", "lisdexamfetamine",
                                 "atomoxetine", "clonidine", "guanfacine"
))
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/abcd5-meds-matrix.rds") %>%
  as.data.frame() %>%
  select(c(1:2), any_of(all.adhd.meds$drug)) %>%
  filter(!(grepl("3", eventname) | grepl("4", eventname))) %>%
  mutate(methylphenidate = ifelse((methylphenidate+ritalin)>=1 & concerta == 0,1,0),
         amphetamine = ifelse((adderall+amphetamine)>=1,1,0),
         lisdexamfetamine = ifelse((vyvanse+lisdexamfetamine)>=1,1,0),
         stim = ifelse((methylphenidate+adderall+amphetamine+concerta+ritalin+vyvanse+lisdexamfetamine)>=1,1,0),
         non_stim = ifelse((intuniv+strattera+tenex+atomoxetine+clonidine+guanfacine)>=1,1,0),
         guanfacine = ifelse((guanfacine+tenex)>=1 & intuniv == 0,1,0),
         atomoxetine = ifelse((atomoxetine+strattera)>=1,1,0)) %>% 
  select(-c(ritalin, adderall, tenex, strattera,intuniv, vyvanse))
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
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
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
# CBCL deltas by drug status for participants -----------------------------
registerDoMC(cores = 3)
cbcl.meds.deltas <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d <- adhd.meds$drug[i]
  # print(i)
  # get number of participants that match this criteria
  # keep participants of 2 events and 2 distinct drug status at least
  n <- nrow(cbcl.all %>%
              rename(drug = which(colnames(cbcl.all)==d)) %>%
              group_by(IID) %>%
              filter(n_distinct(eventname) >= 2 & n_distinct(drug) > 1) %>%
              ungroup() %>%
              distinct(IID))
  # calculate the delta of each cbcl subscale
  # delta defines as (score_on_the_drug - score_off_the_drug)
  t <- cbcl.all %>%
    rename(drug = which(colnames(cbcl.all)==d)) %>%
    group_by(IID) %>%
    filter(n_distinct(eventname) >= 2 & n_distinct(drug) > 1) %>%
    ungroup() %>%
    # some participants have multiple cbcl datapoints for the same drug status
    # get the average by drug status after being grouped based on IID and drug binary status
    group_by(IID, drug) %>%
    mutate_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), .funs = function(x) mean(x)) %>%
    ungroup() %>%
    # make sure to keep one data point per distinct IID, drug status
    distinct(IID, drug, .keep_all = T) %>%
    pivot_longer(cols = c(starts_with("syn"), starts_with("dsm5")), 
                 names_to = "question", values_to = "val") %>%
    arrange(IID, question, drug) %>%
    mutate(drug = as.factor(drug)) %>%
    pivot_wider(names_from = drug, values_from = val, 
                id_cols = c(IID, sex, question)) %>%
    # calculate delta
    mutate(delta = `1` - `0`) %>% 
    mutate(drug = d) %>%
    mutate(n_samples = n)
  return(t)
}
####
# models for predicting deltas by PGS, predicted, or combination ----------
pgs.predicted.deltas <- inner_join(cbcl.meds.deltas %>% pivot_wider(names_from = "question", values_from = "delta", id_cols = c(IID,sex,drug, n_samples)), 
                                   inner_join(abcd.pgs, abcd.pred))
registerDoMC(cores = 3)
drug.deltas.rsquared <- foreach (i = 1:length(unique(pgs.predicted.deltas$drug)), .combine = rbind) %dopar% {
  dname <- unique(pgs.predicted.deltas$drug)[i]
  nsamples <- pgs.predicted.deltas %>% filter(drug == dname) %>% distinct(n_samples)
  models.r.squared <- do.call(rbind, lapply(pgs.predicted.deltas %>% filter(drug == dname) %>% select(contains("_as_")), function(x) {
    m1 <- lm(y ~ `ADHD-Demontis`, data = pgs.predicted.deltas %>% filter(drug == dname) %>% mutate(y=x))
    m2 <- lm(y ~ cog_gFactor, data = pgs.predicted.deltas %>% filter(drug == dname) %>% mutate(y=x))
    m3 <- lm(y ~ predicted, data = pgs.predicted.deltas %>% filter(drug == dname) %>% mutate(y=x))
    m12 <- lm(y ~ `ADHD-Demontis` + cog_gFactor, data = pgs.predicted.deltas %>% filter(drug == dname) %>% mutate(y=x))
    m13 <- lm(y ~ `ADHD-Demontis` + predicted, data = pgs.predicted.deltas %>% filter(drug == dname) %>% mutate(y=x))
    m23 <- lm(y ~ cog_gFactor + predicted, data = pgs.predicted.deltas %>% filter(drug == dname) %>% mutate(y=x))
    m232 <- lm(y ~ cog_gFactor + predicted + cog_gFactor:predicted, data = pgs.predicted.deltas %>% filter(drug == dname) %>% mutate(y=x))
    m123 <- lm(y ~ `ADHD-Demontis` + cog_gFactor + predicted, data = pgs.predicted.deltas %>% filter(drug == dname) %>% mutate(y=x))
    mall <- lm(y ~ `ADHD-Demontis` + cog_gFactor + predicted + cog_gFactor:predicted , data = pgs.predicted.deltas %>% filter(drug == dname) %>% mutate(y=x))
    df <- data.frame(m1 = summary(m1)$adj.r.squared,
                     m2 = summary(m2)$adj.r.squared,
                     m3 = summary(m3)$adj.r.squared,
                     m12 = summary(m12)$adj.r.squared,
                     m13 = summary(m13)$adj.r.squared,
                     m23 = summary(m23)$adj.r.squared,
                     m232 = summary(m232)$adj.r.squared,
                     m123 = summary(m123)$adj.r.squared)
    # changed the model from mall to m123
    mm <- m12
    df3 <- summary(mm)$coefficients %>% 
      as.data.frame() %>% 
      rownames_to_column("var") %>% 
      mutate(var = str_remove_all(pattern = "`",string = var), 
             confint_min = confint(mm)[,1], 
             confint_max = confint(mm)[,2]) %>% 
      filter(!grepl("Intercept", var)) %>%
      rename(pval = 5)
    # return(df)
    return(df3)
  }))
  df2 <- models.r.squared %>%
    rownames_to_column("question") %>%
    mutate(drug = dname, n = nsamples$n_samples)
  return(df2)
}
# write_csv(drug.deltas.rsquared, "data/figs-data/predicting-cbcl-deltas-m12-012.csv")
p1 <- drug.deltas.rsquared %>%
  filter(!drug %in% c("lisdexamfetamine", "concerta", "atomoxetine")) %>%
  mutate(question = sub("\\.[0-9]", "", question)) %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  mutate(m_sig = ifelse(pval<0.08, "pval < 0.1", "pval \u2265 0.1")) %>%
  mutate(question = sub("as_", "", question)) %>%
  mutate(var = factor(var, levels = c("ADHD-Demontis", "cog_gFactor", "predicted", "cog_gFactor:predicted"))) %>%
  ggplot(aes(x=Estimate, y = var)) +
  # ggplot(aes(x=Estimate, y = question)) +
  geom_point(aes(alpha = sig, color = var, shape = m_sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.3)) +
  scale_shape_manual(values = c(8,1)) +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig, color = var), 
                 linewidth = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = six.colors[c(1,4)]) +
  facet_grid2(rows = vars(question), cols = vars(drug), remove_labels = "y")+
  # facet_grid2(cols = vars(drug), remove_labels = "y")+
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0)+
  annotate("segment", x=-Inf, xend=-Inf, y=0, yend=0) +
  theme(strip.text.y.right = element_text(angle = 0),
        axis.line.x.bottom = element_line(colour = "black")) +
  labs(x = "Estimate", 
       y="")
p2<- drug.deltas.rsquared %>%
  filter(!drug %in% c("lisdexamfetamine", "concerta", "atomoxetine")) %>%
  distinct(drug, n) %>%
  ggplot(aes(x=drug, y=n, fill = drug, label = n))+
  geom_bar(stat = "identity", show.legend = F, width = 0.3) +
  geom_text(size=3, nudge_y = -10) +
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5), 
        axis.text.y.left = element_blank()) +
  labs(x="", y="n(samples)")
patchwork::wrap_plots(p1,p2, heights = c(12,1)) # predicting CBCL deltas
####
# supplementary figures ---------------------------------------------------
# histogram for PGS distribution
inner_join(cbcl.meds.deltas %>% 
             filter(drug == "methylphenidate") %>% 
             pivot_wider(names_from = question, 
                         values_from = delta, 
                         id_cols = c(IID, sex)),
           abcd.pgs) %>%
  pivot_longer(colnames(abcd.pgs)[-1], names_to = "PGS", values_to = "score") %>%
  filter(grepl("ADHD", PGS) | grepl("gFa", PGS)) %>%
  ggplot(aes(x=score))+
  geom_histogram()+
  facet_wrap("PGS", scales = "free")+
  labs(x="PGS", 
       title = "distribution of PGS for MPH samples in CBCL", 
       caption = paste0("n(samples): ", nrow(cbcl.meds.deltas %>% filter(drug == "methylphenidate") %>% distinct(IID))))
# deltas boxplot for distribution
cbcl.meds.deltas %>% 
  filter(drug == "methylphenidate", grepl("_as_", question), !grepl("totprob", question)) %>%
  ggplot(aes(x=question, y=delta, fill = question))+
  geom_boxplot(outlier.size = 0.3, show.legend = F, width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  labs(y="score delta")+
  labs(title = "distribution of CBCL deltas for MPH samples", x="",
       caption = paste0("n(samples): ", nrow(cbcl.meds.deltas %>% filter(drug=="methylphenidate")%>% distinct(IID))))
# on and off mph fig
p1 <- cbcl.meds.deltas %>% 
  filter(grepl("_as_", question), !grepl("totp", question), !drug %in% c("atomoxetine", "lisdexamfetamine", "concerta")) %>%
  filter(grepl("adhd", question) | grepl("attention", question)) %>%
  mutate(question = sub("_as", "", question)) %>%
  ggplot(aes(x=`0`, y=`1`)) +
  geom_point(size=0.3)+
  # geom_smooth(method = "lm", se = F)+
  facet_grid2(rows = vars(question), cols = vars(drug))+
  stat_cor(method = "spearman", size = 2.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  theme(strip.text.y.right = element_text(angle = 0))+
  labs(x="score when drug status = 0", y="score when drug status = 1",
       caption = "scores are corrected for age, sex, and their interaction")
p2 <- cbcl.meds.deltas %>%
  filter(grepl("_as_", question), !grepl("totp", question), !drug %in% c("atomoxetine", "lisdexamfetamine", "concerta")) %>%
  mutate(question = sub("_as", "", question)) %>%
  group_by(question, drug) %>% 
  summarise(count = n()) %>%
  ungroup() %>%
  distinct(drug, count) %>%
  ggplot(aes(x=drug, y=count, fill = drug, label = count)) +
  geom_bar(stat = "identity", show.legend = F, width = 0.3) +
  geom_text(size=3, nudge_y = -10) +
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5),
        axis.text.y.left = element_blank()) +
  labs(x="", y="n(samples)")

patchwork::wrap_plots(p1,p2,heights = c(3,1))
####
####
# supplementary tables ----------------------------------------------------
# samples included per drug
cbcl.mph.demo <- inner_join(demo, 
                            inner_join(cbcl.meds.deltas%>%distinct(IID,drug),
                                       cbcl.all%>%select(IID, eventname))) %>%
  filter(IID %in% abcd.pgs$IID)
write_csv(inner_join(cbcl.mph.demo %>% 
                       group_by(sex,drug) %>%
                       summarise(avg = mean(interview_age), 
                                 sd = sd(interview_age), 
                                 min = min(interview_age), 
                                 max = max(interview_age)),
                     cbcl.mph.demo %>% 
                       distinct(IID ,sex, drug)%>% 
                       group_by(sex, drug) %>%
                       summarise(count = n())) %>%
            mutate(measure = "CBCL", data = "ABCD"),
          file = "figs/paper/tmp/abcd-cbcl-data-stats.csv")
####
