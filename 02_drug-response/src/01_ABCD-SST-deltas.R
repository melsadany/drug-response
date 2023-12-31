################################################################################
#                  check SST change by status of taking ADHD drugs             #
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
missing.samples <- read_csv("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/subjects-missing-med-name.csv") %>%
  mutate(drop = T)
abcd.meds.r <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/abcd5-meds-matrix.rds") %>%
  as.data.frame() %>%
  select(c(1:2), any_of(all.adhd.meds$drug)) %>%
  # filter(!(grepl("3", eventname) | grepl("4", eventname))) %>%
  mutate(methylphenidate = ifelse((methylphenidate+ritalin)>=1 & concerta == 0,1,0),
         amphetamine = ifelse((adderall+amphetamine)>=1,1,0),
         lisdexamfetamine = ifelse((vyvanse+lisdexamfetamine)>=1,1,0),
         stim = ifelse((methylphenidate+adderall+amphetamine+concerta+ritalin+vyvanse+lisdexamfetamine)>=1,1,0),
         non_stim = ifelse((intuniv+strattera+tenex+atomoxetine+clonidine+guanfacine)>=1,1,0),
         guanfacine = ifelse((guanfacine+tenex)>=1 & intuniv == 0,1,0),
         atomoxetine = ifelse((atomoxetine+strattera)>=1,1,0)) %>% 
  select(-c(ritalin, adderall, tenex, strattera,intuniv, vyvanse))
abcd.meds <- left_join(abcd.meds.r, missing.samples) %>%
  filter(is.na(drop)) %>%
  select(-drop)
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
# ABCD SST data -----------------------------------------------------------
sst.raw <- read_csv(paste0(abcd.raw.dir, "/imaging/mri_y_tfmr_sst_beh.csv"))
# only keeping 4 questions of interest for the first run of the task only
sst.r1 <- left_join(sst.raw %>%
                      filter(tfmri_sst_beh_glitchflag ==0) %>%
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
# SST deltas by drug status for participants -----------------------------
registerDoMC(cores = 3)
sst.meds.deltas <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d <- adhd.meds$drug[i]
  # print(i)
  # get number of participants that match this criteria
  # keep participants of 2 events and 2 distinct drug status at least
  n <- nrow(sst.all %>%
              rename(drug = which(colnames(sst.all)==d)) %>%
              group_by(IID) %>%
              filter(n_distinct(eventname) >= 2 & n_distinct(drug) > 1) %>%
              ungroup() %>%
              distinct(IID))
  # calculate the delta of each sst question
  # delta defines as (score_on_the_drug - score_off_the_drug)
  t <- sst.all %>%
    rename(drug = which(colnames(sst.all)==d)) %>%
    group_by(IID) %>%
    filter(n_distinct(eventname) >= 2 & n_distinct(drug) > 1) %>%
    ungroup() %>%
    # some participants have multiple sst datapoints for the same drug status
    # get the average by drug status after being grouped based on IID and drug binary status
    group_by(IID, drug) %>%
    mutate_at(.vars = vars(starts_with("e_")), .funs = function(x) mean(x)) %>%
    ungroup() %>%
    # make sure to keep one data point per distinct IID, drug status
    distinct(IID, drug, .keep_all = T) %>%
    pivot_longer(cols = c(starts_with("e_")), 
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
pgs.predicted.deltas <- inner_join(sst.meds.deltas %>% pivot_wider(names_from = "question", values_from = "delta", 
                                                                   id_cols = c(IID,sex,drug, n_samples)), 
                                   inner_join(abcd.pgs, abcd.pred)) %>% 
  mutate(group = ifelse(cog_gFactor > 0 & `ADHD-Demontis` > 0, 1, 
                        ifelse(cog_gFactor > 0 & `ADHD-Demontis` < 0, 2, 
                               ifelse(cog_gFactor < 0 & `ADHD-Demontis` > 0, 4, 
                                      3)))) %>%
  filter(group == 4)

registerDoMC(cores = 3)
drug.deltas.rsquared <- foreach (i = 1:length(unique(pgs.predicted.deltas$drug)), .combine = rbind) %dopar% {
  dname <- unique(pgs.predicted.deltas$drug)[i]
  nsamples <- nrow(pgs.predicted.deltas %>% filter(drug == dname))
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
    # changed the model based on the wanted figure
    mm <- m12
    # mm <- m3
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
    mutate(drug = dname, n = nsamples)
  return(df2)
}
# write_csv(drug.deltas.rsquared, "data/figs-data/predicting-sst-deltas-m12-012.csv")
p1 <- drug.deltas.rsquared %>%
  filter(!drug %in% c("lisdexamfetamine", "concerta", "atomoxetine"), !grepl("correct_go", question)) %>%
  mutate(drug = factor(drug, levels = c("clonidine", "guanfacine", "atomoxetine", 
                                        "amphetamine", "methylphenidate", "concerta", "lisdexamfetamine",
                                        "non_stim", "stim"))) %>%
  mutate(question = sub("\\.[0-9]", "", question)) %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  mutate(m_sig = ifelse(fdr<0.05, "FDR < 0.05", "FDR \u2265 0.05")) %>%
  mutate(question = sub("e_as_", "", question)) %>%
  mutate(var = factor(var, levels = c("ADHD-Demontis", "cog_gFactor", "predicted", "cog_gFactor:predicted"))) %>%
  ggplot(aes(x=Estimate, y = var)) +
  # ggplot(aes(x=Estimate, y = question)) +
  geom_point(aes(alpha = sig, color = var, shape = m_sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name = "") +
  scale_shape_manual(values = c("FDR < 0.05" = 8, "FDR \u2265 0.05" = 1), name = "") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig, color = var), 
                 linewidth = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("ADHD-Demontis" = six.colors[1], 
                                "cog_gFactor" = six.colors[4], 
                                "predicted" = six.colors[3]), name = "") +
  facet_grid2(rows = vars(question), cols = vars(drug), remove_labels = "y")+
  # facet_grid2(cols = vars(drug), remove_labels = "y")+
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0)+
  annotate("segment", x=-Inf, xend=-Inf, y=0, yend=0) +
  theme(
    strip.text.y.right = element_text(angle = 0),
    # strip.text.y.right = element_blank(),
    # strip.text.y.left = element_blank(),
    # axis.text.y.left = element_blank(), 
    axis.line.x.bottom = element_line(colour = "black")
    # title = element_text(colour = six.colors[4])
    ) +
  labs(x = "Estimate", 
       y="")
p2<- drug.deltas.rsquared %>%
  filter(!drug %in% c("lisdexamfetamine", "concerta", "atomoxetine")) %>%
  mutate(drug = factor(drug, levels = c("clonidine", "guanfacine", "atomoxetine", 
                                        "amphetamine", "methylphenidate", "concerta", "lisdexamfetamine",
                                        "non_stim", "stim"))) %>%
  distinct(drug, n) %>%
  ggplot(aes(x=drug, y=n, fill = drug, label = n))+
  geom_bar(stat = "identity", show.legend = F, width = 0.3) +
  geom_text(size=3, nudge_y = -10) +
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5), 
        axis.text.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.line.y.left = element_blank()) +
  labs(x="", y="n(samples)")
patchwork::wrap_plots(p1,p2, heights = c(6,1)) # predicting SST deltas
# hh <- patchwork::wrap_plots(p1+labs(title = "HIGH cognitive PGS & HIGH ADHD PGS"),p2, heights = c(4,1))
# hh3 <- patchwork::wrap_plots(p1+labs(title = "HIGH cognitive PGS & HIGH ADHD PGS"),p2, heights = c(2,1))
####
# try this
# patchwork::wrap_plots(hh1, hh, ncol = 1, heights = c(1,8))
patchwork::wrap_plots(hl,hh,ll,lh, ncol = 2, nrow = 2)
patchwork::wrap_plots(hl3,hh3,ll3,lh3, ncol = 2, nrow = 2)

####






####
# supplementary figs ------------------------------------------------------
# histogram for PGS distribution
inner_join(sst.meds.deltas %>% 
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
       title = "distribution of PGS for SST samples", 
       caption = paste0("n(samples): ", nrow(sst.meds.deltas %>% filter(drug == "methylphenidate") %>% distinct(IID))))
# deltas boxplot for distribution
sst.meds.deltas %>% 
  filter(drug == "methylphenidate", grepl("_as_", question)) %>%
  mutate(question = sub("e_", "", question)) %>%
  ggplot(aes(x=question, y=delta, fill = question))+
  geom_boxplot(outlier.size = 0.3, show.legend = F, width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  labs(y="score delta")+
  labs(title = "distribution of SST deltas for MPH samples", x="",
       caption = paste0("n(samples): ", nrow(sst.meds.deltas %>% filter(drug=="methylphenidate")%>% distinct(IID))))
# on and off mph fig
p1 <- sst.meds.deltas %>% 
  filter(grepl("_as_", question), !drug %in% c("atomoxetine", "lisdexamfetamine", "concerta")) %>%
  filter(!grepl("correct_go", question)) %>%
  mutate(question = sub("e_as_", "", question)) %>%
  ggplot(aes(x=`0`, y=`1`)) +
  geom_point(size=0.3)+
  facet_grid2(rows = vars(question), cols = vars(drug))+
  stat_cor(method = "spearman", size = 2.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  theme(strip.text.y.right = element_text(angle = 0))+
  labs(x="score when drug status = 0", y="score when drug status = 1",
       caption = "scores are corrected for age, sex, and their interaction")
p2 <- sst.meds.deltas %>%
  filter(grepl("_as_", question), !drug %in% c("atomoxetine", "lisdexamfetamine", "concerta")) %>%
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
# supplementary tables ----------------------------------------------------
# samples included per drug
sst.mph.demo <- inner_join(demo, 
                            inner_join(sst.meds.deltas%>%distinct(IID,drug),
                                       sst.all%>%select(IID, eventname))) %>%
  filter(IID %in% abcd.pgs$IID)
write_csv(inner_join(sst.mph.demo %>% 
                       group_by(sex,drug) %>%
                       summarise(avg = mean(interview_age), 
                                 sd = sd(interview_age), 
                                 min = min(interview_age), 
                                 max = max(interview_age)),
                     sst.mph.demo %>% 
                       distinct(IID ,sex, drug)%>% 
                       group_by(sex, drug) %>%
                       summarise(count = n())) %>%
            mutate(measure = "SST", data = "ABCD"),
          file = "figs/paper/tmp/abcd-sst-data-stats.csv")
####
