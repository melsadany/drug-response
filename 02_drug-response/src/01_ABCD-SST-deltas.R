################################################################################
#                  check SST change by status of taking ADHD drugs             #
################################################################################
# packages setup ----------------------------------------------------------
# there's an "Extras" section at the end of this file.
# Make sure to run it after loading the needed packages
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
                                     "dextroamphetamine",  "ritalin", "intuniv", "strattera",
                                     "tenex", "amphetamine", "dexmethylphenidate", "lisdexamfetamine",
                                     "atomoxetine", "clonidine", "guanfacine"
))
adhd.meds <- data.frame(drug = c("methylphenidate",
                                 "concerta",
                                 "stim", "non_stim",
                                 # "dextroamphetamine", "dexmethylphenidate", 
                                 "amphetamine", "lisdexamfetamine",
                                 "atomoxetine", "clonidine", "guanfacine"
))
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/abcd5-meds-matrix.rds") %>%
  as.data.frame() %>%
  select(c(1:2), any_of(all.adhd.meds$drug)) %>%
  # mutate(methylphenidate = ifelse((methylphenidate+ritalin+concerta)>=1,1,0),
  mutate(methylphenidate = ifelse((methylphenidate+ritalin)>=1,1,0),
         amphetamine = ifelse((adderall+amphetamine)>=1,1,0),
         lisdexamfetamine = ifelse((vyvanse+lisdexamfetamine)>=1,1,0),
         stim = ifelse((methylphenidate+adderall+amphetamine+concerta+ritalin)>=1,1,0),
         non_stim = ifelse((vyvanse+intuniv+strattera+tenex+lisdexamfetamine+atomoxetine+clonidine+guanfacine)>=1,1,0),
         # guanfacine = ifelse((guanfacine+tenex+intuniv)>=1,1,0),
         guanfacine = ifelse((guanfacine+tenex)>=1,1,0),
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
  # filter(IID %in% sst.r1.deltas$IID) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
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
# correlation between taking drug and SST --------------------------------
corr.table(sst.all %>% select(adhd.meds$drug),
           sst.all %>% select(starts_with("e_")),
           method = "spearman") %>%
  filter(V1 %in% adhd.meds$drug, ! V2 %in% adhd.meds$drug) %>%
  mutate(value_type = factor(ifelse(grepl("raw", V2), "raw data", 
                                    "corrected for age, sex, and interaction"), 
                             levels = c("raw data", 
                                        "corrected for age, sex, and interaction"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", nrow(sst.all), "\n",
                        "corrected for:", "\n",
                        "\tinterview_age + sex + interview_age:sex", "\n",
                        "* pval < 0.05 & not FDR sig", "\n", 
                        "** pval < 0.01 & not FDR sig", "\n", 
                        "*** FDR < 0.05"), 
       title = "correlation of taking a medication with sst scores")
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
# SST deltas average by drug ---------------------------------------------
sst.meds.deltas %>%
  group_by(question, drug) %>%
  summarise_at(.vars = vars(delta, n_samples), 
               .funs = function(x) mean(x)) %>%
  mutate(value_type = factor(ifelse(grepl("raw", question), 
                                    "raw data", "corrected for age, sex, and interaction"),
                             levels = c("raw data", 
                                        "corrected for age, sex, and interaction"))) %>%
  ggplot(aes(x=drug, y=question, fill = delta)) +
  geom_tile()+
  my.guides+null_labs+
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], 
                       name = "average delta (on_drug - off_drug)")+
  labs(caption = paste0("n(samples): ", "\n\t",
                        paste(apply(sst.meds.deltas%>%ungroup()%>%select(drug, n_samples)%>%distinct(), 
                                    1, function(x) paste(x[1], ":", x[2])), collapse = "\n\t")), 
       title = "average deltas for sst performance measures per drug")
# sst deltas per drug in a boxplot for raw sst
sst.meds.deltas %>%
  filter(!grepl("totprob", question)) %>%
  mutate(value_type = factor(ifelse(grepl("raw", question), 
                                    "raw data", "corrected for age, sex, and interaction"),
                             levels = c("raw data", 
                                        "corrected for age, sex, and interaction"))) %>%
  filter(grepl("raw", question)) %>%
  ggplot(aes(x=question, y=delta, fill=question))+
  geom_boxplot(outlier.size = 0.5)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  my.guides+null_labs+ labs(y="score delta")+
  facet_wrap("drug", scales = "free_y", ncol = 3)+
  labs(caption = paste0("n(samples): ", "\n\t",
                        paste(apply(sst.meds.deltas%>%ungroup()%>%select(drug, n_samples)%>%distinct(), 
                                    1, function(x) paste(x[1], ":", x[2])), collapse = "\n\t")), 
       title = "sst raw score change (delta) per drug")
# sst deltas per drug in a boxplot for age-sex corrected sst
sst.meds.deltas %>%
  filter(!grepl("totprob", question)) %>%
  mutate(value_type = factor(ifelse(grepl("raw", question), 
                                    "raw data", "corrected for age, sex, and interaction"),
                             levels = c("raw data", 
                                        "corrected for age, sex, and interaction"))) %>%
  filter(grepl("as", question)) %>%
  ggplot(aes(x=question, y=delta, fill=question))+
  geom_boxplot(outlier.size = 0.5)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  my.guides+null_labs+ labs(y="score delta")+
  facet_wrap("drug", scales = "free_y", ncol = 3)+
  labs(caption = paste0("n(samples): ", "\n\t",
                        paste(apply(sst.meds.deltas%>%ungroup()%>%select(drug, n_samples)%>%distinct(), 
                                    1, function(x) paste(x[1], ":", x[2])), collapse = "\n\t")), 
       title = "sst age-sex corrected score change (delta) per drug")
####
# SST deltas per drug with PGS -------------------------------------------
sst.meds.deltas.pgs <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d <- adhd.meds$drug[i]
  t <- inner_join(sst.meds.deltas %>% 
                    filter(drug == d) %>% 
                    pivot_wider(names_from = question, 
                                values_from = delta, 
                                id_cols = c(IID, sex)),
                  abcd.pgs)
  ret <- corr.table(t %>% select(colnames(abcd.pgs)[-1]),
                    t %>% select(starts_with("e_")),
                    method = "spearman") %>%
    filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
    mutate(value_type = factor(ifelse(grepl("raw_", V2), 
                                      "raw data", 
                                      "corrected for age, sex, and interaction"), 
                               levels = c("raw data", 
                                          "corrected for age, sex, and interaction"))) %>%
    group_by(value_type) %>%
    mutate(FDR = p.adjust(pval, method = "fdr")) %>%
    mutate(drug = d) %>%
    mutate(n_samples = nrow(t))
  return(ret)
}
sst.meds.deltas.pgs %>%
  filter(drug == "methylphenidate") %>%
  # filter(drug == "stim") %>%
  # filter(drug == "atomoxetine") %>%
  filter(grepl("as", V2)) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), cols = vars(drug), scales = "free_y") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", "\n\t",
                        # paste(apply(sst.meds.deltas.pgs%>%ungroup()%>%select(drug, n_samples)%>%distinct(), 
                        # "methylphenidate: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "methylphenidate/ritalin/concerta: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        "methylphenidate/ritalin: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "concerta: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="concerta")%>%distinct(n_samples), "\n",
                        # "guanfacine/tenex/intuniv: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="guanfacine")%>%distinct(n_samples), "\n",
                        # "guanfacine/tenex: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="guanfacine")%>%distinct(n_samples), "\n",
                        # "clonidine: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="clonidine")%>%distinct(n_samples), "\n",
                        # "atomoxetine/strattera: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="atomoxetine")%>%distinct(n_samples), "\n",
                        # "lisdexamfetamine/vyvanse: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="lisdexamfetamine")%>%distinct(n_samples), "\n",
                        # "amphetamine/adderall: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="amphetamine")%>%distinct(n_samples), "\n",
                        # "stim = methylphenidate/ritalin/concerta/amphetamine/adderall: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="stim")%>%distinct(n_samples), "\n",
                        # "non_stim = vyvanse/intuniv/strattera/tenex/lisdexamfetamine/atomoxetine/clonidine/guanfacine: ", sst.meds.deltas.pgs %>%ungroup()%>% filter(drug=="non_stim")%>%distinct(n_samples), "\n",
                        "* pval < 0.05 & not FDR sig", "\n", 
                        "** pval < 0.01 & not FDR sig", "\n", 
                        "*** FDR < 0.05"), 
       # title = "correlation of sst score change (delta) per drug with PGS")
       title = "correlation of sst score change (delta) with PGS")
####
# SST deltas of MPH and predicted MPH response ---------------------------
sst.pred <- inner_join(sst.meds.deltas %>% 
                          filter(drug == "methylphenidate") %>% 
                          pivot_wider(names_from = question, 
                                      values_from = delta, 
                                      id_cols = c(IID, sex)),
                        abcd.pred)
corr.table(sst.pred %>% select(predicted),
           sst.pred %>% select(starts_with("e_")),
           method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", V2), 
                                    "raw data", 
                                    "corrected for age, sex, and interaction"), 
                             levels = c("raw data", 
                                        "corrected for age, sex, and interaction"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  # filter(grepl("as", V2)) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, paste0(rho, ": ", round(r, 3),
                                                                    ", p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), scales = "free_y") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", nrow(sst.pred), "\n",
                        "\tmethylphenidate"), 
       title = "correlation of sst score change (delta) of MPH with predicted MPH response") + 
  theme(plot.title = element_text(hjust = 1, size = 9),axis.text.x.bottom = element_text(angle = 0))
####


####
# supplementary figs ------------------------------------------------------
# scatterplot for SST performance on and off MPH
sst.meds.deltas %>%
  filter(grepl("raw", question), drug == "methylphenidate") %>%
  mutate(question = sub("e_","", question)) %>%
  ggplot(aes(x=`0`, y=`1`)) +
  geom_point(size=0.6)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap("question", scales = "free")+
  labs(x="off-MPH", y="on-MPH", 
       caption = paste0("n(samples): ", nrow(sst.meds.deltas %>% filter(drug == "methylphenidate") %>% distinct(IID)), "\n",
                        "the red dashed line represents a line with a slope of 1 and intercepts with the origin (0,0)"))
# histogram for PGS distribution
inner_join(sst.meds.deltas %>% 
             filter(drug == "methylphenidate") %>% 
             pivot_wider(names_from = question, 
                         values_from = delta, 
                         id_cols = c(IID, sex)),
           abcd.pgs) %>%
  pivot_longer(colnames(abcd.pgs)[-1], names_to = "PGS", values_to = "score") %>%
  ggplot(aes(x=score))+
  geom_histogram()+
  facet_wrap("PGS", scales = "free")+
  labs(x="PGS", 
       title = "distribution of PGS for SST samples", 
       caption = paste0("n(samples): ", nrow(sst.meds.deltas %>% filter(drug == "methylphenidate") %>% distinct(IID))))
# histogram for predicted response distribution for SST samples
sst.pred %>% 
  ggplot(aes(x=predicted)) +
  geom_histogram()+
  labs(title = "distribution of predicted MPH response for SST samples", 
       caption = paste0("n(samples): ", nrow(sst.pred  %>% distinct(IID))))
# deltas boxplot for distribution
sst.meds.deltas %>% 
  filter(drug == "methylphenidate", grepl("raw", question)) %>%
  mutate(question = sub("e_", "", question)) %>%
  ggplot(aes(x=question, y=delta, fill = question))+
  geom_boxplot(outlier.size = 0.3, show.legend = F) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  labs(y="score delta")+
  labs(title = "distribution of SST deltas for MPH samples", x="",
       caption = paste0("n(samples): ", nrow(sst.meds.deltas %>% filter(drug=="methylphenidate")%>% distinct(IID))))
# correlation between predicted MPH response and PGS
pgs.pred <- inner_join(abcd.pred, abcd.pgs)
# pgs.pred <- inner_join(inner_join(abcd.pred, abcd.pgs), sst.meds.deltas %>% select(IID, drug)%>%filter(drug=="methylphenidate"))%>%distinct(IID ,.keep_all = T)
pgs.pred %>% 
  pivot_longer(cols = colnames(abcd.pgs)[-1], names_to = "PGS", values_to = "score") %>%
  ggplot(aes(x=predicted, y=score))+
  geom_point(size=0.3)+
  geom_smooth(method = "glm") +
  stat_cor() +
  facet_wrap("PGS", scales = "free_y") +
  labs(caption = paste0("n(samples): ", nrow(pgs.pred)
                        ,"\n\tmethylphenidate/ritalin"),
       title = "correlation between predicted MPH response and PGS")
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
# Extras ------------------------------------------------------------------
as.cap <- paste0("delta per question = score_on_the_drug - score_off_the_drug", "\n", 
                 "corrected for:", "\n",
                 "\tinterview_age + sex + interview_age:sex")
####
