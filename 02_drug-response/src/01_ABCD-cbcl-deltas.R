################################################################################
#                  check CBCL change by status of taking ADHD drugs            #
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
                                 "dextroamphetamine", 
                                 "amphetamine", "dexmethylphenidate", "lisdexamfetamine",
                                 "atomoxetine", "clonidine", "guanfacine"
))
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/abcd5/abcd5-meds-matrix.rds") %>%
  as.data.frame() %>%
  select(c(1:2), any_of(all.adhd.meds$drug)) %>%
  # mutate(methylphenidate = ifelse((methylphenidate+ritalin+concerta)>=1,1,0),
  mutate(methylphenidate = ifelse((methylphenidate+ritalin)>=1,1,0),
         amphetamine = ifelse((adderall+amphetamine)>=1,1,0),
         lisdexamfetamine = ifelse((vyvanse+lisdexamfetamine)>=1,1,0),
         guanfacine = ifelse((guanfacine+tenex+intuniv)>=1,1,0),
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
# correlation between taking drug and CBCL --------------------------------
corr.table(cbcl.all %>% select(adhd.meds$drug),
           cbcl.all %>% select(starts_with("syn"), starts_with("dsm5")),
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
  labs(caption = paste0("n(samples): ", nrow(cbcl.all), "\n",
                        "corrected for:", "\n",
                        "\tinterview_age + sex + interview_age:sex", "\n",
                        "* pval < 0.05 & not FDR sig", "\n", 
                        "** pval < 0.01 & not FDR sig", "\n", 
                        "*** FDR < 0.05"), 
       title = "correlation of taking a medication with cbcl scores")
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
# CBCL deltas average by drug ---------------------------------------------
cbcl.meds.deltas %>%
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
                        paste(apply(cbcl.meds.deltas%>%ungroup()%>%select(drug, n_samples)%>%distinct(), 
                                    1, function(x) paste(x[1], ":", x[2])), collapse = "\n\t")), 
       title = "average deltas for cbcl problems per drug")

# cbcl deltas per drug in a boxplot for raw cbcl
cbcl.meds.deltas %>%
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
                        paste(apply(cbcl.meds.deltas%>%ungroup()%>%select(drug, n_samples)%>%distinct(), 
                                    1, function(x) paste(x[1], ":", x[2])), collapse = "\n\t")), 
       title = "cbcl raw score change (delta) per drug")
# cbcl deltas per drug in a boxplot for age-sex corrected cbcl
cbcl.meds.deltas %>%
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
                        paste(apply(cbcl.meds.deltas%>%ungroup()%>%select(drug, n_samples)%>%distinct(), 
                                    1, function(x) paste(x[1], ":", x[2])), collapse = "\n\t")), 
       title = "cbcl age-sex corrected score change (delta) per drug")
####
# CBCL deltas per drug with PGS -------------------------------------------
cbcl.meds.deltas.pgs <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d <- adhd.meds$drug[i]
  t <- inner_join(cbcl.meds.deltas %>% 
                    filter(drug == d) %>% 
                    pivot_wider(names_from = question, 
                                values_from = delta, 
                                id_cols = c(IID, sex)),
                  abcd.pgs)
  ret <- corr.table(t %>% select(colnames(abcd.pgs)[-1]),
                    t %>% select(starts_with("syn"), starts_with("dsm5")),
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
cbcl.meds.deltas.pgs %>%
  filter(drug == "methylphenidate") %>%
  filter(grepl("as", V2)) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), cols = vars(drug), scales = "free_y") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", "\n\t",
                        # paste(apply(cbcl.meds.deltas.pgs%>%ungroup()%>%select(drug, n_samples)%>%distinct(), 
                        "methylphenidate: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "methylphenidate/ritalin/concerta: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "methylphenidate/ritalin: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "concerta: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="concerta")%>%distinct(n_samples), "\n",
                        # "guanfacine/tenex/intuniv: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="guanfacine")%>%distinct(n_samples), "\n",
                        # "amphetamine/adderall: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="amphetamine")%>%distinct(n_samples), "\n",
                        "* pval < 0.05 & not FDR sig", "\n", 
                        "** pval < 0.01 & not FDR sig", "\n", 
                        "*** FDR < 0.05"), 
       # title = "correlation of cbcl score change (delta) per drug with PGS")
       title = "correlation of cbcl score change (delta) with PGS")
####
# CBCL deltas of MPH and predicted MPH response ---------------------------
cbcl.pred <- inner_join(cbcl.meds.deltas %>% 
                          filter(drug == "methylphenidate") %>% 
                          pivot_wider(names_from = question, 
                                      values_from = delta, 
                                      id_cols = c(IID, sex)),
                        abcd.pred)
corr.table(cbcl.pred %>% select(predicted),
           cbcl.pred %>% select(starts_with("syn"), starts_with("dsm5")),
           method = "spearman") %>%
  filter(V1 == "predicted", V1 != V2) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", V2), 
                                    "raw data", 
                                    "corrected for age, sex, and interaction"), 
                             levels = c("raw data", 
                                        "corrected for age, sex, and interaction"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  filter(grepl("as", V2)) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.1, paste0(rho, ": ", round(r, 3),
                                                                    ", p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), scales = "free_y") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", nrow(cbcl.pred), "\n",
                        "\tmethylphenidate"), 
       title = "correlation of cbcl score change (delta) of MPH with predicted MPH response") + 
  theme(plot.title = element_text(hjust = 1, size = 9),axis.text.x.bottom = element_text(angle = 0))
####


####
# supplementary figures ---------------------------------------------------
# scatterplot for CBCL scores on and off MPH
cbcl.meds.deltas %>%
  filter(grepl("raw", question), drug == "methylphenidate", !grepl("totprob", question)) %>%
  ggplot(aes(x=`0`, y=`1`)) +
  geom_point(size=0.6)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap("question", scales = "free")+
  labs(x="off-MPH", y="on-MPH", 
       caption = paste0("n(samples): ", nrow(cbcl.meds.deltas %>% filter(drug == "methylphenidate") %>% distinct(IID)), "\n",
                        "the red dashed line represents a line with a slope of 1 and intercepts with the origin (0,0)"))
# histogram for PGS distribution
inner_join(cbcl.meds.deltas %>% 
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
       title = "distribution of PGS for CBCL samples", 
       caption = paste0("n(samples): ", nrow(cbcl.meds.deltas %>% filter(drug == "methylphenidate") %>% distinct(IID))))
# histogram for predicted response distribution for CBCL samples
cbcl.pred %>% 
  ggplot(aes(x=predicted)) +
  geom_histogram()+
  labs(title = "distribution of predicted MPH response for CBCL samples", 
       caption = paste0("n(samples): ", nrow(cbcl.pred  %>% distinct(IID))))
# deltas boxplot for distribution
cbcl.meds.deltas %>% 
  filter(drug == "methylphenidate", grepl("raw", question), !grepl("totprob", question)) %>%
  ggplot(aes(x=question, y=delta, fill = question))+
  geom_boxplot(outlier.size = 0.3, show.legend = F) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  labs(y="score delta")+
  labs(title = "distribution of CBCL deltas for MPH samples", x="",
       caption = paste0("n(samples): ", nrow(cbcl.meds.deltas %>% filter(drug=="methylphenidate")%>% distinct(IID))))
####
# Extras ------------------------------------------------------------------
as.cap <- paste0("delta per question = score_on_the_drug - score_off_the_drug", "\n", 
                 "corrected for:", "\n",
                 "\tinterview_age + sex + interview_age:sex")
####
