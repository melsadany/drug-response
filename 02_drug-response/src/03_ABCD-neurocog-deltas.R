################################################################################
#                 change in other neurocog tasks based on drug status          #
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
# ABCD predicted MPH response ---------------------------------------------
abcd.pred <- read_rds("../data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(predicted = m) %>%
  # filter(IID %in% sst.r1.deltas$IID) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
####
# nihtbx task ---------------------------------------------------------------
abcd.nihtbx <- read_csv(paste0(abcd.raw.dir, "/neurocognition/nc_y_nihtb.csv"))
# keep variables of interest 
# only keeping the anxiety/frustration/irritability/happiness statuses pre and post task
abcd.nihtbx.filt <- abcd.nihtbx %>%
  select(IID = src_subject_id, eventname, 
         contains("agecorrected")) %>%
  filter(rowSums(is.na(abcd.nihtbx%>%select(contains("agecorrected"))))<=5)
  # drop_na()
# combine nihtbx, age, sex, meds data 
abcd.nihtbx.filt <- inner_join(inner_join(demo, abcd.nihtbx.filt), abcd.meds) %>% 
  # keep this order for OCD pref
  select(IID, eventname, interview_age,sex, 
         starts_with("nih"), 
         colnames(abcd.meds)) %>%
  drop_na()
####
# nihtbx performance based on PGS and categorized by drug -----------------
pgs.nihtbx.drug <- inner_join(abcd.nihtbx.filt, abcd.pgs)
pgs.nihtbx.drug %>% 
  pivot_longer(cols = starts_with("nih"), names_to = "task", values_to = "task_score") %>%
  pivot_longer(cols = colnames(abcd.pgs)[-1], names_to = "pgs", values_to = "pg_score") %>%
  pivot_longer(cols = adhd.meds$drug, names_to = "drug", values_to = "drug_status") %>%
  filter(drug_status == 1) %>% 
  filter(drug %in% c("stim", "non_stim")) %>%
  filter(grepl("ADHD", pgs) | grepl("gFa", pgs)) %>%
  mutate(task =sub("_agecorrected", "", task)) %>% mutate(task =sub("nihtbx_", "", task)) %>%
  ggplot(aes(x=pg_score, y=task_score, color = drug)) +
  # geom_point(size=0.1)+
  geom_smooth(method = "glm", se = F) +
  stat_cor(show.legend = F)+
  # facet_wrap("task", scales = "free")
  facet_grid2(rows = vars(task), cols = vars(pgs), scales = "free_y") +
  labs(x="PGS", caption = paste0("n(samples):\n",
                                 "\tstim: ", sum(pgs.nihtbx.drug$stim), "\n",
                                 "\tnon_stim: ", sum(pgs.nihtbx.drug$non_stim), "\n"))



####
# nihtbx deltas by drug status for participants -----------------------------
registerDoMC(cores = 3)
nihtbx.meds.deltas <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d <- adhd.meds$drug[i]
  # print(i)
  # get number of participants that match this criteria
  # keep participants of 2 events and 2 distinct drug status at least
  n <- nrow(abcd.nihtbx.filt %>%
              rename(drug = which(colnames(abcd.nihtbx.filt)==d)) %>%
              group_by(IID) %>%
              filter(n_distinct(eventname) >= 2 & n_distinct(drug) > 1) %>%
              ungroup() %>%
              distinct(IID))
  if (n>0) {
    # calculate the delta of each nihtbx variable
    # delta defines as (score_on_the_drug - score_off_the_drug)
    t <- abcd.nihtbx.filt %>%
      rename(drug = which(colnames(abcd.nihtbx.filt)==d)) %>%
      group_by(IID) %>%
      filter(n_distinct(eventname) >= 2 & n_distinct(drug) > 1) %>%
      ungroup() %>%
      # some participants have multiple nihtbx datapoints for the same drug status
      # get the average by drug status after being grouped based on IID and drug binary status
      group_by(IID, drug) %>%
      mutate_at(.vars = vars(starts_with("nihtbx")), .funs = function(x) mean(x)) %>%
      ungroup() %>%
      # make sure to keep one data point per distinct IID, drug status
      distinct(IID, drug, .keep_all = T) %>%
      pivot_longer(cols = c(starts_with("nihtbx")), 
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
  }else {
    return(NULL)
  }
}
####
# nihtbx deltas per drug with PGS -------------------------------------------
# heatmaps
nihtbx.meds.deltas.pgs <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d <- adhd.meds$drug[i]
  t <- inner_join(nihtbx.meds.deltas %>% 
                    filter(drug == d) %>% 
                    pivot_wider(names_from = question, 
                                values_from = delta, 
                                id_cols = c(IID, sex)),
                  abcd.pgs)
  if (nrow(t)>4) {
    ret <- corr.table(t %>% select(colnames(abcd.pgs)[-1]),
                      t %>% select(starts_with("nih")),
                      method = "spearman") %>%
      filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
      mutate(value_type = factor(ifelse(grepl("raw_", V2), 
                                        "raw data", 
                                        "corrected for age"), 
                                 levels = c("raw data", 
                                            "corrected for age"))) %>%
      group_by(value_type) %>%
      mutate(FDR = p.adjust(pval, method = "fdr")) %>%
      mutate(drug = d) %>%
      mutate(n_samples = nrow(t))
    return(ret)
  }else {
    return(NULL)
  }
  
}
nihtbx.meds.deltas.pgs %>%
  # filter(drug == "methylphenidate") %>%
  filter(drug %in% c("stim", "non_stim")) %>%
  # filter(drug == "guanfacine") %>%
  # filter(grepl("as", V2)) %>% mutate(V2 = sub("as_", "", V2)) %>%
  filter(grepl("ADHD", V1) | grepl("cog_gFa", V1)) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), cols = vars(drug), scales = "free_y") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", "\n\t",
                        # paste(apply(nihtbx.meds.deltas.pgs%>%ungroup()%>%select(drug, n_samples)%>%distinct(), 
                        # "methylphenidate: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "methylphenidate/ritalin/concerta: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "methylphenidate/ritalin: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "concerta: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="concerta")%>%distinct(n_samples), "\n",
                        # "guanfacine/tenex/intuniv: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="guanfacine")%>%distinct(n_samples), "\n",
                        # "guanfacine/tenex: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="guanfacine")%>%distinct(n_samples), "\n",
                        # "lisdexamfetamine/vyvanse: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="lisdexamfetamine")%>%distinct(n_samples), "\n",
                        # "amphetamine/adderall: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="amphetamine")%>%distinct(n_samples), "\n",
                        # "clonidine: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="clonidine")%>%distinct(n_samples), "\n",
                        # "atomoxetine/strattera: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="atomoxetine")%>%distinct(n_samples), "\n",
                        "stim = methylphenidate/ritalin/concerta/amphetamine/adderall/vyvanse/lisdexamfetamine: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="stim")%>%distinct(n_samples), "\n",
                        "\tnon_stim = intuniv/strattera/tenex/atomoxetine/clonidine/guanfacine: ", nihtbx.meds.deltas.pgs %>%ungroup()%>% filter(drug=="non_stim")%>%distinct(n_samples), "\n",
                        "* pval < 0.05 & not FDR sig", "\n", 
                        "** pval < 0.01 & not FDR sig", "\n", 
                        "*** FDR < 0.05"), 
       # title = "correlation of nihtbx score change (delta) per drug with PGS")
       title = "correlation of nihtbx score change (delta) with PGS")
####

####



####
