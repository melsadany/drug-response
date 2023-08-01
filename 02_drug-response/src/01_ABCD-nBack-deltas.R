################################################################################
#                  check nback change by status of taking ADHD drugs             #
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
  mutate(methylphenidate = ifelse((methylphenidate+ritalin)>=1 & concerta == 0,1,0),
         amphetamine = ifelse((adderall+amphetamine)>=1,1,0),
         lisdexamfetamine = ifelse((vyvanse+lisdexamfetamine)>=1,1,0),
         stim = ifelse((methylphenidate+adderall+amphetamine+concerta+ritalin+vyvanse+lisdexamfetamine)>=1,1,0),
         non_stim = ifelse((intuniv+strattera+tenex+atomoxetine+clonidine+guanfacine)>=1,1,0),
         # guanfacine = ifelse((guanfacine+tenex+intuniv)>=1,1,0),
         guanfacine = ifelse((guanfacine+tenex)>=1 & intuniv == 0,1,0),
         atomoxetine = ifelse((atomoxetine+strattera)>=1,1,0)) %>% 
  select(-c(ritalin, adderall, tenex, strattera,intuniv, vyvanse))
####
# ABCD PGS file -----------------------------------------------------------
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, contains("cog")&contains("gFa"), contains("PGC")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", "")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-PGC-20[0-9]+", ""))
####
# ABCD predicted MPH response ---------------------------------------------
abcd.pred <- read_rds("../data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(predicted = m) %>%
  # filter(IID %in% nback.r1.deltas$IID) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
####
# ABCD nback data -----------------------------------------------------------
nback.raw <- read_csv(paste0(abcd.raw.dir, "/imaging/mri_y_tfmr_nback_beh.csv"))
# only keeping questions of interest for the first run of the task only
nback.r1 <- left_join(nback.raw %>%
                        filter(tfmri_nback_beh_performflag == 1) %>%
                        select(IID = src_subject_id, eventname,
                               tfmri_nb_r1_beh_cplace_nt,
                               tfmri_nb_r1_beh_cpf_nt,
                               tfmri_nb_r1_beh_cngf_nt,
                               tfmri_nb_r1_beh_cnf_nt),
                      left_join(abcd.meds, demo)) %>% drop_na()
# combine nback, age, sex, meds data 
abcd.nback.filt <- inner_join(inner_join(demo, nback.r1), abcd.meds) %>% 
  # keep this order for OCD pref
  select(IID, eventname, interview_age,sex, 
         ends_with("nt"), 
         colnames(abcd.meds)) %>%
  drop_na()
####
# nback correction for age and sex -----------------------------------------
nback.as.corrected <- cbind(abcd.nback.filt %>% 
                            select(IID, interview_age, sex, eventname), 
                          apply(abcd.nback.filt %>% 
                                  select(ends_with("nt")), 
                                MARGIN = 2, FUN = function(x) {
                                  residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                data = abcd.nback.filt %>% mutate(y = x), 
                                                family = poisson()))
                                 }))
# combine nback raw with nback corrected
nback.all <- inner_join(abcd.nback.filt,
                      nback.as.corrected %>% 
                        rename_at(.vars = vars(ends_with("nt")), 
                                  .funs = function(x) sub("tfmri", "tfmri_as", x)))
####
# nback deltas by drug status for participants -----------------------------
registerDoMC(cores = 3)
nback.meds.deltas <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d <- adhd.meds$drug[i]
  # print(i)
  # get number of participants that match this criteria
  # keep participants of 2 events and 2 distinct drug status at least
  n <- nrow(nback.all %>%
              rename(drug = which(colnames(nback.all)==d)) %>%
              group_by(IID) %>%
              filter(n_distinct(eventname) >= 2 & n_distinct(drug) > 1) %>%
              ungroup() %>%
              distinct(IID))
  # calculate the delta of each nback question
  # delta defines as (score_on_the_drug - score_off_the_drug)
  t <- nback.all %>%
    rename(drug = which(colnames(nback.all)==d)) %>%
    group_by(IID) %>%
    filter(n_distinct(eventname) >= 2 & n_distinct(drug) > 1) %>%
    ungroup() %>%
    # some participants have multiple nback datapoints for the same drug status
    # get the average by drug status after being grouped based on IID and drug binary status
    group_by(IID, drug) %>%
    mutate_at(.vars = vars(ends_with("nt")), .funs = function(x) mean(x)) %>%
    ungroup() %>%
    # make sure to keep one data point per distinct IID, drug status
    distinct(IID, drug, .keep_all = T) %>%
    pivot_longer(cols = c(ends_with("nt")), 
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
# nback deltas per drug with PGS -------------------------------------------
# scatterplot
t <- inner_join(nback.meds.deltas %>% 
                  filter(drug == "stim") ,
                abcd.pgs) %>%
  filter(grepl("as", question)) %>% 
  mutate(question = sub("tfmri_as_nb_r1_beh_", "", sub("_nt", "", question))) %>%
  pivot_longer(cols = c("ADHD-Demontis", "cog_gFactor"), names_to = "PGS", values_to = "score") %>%
  select(IID, sex, question, delta, n_samples, PGS, score, drug)
t %>%
  ggplot(aes(x=delta, y=score))+
  geom_point(size=1) +
  facet_grid2(rows = vars(PGS), cols = vars(question), scales = "free") +
  stat_cor(method = "spearman") + geom_smooth(method = "glm") +
  labs(x="nBack delta (on-drug - off-drug)", y = "polygenic score",
       title = "correlation between PGS and nback deltas",
       caption = paste0("n(samples): ", length(unique(t$IID)), "\n",
                        unique(t$drug)
                        # "guanfacine = guanfacine/tenex"
                        # "methylphenidate = methylphenidate/ritalin"
                        # "non_stim =intuniv/strattera/tenex/atomoxetine/clonidine/guanfacine"
                        # "stim = methylphenidate/ritalin/concerta/amphetamine/adderall/vyvanse/lisdexamfetamine"
       ))
# heatmaps
nback.meds.deltas.pgs <- foreach (i = 1:length(adhd.meds$drug), .combine = rbind) %dopar% {
  d <- adhd.meds$drug[i]
  t <- inner_join(nback.meds.deltas %>% 
                    filter(drug == d) %>% 
                    pivot_wider(names_from = question, 
                                values_from = delta, 
                                id_cols = c(IID, sex)),
                  abcd.pgs)
  ret <- corr.table(t %>% select(colnames(abcd.pgs)[-1]),
                    t %>% select(starts_with("tfmri")),
                    method = "spearman") %>%
    filter(V1 %in% colnames(abcd.pgs)[-1], !V2 %in% colnames(abcd.pgs)[-1]) %>%
    mutate(value_type = factor(ifelse(grepl("tfmri_nb", V2), 
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
nback.meds.deltas.pgs %>%
  # filter(drug == "methylphenidate") %>%
  # filter(drug %in% c("stim", "non_stim")) %>%
  # filter(drug == "guanfacine") %>%
  filter(grepl("as", V2)) %>% 
  mutate(V2 = sub("tfmri_as_nb_r1_beh_", "", sub("_nt", "", V2))) %>%
  # filter(grepl("ADHD", V1) | grepl("cog_gFa", V1)) %>%
  ggplot(aes(y=V1, x=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), cols = vars(drug), scales = "free_x", independent = "x") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", "\n",
                        "\tmethylphenidate/ritalin: ", nback.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        "\tconcerta: ", nback.meds.deltas.pgs %>%ungroup()%>% filter(drug=="concerta")%>%distinct(n_samples), "\n",
                        "\tguanfacine/tenex: ", nback.meds.deltas.pgs %>%ungroup()%>% filter(drug=="guanfacine")%>%distinct(n_samples), "\n",
                        "\tclonidine: ", nback.meds.deltas.pgs %>%ungroup()%>% filter(drug=="clonidine")%>%distinct(n_samples), "\n",
                        "\tatomoxetine/strattera: ", nback.meds.deltas.pgs %>%ungroup()%>% filter(drug=="atomoxetine")%>%distinct(n_samples), "\n",
                        "\tlisdexamfetamine/vyvanse: ", nback.meds.deltas.pgs %>%ungroup()%>% filter(drug=="lisdexamfetamine")%>%distinct(n_samples), "\n",
                        "\tamphetamine/adderall: ", nback.meds.deltas.pgs %>%ungroup()%>% filter(drug=="amphetamine")%>%distinct(n_samples), "\n",
                        "\tstim = methylphenidate/ritalin/concerta/amphetamine/adderall: ", nback.meds.deltas.pgs %>%ungroup()%>% filter(drug=="stim")%>%distinct(n_samples), "\n",
                        "\tnon_stim = vyvanse/intuniv/strattera/tenex/lisdexamfetamine/atomoxetine/clonidine/guanfacine: ", nback.meds.deltas.pgs %>%ungroup()%>% filter(drug=="non_stim")%>%distinct(n_samples), "\n",
                        "variables description:\n",
                        "\tcnf: correct neutral face\n",
                        "\tcngf: correct negative face\n",
                        "\tcpf: correct positive face\n",
                        "\tcplace: correct place\n",
                        "* pval < 0.05 & not FDR sig", "\n", 
                        "** pval < 0.01 & not FDR sig", "\n", 
                        "*** FDR < 0.05"), 
       title = "correlation of nback score change (delta) with PGS")
####

####


####
# supplementary figs ------------------------------------------------------

####
# supplementary tables ----------------------------------------------------

####
# Extras ------------------------------------------------------------------

####
