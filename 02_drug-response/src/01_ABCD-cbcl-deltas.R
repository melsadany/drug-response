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
# CBCL deltas per drug with PGS -------------------------------------------
# scatterplot
t <- inner_join(cbcl.meds.deltas %>% 
                  filter(drug == "methylphenidate") ,
                abcd.pgs) %>%
  # filter(grepl("attention", question) | grepl("adhd", question)) %>%
  filter(grepl("attention", question) | grepl("adhd", question) | grepl("social", question) | grepl("conduct", question)) %>%
  filter(grepl("as", question)) %>%  mutate(question = sub("as_", "", question)) %>%
  pivot_longer(cols = c("ADHD-Demontis", "cog_gFactor"), names_to = "PGS", values_to = "score") %>%
  select(IID, sex, question, delta, n_samples, PGS, score, drug)
t %>%
  ggplot(aes(y=delta, x=score))+
  geom_point(size=0.3) +
  facet_grid2(cols = vars(PGS), rows = vars(question), scales = "free") +
  stat_cor(method = "spearman") + geom_smooth(method = "glm") +
  labs(y="CBCL delta (on-drug - off-drug)", x = "polygenic score",
       title = "correlation between PGS and CBCL deltas",
       caption = paste0("n(samples): ", length(unique(t$IID)), "\n",
                        # unique(t$drug)
                        # "\tguanfacine = guanfacine / tenex"
                        "\tmethylphenidate = methylphenidate / ritalin"
                        # "non_stim =intuniv/strattera/tenex/atomoxetine/clonidine/guanfacine"
                        # "stim = methylphenidate/ritalin/concerta/amphetamine/adderall/vyvanse/lisdexamfetamine"
                        ))
# heatmaps
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
  # filter(drug == "methylphenidate") %>%
  filter(drug %in% c("stim", "non_stim")) %>%
  # filter(drug == "guanfacine") %>%
  filter(grepl("as", V2)) %>% mutate(V2 = sub("as_", "", V2)) %>%
  filter(grepl("ADHD", V1) | grepl("cog_gFa", V1)) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3)+
  facet_grid2(rows = vars(value_type), cols = vars(drug), scales = "free_y") +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", "\n\t",
                        # paste(apply(cbcl.meds.deltas.pgs%>%ungroup()%>%select(drug, n_samples)%>%distinct(), 
                        # "methylphenidate: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "methylphenidate/ritalin/concerta: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "methylphenidate/ritalin: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="methylphenidate")%>%distinct(n_samples), "\n",
                        # "concerta: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="concerta")%>%distinct(n_samples), "\n",
                        # "guanfacine/tenex/intuniv: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="guanfacine")%>%distinct(n_samples), "\n",
                        # "guanfacine/tenex: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="guanfacine")%>%distinct(n_samples), "\n",
                        # "lisdexamfetamine/vyvanse: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="lisdexamfetamine")%>%distinct(n_samples), "\n",
                        # "amphetamine/adderall: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="amphetamine")%>%distinct(n_samples), "\n",
                        # "clonidine: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="clonidine")%>%distinct(n_samples), "\n",
                        # "atomoxetine/strattera: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="atomoxetine")%>%distinct(n_samples), "\n",
                        "stim = methylphenidate/ritalin/concerta/amphetamine/adderall/vyvanse/lisdexamfetamine: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="stim")%>%distinct(n_samples), "\n",
                        "\tnon_stim = intuniv/strattera/tenex/atomoxetine/clonidine/guanfacine: ", cbcl.meds.deltas.pgs %>%ungroup()%>% filter(drug=="non_stim")%>%distinct(n_samples), "\n",
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
cbcl.pred %>%
  ggplot(aes(x=predicted, y=dsm5_as_adhd))+
  geom_point()+
  geom_smooth(method = "glm")+stat_cor(method = "spearman") +
  labs(y = "delta of DSM5_ADHD_symptoms", x = "predicted MPH response",
       title = "correlation of cbcl adhd change (delta) of MPH with predicted MPH response",
       caption = paste0("n(samples): ", length(unique(cbcl.pred$IID)), "\n",
                        "\tmethylphenidate = methylphenidate / ritalin"))
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
  # filter(grepl("as", V2)) %>%
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
    mm <- m232
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
drug.deltas.rsquared %>%
  pivot_longer(cols = starts_with("m"), names_to = "model", values_to = "rsquared") %>%
  mutate(model = factor(model, levels = c("m1", "m2", "m3", "m12", "m13", "m23", "m232", "m123"))) %>%
  mutate(question = sub("_as", "", question)) %>%
  # filter(!(grepl("3", model)&drug %in% c("amphetamine", "atomoxetine", "clonidine", "guanfacine", "lisdexamfetamine", "non_stim"))) %>%
  ggplot(aes(x=model, y=rsquared, fill=model))+
  geom_bar(stat = "identity", show.legend = F)+
  facet_grid2(rows = vars(question), cols = vars(drug), scales = "free_y", independent = "y")+
  theme(strip.text.y.right = element_text(angle = 0))+
  scale_fill_manual(values = c(six.colors, boxplot.colors)) +
  labs(title = "predicting CBCL deltas by PGS or predicted MPH response",
       caption = paste0("all models are for predicting the delta using different variables\n",
                        "\tm1\t\tADHD_PGS\n",
                        "\tm2\t\tgFactor_PGS\n",
                        "\tm3\t\tpredcited MPH response\n",
                        "\tm12\t\tADHD_PGS + gFactor\n",
                        "\tm13\t\tADHD_PGS + predicted\n",
                        "\tm23\t\tgFactor_PGS + predicted\n",
                        "\tm232\tgFactor_PGS + predicted + gFactor_PGS:predicted\n",
                        "\tm123\tADHD_PGS + gFactor_PGS + predicted\n"))
drug.deltas.rsquared %>%
  mutate(question = sub("\\.[0-9]", "", question)) %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  mutate(question = sub("as_", "", question)) %>%
  mutate(var = factor(var, levels = c("ADHD-Demontis", "cog_gFactor", "predicted", "cog_gFactor:predicted"))) %>%
  ggplot(aes(x=Estimate, y = var)) +
  geom_point(aes(alpha = sig, color = var),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.5)) +
  # scale_shape_manual(values = c(1, 2)) + 
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig, color = var), 
                 linewidth = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = six.colors[c(2,3,5)]) +
  facet_grid2(rows = vars(question), cols = vars(drug))+
  theme(strip.text.y.right = element_text(angle = 0)) +
  labs(x = "Estimate", y="", 
       # title = "predicting CBCL deltas by PGS or predicted MPH response",
       # title = "predicting CBCL deltas by PGS",
       title = "predicting CBCL deltas",
       # title = "predicting CBCL deltas by predicted MPH response",
       caption = paste0("n(samples): \n", 
                        "\tmethylphenidate/ritalin: ", drug.deltas.rsquared%>% filter(drug=="methylphenidate")%>%distinct(n), "\n",
                        "\tconcerta: ", drug.deltas.rsquared%>% filter(drug=="concerta")%>%distinct(n), "\n",
                        "\tguanfacine/tenex: ", drug.deltas.rsquared%>% filter(drug=="guanfacine")%>%distinct(n), "\n",
                        "\tlisdexamfetamine/vyvanse: ", drug.deltas.rsquared%>% filter(drug=="lisdexamfetamine")%>%distinct(n), "\n",
                        "\tamphetamine/adderall: ", drug.deltas.rsquared%>% filter(drug=="amphetamine")%>%distinct(n), "\n",
                        "\tclonidine: ", drug.deltas.rsquared%>% filter(drug=="clonidine")%>%distinct(n), "\n",
                        "\tatomoxetine/strattera: ", drug.deltas.rsquared%>% filter(drug=="atomoxetine")%>%distinct(n), "\n",
                        "\tstim = methylphenidate/ritalin/concerta/amphetamine/adderall/vyvanse/lisdexamfetamine: ", drug.deltas.rsquared%>% filter(drug=="stim")%>%distinct(n), "\n",
                        "\tnon_stim = intuniv/strattera/tenex/atomoxetine/clonidine/guanfacine: ", drug.deltas.rsquared%>% filter(drug=="non_stim")%>%distinct(n)
                        ))
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
# correlation between predicted MPH response and PGS
pgs.pred <- inner_join(abcd.pred, abcd.pgs)
# pgs.pred <- inner_join(inner_join(abcd.pred, abcd.pgs), cbcl.meds.deltas %>% select(IID, drug)%>%filter(drug=="methylphenidate"))%>%distinct(IID ,.keep_all = T)
pgs.pred %>% 
  pivot_longer(cols = colnames(abcd.pgs)[-1], names_to = "PGS", values_to = "score") %>%
  filter(grepl("gF", PGS)) %>%
  ggplot(aes(x=predicted, y=score))+
  geom_point(size=0.3)+
  geom_smooth(method = "glm") +
  stat_cor() +
  # facet_wrap("PGS", scales = "free_y") +
  labs(y="cognitive gFactor PGS", x="predicted MPH response",
       caption = paste0("n(samples): ", nrow(pgs.pred)
                        ,"\n\tmethylphenidate = methylphenidate / ritalin"),
       title = "correlation between predicted MPH response and PGS")
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
# Extras ------------------------------------------------------------------
as.cap <- paste0("delta per question = score_on_the_drug - score_off_the_drug", "\n", 
                 "corrected for:", "\n",
                 "\tinterview_age + sex + interview_age:sex")
####
