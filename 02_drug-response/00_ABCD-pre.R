################################################################################
#                                pre figures for ABCD                          #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(heatmap3)
library(ggpubr);library(ggExtra);library(ggh4x)
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response"
setwd(project.dir)
################################################################################
################################################################################
################################################################################
abcd.raw.dir <- "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core"
abcd.deriv.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/abcd5"
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv")
colnames(abcd.pgs) <- sub("corrected_", "", colnames(abcd.pgs))
demo <- read_csv(paste0(abcd.deriv.dir, "/age-sex-by-eventname.csv"))
abcd.meds <- read_rds(paste0(abcd.deriv.dir, "/../meds/abcd5/abcd5-meds-matrix.rds"))
# abcd.meds <- inner_join(abcd.meds, abcd.pgs) %>%
#   mutate_at(.vars = vars(colnames(abcd.pgs)[-1]), .funs = function(x) scale(x, scale = T, center = T))
tmp <- abcd.meds %>%
  pivot_longer(cols = -c("IID", "eventname"), names_to = "med", values_to = "val") %>%
  filter(val != 0) %>%
  group_by(med) %>%
  summarise(count = n()) %>%
  filter(count>50) %>%
  filter(! med %in% c("sodium", "fluoride", "vitamin", "1151133 pill"))
corr <- cor(abcd.meds%>%select(tmp$med), method = "spearman")
heatmap3(corr%>%as.data.frame()%>%mutate_all(.funs = function(x) ifelse(x==1, NA, x)), 
         scale = "none", symm = T, na.rm = T, col = colorRampPalette(c("white", redblu.col[1]))(100))
p<- heatmap3(corr%>%as.data.frame()%>%mutate_all(.funs = function(x) ifelse(x==1, NA, x)), 
             scale = "none", symm = T, na.rm = T, col = colorRampPalette(c("white", redblu.col[1]))(100), 
             TreturnDistMatrix = T)
hc <- p$hcr
clusters <- data.frame(cluster = hc$order, name = hc$labels[hc$order])
################################################################################
# culture and environemnt
ce <- read_rds(paste0(abcd.deriv.dir, "/culture-environment_c.rds"))
ce.c <- inner_join(ce%>%rename(IID=src_subject_id), abcd.meds) %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(IID, eventname, 
         tmp$med, 
         wps_ss_sum, fes_p_ss_exp_sum, fes_p_ss_cohesion_sum, fes_p_ss_org_sum, 
         fes_p_ss_act_rec_sum, fes_p_ss_int_cult_sum, fes_y_ss_fc,
         pmq_y_ss_mean,crpbi_y_ss_parent, crpbi_y_ss_caregiver,
         pnh_ss_protective_scale,pbp_ss_rule_break, pbp_ss_prosocial_peers,
         sag_days_skip_school, sag_grades_last_yr,
         srpf_y_ss_ses, srpf_y_ss_iiss, psb_p_ss_mean, psb_y_ss_mean
         ) %>%
  mutate_all(.funs = function(x) ifelse(x %in% c(999,666,777), NA, x))
left_join(corr.table(ce.c %>% select(tmp$med), ce.c%>%select(-c(IID, eventname, tmp$med))) %>%
  filter(V1 %in% tmp$med, !V2 %in% tmp$med, V1 != "1151133 pill"),
  data.frame(V2 = c("wps_ss_sum", "fes_p_ss_exp_sum", "fes_p_ss_cohesion_sum", "fes_p_ss_org_sum", 
                    "fes_p_ss_act_rec_sum", "fes_p_ss_int_cult_sum", "fes_y_ss_fc",
                    "parental_monitor_ss_mean", "pmq_y_ss_mean",
                    "mnbs_ss_monitor_supervision", "mnbs_ss_ed_support",
                    "crpbi_y_ss_parent", "crpbi_y_ss_caregiver","peerinfluence_ss_mean",
                    "pnh_ss_protective_scale","pbp_ss_rule_break", "pbp_ss_prosocial_peers",
                    "sag_days_skip_school", "sag_grades_last_yr",
                    "srpf_y_ss_ses", "srpf_y_ss_iiss", "psb_p_ss_mean","psb_y_ss_mean"),
             exp = c("Wills Problem Solving Sum","ExpressIon subscale From the famIly EnvIronment Scale Sum oF Parent Report",
                     "CohesIon subscale from the FamIly EnvIronment Scale Sum oF Parent Report",
                     "OrganIzatIon subscale From the famIly EnvIronment Scale Sum oF Parent Report",
                     "ActIvIty-RecreatIonal subscale from the FamIly EnvIronment Scale Sum oF Parent Report",
                     "Intellectual-Cultural subscale from the FamIly EnvIronment Scale Sum oF Parent Report",
                     "Conflict Subscale from the Family Environment Scale Sum of Youth Report",
                     "Parental Monitoring Summary Score","Parental Monitoring: Mean","Monitoring/Supervision: mean", "Educational/support: mean",
                     "CRPBI - Acceptance Subscale Mean of Report by Parent Completing Protocol by youth",
                     "CRPBI - Acceptance Subscale Mean of Report by Secondary Caregiver by youth",
                     "Peer Influence Summary Score","Peer Network Health: Protective Scale Score",
                     "Involvement with Rule Breaking/Delinquent Peers sum","Involvement with Prosocial Peers sum",
                     "How many days did you skip school without an excuse during the past 4 weeks",
                     "What grades did you receive in school last year?","SRPF School Environment Subscale, Sum", "SRPF School Involvement Subscale, Sum",
                     "Prosocial Behavior Subscale Mean of Parent Report on Youth",
                     "Prosocial Behavior Subscale Mean of Youth Self Report"))) %>%
  drop_na(r) %>%
  ggplot(aes(x=factor(V1,levels = clusters$name), y=V2, fill=r, label=ifelse(pval<0.05, ifelse(pval<0.001, "**","*"), "")))+
  geom_tile()+
  geom_text(size=2)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid(rows = vars(reorder(exp, desc(V2))), space = "free", scales = "free")+
  theme(strip.text.y.right = element_text(angle = 0))+
  my.guides+labs(x="", y="", title =  "correlation between taking meds and culture & environment data",
                 caption = paste0("n(samples): ", nrow(ce.c), "\n",
                                  "*: p-value < 0.05", "\n",
                                  "**: p-value < 0.001"))
################################################################################
# gender identity sexual health
gish <- read_rds(paste0(abcd.deriv.dir, "/gender-identity-sexual-health_c.rds"))
gish.c <- inner_join(gish%>%rename(IID=src_subject_id), abcd.meds) %>%
  # filter(eventname == "baseline_year_1_arm_1") %>%
  select(IID, eventname, 
         tmp$med, 
         gish_p_ss_m_exp_avg, gish_p_ss_f_exp_avg, gish_p_ss_m_dys_avg, 
         gish_p_ss_f_dys_avg, gish_p_ss_m_avg, gish_p_ss_f_avg,
         gish_y_ss_m_avg, gish_y_ss_f_avg
         )%>%
  mutate_all(.funs = function(x) ifelse(x %in% c(999,666,777), NA, x))
left_join(corr.table(gish.c %>% select(tmp$med), gish.c%>%select(-c(IID, eventname, tmp$med))) %>%
  filter(V1 %in% tmp$med, !V2 %in% tmp$med, V1 != "1151133 pill"),
  data.frame(V2 = c("gish_p_ss_m_exp_avg", "gish_p_ss_f_exp_avg", "gish_p_ss_m_dys_avg", "gish_p_ss_f_dys_avg",
                  "gish_p_ss_m_avg", "gish_p_ss_f_avg","gish_y_ss_m_avg", "gish_y_ss_f_avg"),
           exp = c("GISH Male Expression Average", "GISH Female Expression Average",
                   "GISH Male Dysphoria Average", "GISH Female Dysphoria Average",
                   "GISH Male GIQ Average", "GISH Female GIQ Average",
                   "GISH Male Gender Averaged", "GISH Female Gender Averaged"))) %>%
  ggplot(aes(x=factor(V1,levels = clusters$name), y=V2, fill=r, label=ifelse(pval<0.05, ifelse(pval<0.001, "**","*"), "")))+
  geom_tile()+
  geom_text(size=2)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid(rows = vars(reorder(exp, desc(V2))), space = "free", scales = "free")+
  theme(strip.text.y.right = element_text(angle = 0))+
  my.guides+labs(x="", y="", title =  "correlation between taking meds and gender identity & sexual health data",
                 caption = paste0("n(samples): ", nrow(gish.c), "\n",
                                  "*: p-value < 0.05", "\n",
                                  "**: p-value < 0.001"))
################################################################################
# latent factors
lf <- read_rds(paste0(abcd.deriv.dir, "/abcd-general_c.rds"))
lf.c <- inner_join(lf%>%rename(IID=src_subject_id), abcd.meds) %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(IID, eventname, 
         tmp$med, 
         starts_with("latent_factor")
  )%>%
  drop_na() %>%
  mutate_all(.funs = function(x) ifelse(x %in% c(999,666,777), NA, x))
corr.table(lf.c %>% select(tmp$med), lf.c%>%select(-c(IID, eventname, tmp$med))) %>%
  filter(V1 %in% tmp$med, !V2 %in% tmp$med, V1 != "1151133 pill") %>%
  mutate(exp = ifelse(V2 == "latent_factor_ss_general_ses", "economic, social, and physiological well-being",
                      ifelse(V2 == "latent_factor_ss_social", "youth perceived social support", 
                             "perinatal health"))) %>%
  ggplot(aes(x=factor(V1,levels = clusters$name), y=V2, fill=r, label=ifelse(pval<0.05, ifelse(pval<0.001, "**","*"), "")))+
  geom_tile()+
  geom_text(size=2)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid(rows = vars(reorder(exp, desc(V2))), space = "free", scales = "free", )+
  theme(strip.text.y.right = element_text(angle = 0))+
  my.guides+labs(x="", y="", title =  "correlation between taking meds and latent factors data",
                 caption = paste0("n(samples): ", nrow(lf.c), "\n",
                                  "*: p-value < 0.05", "\n",
                                  "**: p-value < 0.001"))
################################################################################
# mental health
mh <- read_rds(paste0(abcd.deriv.dir, "/mental-health_c.rds"))
mh.q <- data.frame(V2 = c("ssrs_p_ss_sum", 
                          "ple_p_ss_affected_good_sum","ple_p_ss_affected_bad_sum",
                          "ple_y_ss_affected_bad_sum","ple_y_ss_affected_good_sum",
                          "abcl_scr_prob_total_r","abcl_scr_prob_external_r",
                          "abcl_scr_prob_internal_r","abcl_scr_prob_other_r",
                          "abcl_scr_prob_intrusive_r","abcl_scr_prob_rulebreak_r",
                          "abcl_scr_prob_aggressive_r","abcl_scr_prob_attention_r",
                          "abcl_scr_prob_thought_r","abcl_scr_prob_somatic_r",
                          "abcl_scr_prob_withdrawn_r","abcl_scr_prob_anxious_r",
                          "abcl_scr_adapt_friends_r","abcl_scr_sub_use_drugs_r",
                          "abcl_scr_sub_use_alcohol_r","abcl_scr_sub_use_tobacco_r",
                          "asr_scr_perstr_r", "asr_scr_anxdep_r","asr_scr_withdrawn_r",
                          "asr_scr_somatic_r", "asr_scr_thought_r", "asr_scr_attention_r",
                          "asr_scr_aggressive_r", "asr_scr_rulebreak_r", "asr_scr_intrusive_r",
                          "asr_scr_internal_r","asr_scr_external_r","asr_scr_totprob_r",
                          "asr_scr_depress_r","asr_scr_anxdisord_r","asr_scr_somaticpr_r",
                          "asr_scr_avoidant_r", "asr_scr_adhd_r", "asr_scr_antisocial_r",
                          "asr_scr_inattention_r","asr_scr_hyperactive_r",
                          "cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                          "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                          "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r","cbcl_scr_syn_internal_r",
                          "cbcl_scr_syn_external_r","cbcl_scr_syn_totprob_r","cbcl_scr_dsm5_depress_r",
                          "cbcl_scr_dsm5_anxdisord_r","cbcl_scr_dsm5_somaticpr_r","cbcl_scr_dsm5_adhd_r",
                          "cbcl_scr_dsm5_opposit_r","cbcl_scr_dsm5_conduct_r","cbcl_scr_07_sct_r",
                          "cbcl_scr_07_ocd_r","cbcl_scr_07_stress_r",
                          "bpm_t_scr_attention_r","bpm_t_scr_internal_r","bpm_t_scr_external_r","bpm_t_scr_totalprob_r",
                          "bpm_y_scr_attention_r","bpm_y_scr_internal_r","bpm_y_scr_external_r","bpm_y_scr_totalprob_r"
                          ),
                   exp = c("Short Social Responsiveness Scale sum", 
                           "How Much Affected Good Sum - life events -p","How Much Affected Bad Sum - life events -p",
                           "How Much Affected Bad Sum - life events -y","How Much Affected Good Sum - life events -y",
                           "Total Problems-ABCL","Externalizing-ABCL","Internalizing-ABCL","Other Problems-ABCL",
                           "Intrusive Syndrome Scale-ABCL","Rule-Breaking Behavior Syndrome Scale-ABCL",
                           "Aggressive Behavior Syndrome Scale-ABCL","Attention Problems Syndrome Scale-ABCL",
                           "Thought Problems Syndrome Scale-ABCL","Somatic Complaints Syndrome Scale-ABCL",
                           "Withdrawn/Depressed Syndrome Scale-ABCL","Anxious/Depressed Syndrome Scale-ABCL",
                           "Adaptive Scale: Friends-ABCL","Drugs-ABCL","Alcohol-ABCL","Tobacco-ABCL",
                           "Personal Strength ASR Adaptive Functioning Scale","Anxious/Depressed ASR Syndrome Scale",
                           "Withdrawn ASR Syndrome Scale", "Somatic Complaints ASR Syndrome Scale",
                           "Thought Problems ASR Syndrome Scale","Attention Problems ASR Syndrome Scale",
                           "Aggressive Behavior ASR Syndrome Scale","Rule-Breaking Behavior ASR Syndrome Scale",
                           "Intrusive ASR Syndrome Scale","Interalizing Problems ASR Syndrome Scale",
                           "Externalizing Problems ASR Syndrome Scale","Total Problems ASR Syndrome Scale",
                           "Depressive Problems ASR DSM-5-Oriented Scale","Anxiety Problems ASR DSM-5-Oriented Scale",
                           "Somatic Problems ASR DSM-5-Oriented Scale","Avoidant Personality Problems ASR DSM-5-Oriented Scale",
                           "AD/H Problems ASR DSM-5-Oriented Scale","Antisocial Personality Problems ASR DSM-5-Oriented Scale",
                           "Inattention ASR DSM-5-Oriented Scale","Hyperactivity-Impulsivity ASR DSM-5-Oriented Scale",
                           "AnxDep CBCL Syndrome Scale", "WithDep CBCL Syndrome Scale","Somatic CBCL Syndrome Scale","Social CBCL Syndrome Scale","Thought CBCL Syndrome Scale", "Attention CBCL Syndrome Scale", 
                           "RuleBreak CBCL Syndrome Scale", "Aggressive CBCL Syndrome Scale", "Internal CBCL Syndrome Scale ", "External CBCL Syndrome Scale", "TotProb CBCL Syndrome Scale","Depress CBCL DSM5 Scale",
                           "AnxDisord CBCL DSM5 Scale","SomaticPr CBCL DSM5 Scale","ADHD CBCL DSM5 Scale","Opposit CBCL DSM5 Scale","Conduct CBCL DSM5 Scale","Sluggish Cognitive Tempo (SCT) CBCL Scale2007 Scale",
                           "Obsessive-Compulsive Problems (OCD) CBCL Scale2007 Scale", "Stress CBCL Scale2007 Scale",
                           "ATT - teacher", "internal - teacher", "external - teacher", "total prob - teacher",
                           "ATT - y", "internal - y", "external - y", "total prob - y"
                           ))
mh.c <- inner_join(mh%>%rename(IID=src_subject_id), abcd.meds) %>%
  # filter(eventname == "baseline_year_1_arm_1") %>%
  select(IID, eventname, 
         tmp$med, 
         mh.q$V2
         )%>%
  # drop_na() %>%
  mutate_all(.funs = function(x) ifelse(x %in% c(999,666,777), NA, x))
left_join(corr.table(mh.c %>% select(tmp$med), mh.c%>%select(-c(IID, eventname, tmp$med))) %>%
            filter(V1 %in% tmp$med, !V2 %in% tmp$med, V1 != "1151133 pill"),
          mh.q) %>%
  ggplot(aes(x=factor(V1,levels = clusters$name), y=V2, fill=r, label=ifelse(pval<0.05, ifelse(pval<0.001, "**","*"), "")))+
  geom_tile()+
  geom_text(size=2)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  facet_grid(rows = vars(reorder(exp, desc(V2))), space = "free", scales = "free", )+
  theme(strip.text.y.right = element_text(angle = 0))+
  my.guides+labs(x="", y="", title =  "correlation between taking meds and latent factors data",
                 caption = paste0("n(samples): ", nrow(mh.c), "\n",
                                  "*: p-value < 0.05", "\n",
                                  "**: p-value < 0.001"))
################################################################################
# wanted to check correlation between PGS and taking medications
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))

pgs.c <- inner_join(abcd.pgs, abcd.meds) %>%
  # filter(eventname == "baseline_year_1_arm_1") %>%
  select(IID, eventname, 
         tmp$med, 
         colnames(abcd.pgs)[-1]
  )%>%
  # drop_na() %>%
  mutate_all(.funs = function(x) ifelse(x %in% c(999,666,777), NA, x))
corr.table(pgs.c %>% select(tmp$med), 
           pgs.c%>%select(-c(IID, eventname, tmp$med))) %>%
  filter(V1 %in% tmp$med, !V2 %in% tmp$med, V1 != "1151133 pill") %>%
  ggplot(aes(x=factor(V1,levels = clusters$name), y=V2, fill=r, label=ifelse(pval<0.05, ifelse(pval<0.001, "**","*"), "")))+
  geom_tile()+
  geom_text(size=2)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], name = "ρ")+
  # facet_grid(rows = vars(V2), space = "free", scales = "free", )+
  theme(strip.text.y.right = element_text(angle = 0))+
  my.guides+labs(x="", y="", title =  "correlation between taking meds and PGS",
                 caption = paste0("n(samples): ", nrow(pgs.c), "\n",
                                  "*: p-value < 0.05", "\n",
                                  "**: p-value < 0.001"))

################################################################################
# correlation between PGS and cbcl
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))
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

abcd.cbcl.filt <- inner_join(abcd.cbcl.filt, demo) %>% select(IID, eventname, interview_age,sex, 
                                                                     starts_with("syn"), starts_with("dsm5")) %>%
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
# combine raw and corrected sst data
cbcl.all <- inner_join(abcd.cbcl.filt,
                           cbcl.as.corrected %>% rename_at(.vars = vars(starts_with("syn"), starts_with("dsm5")), 
                                                                          .funs = function(x) ifelse(grepl("syn", x), sub("syn_raw_", "syn_as_", x),
                                                                                                     sub("dsm5_raw_", "dsm5_as_", x))))
pgs.cbcl <- left_join(cbcl.all, abcd.pgs)
corr.table(pgs.cbcl %>% select(starts_with("syn"), starts_with("dsm")),
           pgs.cbcl %>% select(colnames(abcd.pgs), -IID),
           method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs), !V2 %in% colnames(abcd.pgs)) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                                "corrected for age, sex, and interaction", 
                                                                                "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  mutate(subscale = ifelse(grepl("dsm5", V2), "dsm5", "syn")) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3) +
  redblu.col.gradient+my.guides+null_labs +
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  labs(caption = paste0("n(samples): ", nrow(pgs.cbcl), "\n",
                        "* pval < 0.05", "\n", 
                        "** pval < 0.01", "\n", 
                        "*** FDR < 0.05"))
#########################
# make same figure, but by item level and not syn
tt <- read_rds("/wdata/msmuhammad/data/ABCD/cbcl-scales.rds")
abcd.cbcl.filt.item <- inner_join(abcd.cbcl, demo) %>% select(IID, eventname, interview_age,sex, 
                                                              starts_with("cbcl_q")) %>%
  drop_na() 
################################################################################
# correlation between PGS and participants performance in SST run 1
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))
demo <- read_csv(paste0(abcd.deriv.dir, "/age-sex-by-eventname.csv"))
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
                    demo) %>% drop_na() %>%
  mutate_at(.vars = vars(3:7, interview_age), .funs = function(x) as.numeric(x))
# correct for age, sex, their interaction
sst.r1.as.corrected <- cbind(sst.r1 %>% select(IID, interview_age, sex, eventname), 
                                 apply(sst.r1 %>% 
                                         select(starts_with("e_")), 
                                       MARGIN = 2, FUN = function(x) {
                                         residuals(glm(y ~ interview_age + sex + interview_age:sex, 
                                                       data = sst.r1 %>% mutate(y = x), 
                                                       family = poisson()))
                                       }))
sst.r1.all <- inner_join(sst.r1,
                             sst.r1.as.corrected %>% rename_at(.vars = vars(starts_with("e_")), 
                                                               .funs = function(x) sub("e_raw_", "e_as_", x)))
pgs.sst <- inner_join(abcd.pgs, sst.r1.all)
corr.table(pgs.sst %>% select(starts_with("e_")),
           pgs.sst %>% select(colnames(abcd.pgs), -IID),
           method = "spearman") %>%
  filter(V1 %in% colnames(abcd.pgs), !V2 %in% colnames(abcd.pgs)) %>%
  mutate(value_type = factor(ifelse(grepl("raw_", V2), "raw data", ifelse(grepl("as_", V2), 
                                                                          "corrected for age, sex, and interaction", 
                                                                          "corrected for age, sex, interaction, and other ADHD meds")), 
                             levels = c("raw data", "corrected for age, sex, and interaction", 
                                        "corrected for age, sex, interaction, and other ADHD meds"))) %>%
  group_by(value_type) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3) +
  redblu.col.gradient+my.guides+null_labs +
  facet_grid2(rows = vars(value_type), scales = "free_y", independent = "y") +
  labs(caption = paste0("n(samples): ", nrow(pgs.sst), "\n",
                        "* pval < 0.05", "\n", 
                        "** pval < 0.01", "\n", 
                        "*** FDR < 0.05"))
################################################################################

################################################################################

################################################################################
