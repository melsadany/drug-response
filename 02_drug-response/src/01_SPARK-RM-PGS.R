################################################################################
#                               SPARK data checks                              #
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
# SPARK research match data -----------------------------------------------
spark.rm <- read_tsv("../data/rawdata/spark-rm/raw/ChildSurvey20181031.csv")
data.dict <- readxl::read_excel("../data/rawdata/spark-rm/data-dict-child_edited.xlsx", sheet = 1)
# keep basic ID questions, age, sex, and questions starting of q95
spark.rm.filt <- spark.rm[,c(1:4, 160:ncol(spark.rm))]
####
# reformat SPARK answers --------------------------------------------------
# some questions needed reformatting for the encoding of their answers
# I changed the y/n questions to have 0 as NO and 1 as YES
# non-text y/n questions to modify their values
yn.flip = list("q95_alias" = T,
               "q95_sq01_alias" = T,
               "q95_sq02_alias" = T,
               "q96_alias" = T,
               "q96_sq02_alias" = T,
               "q96_sq03_alias" = T,
               "q97_alias" = T,
               "q97_sq02_alias" = T,
               "q97_sq03_alias" = T,
               "q98_alias" = T,
               "q98_sq02_alias" = T,
               "q98_sq03_alias" = T,
               "q99_alias" = T,
               "q99_sq02_alias" = T,
               "q99_sq03_alias" = T,
               "q100_alias" = T,
               "q102_alias" = T,
               "q102_sq02_alias" = T,
               "q102_sq03_alias" = T,
               "q103_alias" = T,
               "q103_sq02_alias" = T,
               "q103_sq03_alias" = T)
# make the flips of yn questions so that yes =1, no =0, NA
for (i in names(yn.flip)){
  for (j in 1:nrow(spark.rm.filt)){
    if (spark.rm.filt[j,i] == 2){spark.rm.filt[j,i] = 0}
    else if(spark.rm.filt[j,i] == -666){spark.rm.filt[j,i] = NA} # missing answers
    else if(spark.rm.filt[j,i] == -999){spark.rm.filt[j,i] = NA} # missing answers
  }
}
colnames(spark.rm.filt)[5:ncol(spark.rm.filt)] <- data.dict$human_friendly_name[157:nrow(data.dict)]
####
# correlation between PGS and meds history --------------------------------
spark.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))
pgs.rm <- inner_join(spark.rm.filt %>% select(IID = ParticipantID, ends_with("yn"), ends_with("effect")),
                     spark.pgs) %>%
  filter(!is.na(ch_med_antihist_yn))
samples <- colSums(pgs.rm %>% select(starts_with("ch_med"), -contains("PM")), na.rm = T) %>% as.data.frame() %>% rename(count = 1) %>% rownames_to_column("q")
corr.table(pgs.rm %>% select(colnames(spark.pgs)[-1]), 
           pgs.rm %>% select(ends_with("yn"), ends_with("effect")),
           method = "spearman") %>%
  filter(V1 %in% colnames(spark.pgs)[-1], !V2 %in% colnames(spark.pgs)[-1]) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  mutate(V2 = sub("ch_med_", "", V2)) %>%
  filter(!grepl("PM", V2)) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3) +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", nrow(pgs.rm), "\n\t",
                        paste(apply(samples, 
                                    1, function(x) paste(sub("ch_med_", "", x[1]), " answering yes:", x[2])), collapse = "\n\t"), "\n",
                        "* pval < 0.05 & not FDR sig", "\n", 
                        "** pval < 0.01 & not FDR sig", "\n", 
                        "*** FDR < 0.05"))
# correlation between SPARK methyl effect and PGS
corr.table(pgs.rm %>% select(colnames(spark.pgs)[-1]), 
           pgs.rm %>% select(ch_med_methyl_effect),
           method = "spearman") %>%
  filter(V1 %in% colnames(spark.pgs)[-1], !V2 %in% colnames(spark.pgs)[-1]) %>%
  mutate(V2 = sub("ch_med_", "", V2)) %>%
  ggplot(aes(x=V2, y=V1, fill = r, label = ifelse(pval<0.1, 
                                                  paste0(rho, ": ", round(r, 3), "\n",
                                                         "p: ", round(pval, 5)), ""))) +
  geom_tile()+
  geom_text(size = 3) +
  redblu.col.gradient+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", nrow(pgs.rm%>%drop_na(ch_med_methyl_effect)))) +
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5))
####
# correlation of MPH effect and SCQ ---------------------------------------
# check correlation between mph_effect vector and subject communication?
spark.scq <- read_csv("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12/scq_2022-12-12.csv")
mph.scq <- inner_join(spark.rm.filt%>%select(IID=ParticipantID, 
                                             interview_age=Dependent_age_at_eval_start_months,
                                             methyl_effect = ch_med_methyl_effect)%>%
                        drop_na(methyl_effect), 
                      spark.scq%>%rename(IID=subject_sp_id))
# figure for 
do.call(rbind,
        lapply(mph.scq %>% select(starts_with("q")), 
               function(x) {
                 t <- fisher.test(mph.scq$methyl_effect, x)
                 data.frame(pval = t$p.value,
                            OR = t$estimate[[1]],
                            confint_min = as.data.frame(t$conf.int)[1,1],
                            confint_max = as.data.frame(t$conf.int)[2,1])})) %>%
  rownames_to_column("question") %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  mutate(question = as.factor(question)) %>%
  ggplot(aes(x=OR, y = question)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_shape_manual(values = c(1, 2)) + 
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = six.colors[1:4]) +
  labs(x = "Odds Ratio", y="", caption = paste0("n(samples): ", nrow(mph.scq), "\n", 
                                      "odds ratio from Fisher's Exact test"))
####
# correlation of MPH effect and cog. impairment ---------------------------
spark.iq <- read_csv("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12/predicted_iq_experimental_2022-12-12.csv") %>% 
  select(IID = subject_sp_id, derived_cog_impair) %>%
  drop_na()
spark.cog <- inner_join(spark.rm.filt%>%select(IID=ParticipantID, 
                                               methyl_effect=ch_med_methyl_effect) %>%
                          drop_na(methyl_effect), 
                        spark.iq)
t <- fisher.test(spark.cog$methyl_effect, spark.cog$derived_cog_impair)
spark.cog %>% 
  mutate(methyl_effect = ifelse(methyl_effect==1, "yes", "no"), 
         derived_cog_impair = ifelse(derived_cog_impair==1, "yes", "no")) %>%
  ggplot(aes(methyl_effect, fill=derived_cog_impair))+
  geom_bar() +
  scale_fill_manual(values = boxplot.colors, name = "cognitive impairment")+
  labs(subtitle = paste0("Fisher's Exact test", "\n",
                         "OR: ", round(t$estimate[[1]],4), 
                         "\npvalue: ", round(t$p.value,4))) +
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5))

###
# supplementary figures ---------------------------------------------------
# barplot for SPARM RM answers
pgs.rm %>% 
  pivot_longer(cols = c(starts_with("ch_med"), -contains("PM")), names_to = "q", values_to = "answer") %>%
  select(IID, q, answer) %>%
  drop_na() %>%
  mutate(answer = ifelse(answer==1, "yes", "no")) %>%
  ggplot(aes(x=answer))+
  geom_bar(position = "identity") +
  facet_wrap("q", scales = "free_y", ncol = 4)+
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5))
# histogram for PGS distribution
pgs.rm %>% 
  pivot_longer(cols = c(colnames(spark.pgs)[-1]), names_to = "pgs", values_to = "score") %>%
  select(IID, pgs, score) %>%
  drop_na() %>%
  ggplot(aes(x=score))+
  geom_histogram() +
  facet_wrap("pgs", scales = "free_y", ncol = 3)+
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5))
####
# supplementary tables ----------------------------------------------------
write_csv(inner_join(mph.scq%>%select(IID, interview_age, sex), pgs.rm) %>%
            # pivot_longer(cols = starts_with("ch_med"), names_to = "question", values_to = "answer") %>%
            group_by(sex) %>%
            summarise(avg = mean(interview_age), 
                      sd = sd(interview_age), 
                      min = min(interview_age), 
                      max = max(interview_age),
                      count = n()),
          file = "figs/paper/tmp/spark-mph-scq-data-stats.csv")
write_csv(do.call(rbind,
                  lapply(mph.scq %>% select(starts_with("q")), 
                         function(x) {
                           t <- fisher.test(mph.scq$methyl_effect, x)
                           data.frame(pval = t$p.value,
                                      OR = t$estimate[[1]],
                                      confint_min = as.data.frame(t$conf.int)[1,1],
                                      confint_max = as.data.frame(t$conf.int)[2,1])})) %>%
            rownames_to_column("question"),
          file = "figs/paper/tmp/spark-mph-scq-fisher-data-stats.csv")
####