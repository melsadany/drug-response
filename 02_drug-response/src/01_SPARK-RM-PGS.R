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
                     spark.pgs)
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
  labs(caption = paste0("n(samples): ", nrow(pgs.rm), "\n",
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
  labs(caption = paste0("n(samples): ", nrow(pgs.rm))) +
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5))
####
