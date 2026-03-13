################################################################################
#                                pre figures for SPARK                        #
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(heatmap3)
library(ggpubr);library(ggExtra);library(ggh4x)
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response")
setwd(project.dir)
################################################################################
spark.rm <- read_tsv("../data/rawdata/spark-rm/raw/ChildSurvey20181031.csv")
data.dict <- readxl::read_excel("../data/rawdata/spark-rm/data-dict-child_edited.xlsx", sheet = 1)
# keep basic ID questions, age, sex, and questions starting of q95
spark.rm.filt <- spark.rm[,c(1:4, 160:ncol(spark.rm))]
################################################################################
################################################################################
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
    else if(spark.rm.filt[j,i] == -666){spark.rm.filt[j,i] = NA}
    else if(spark.rm.filt[j,i] == -999){spark.rm.filt[j,i] = NA}
  }
}
spark.rm.filt[spark.rm.filt == -666 | spark.rm.filt == -999] <- NA
colnames(spark.rm.filt)[5:ncol(spark.rm.filt)] <- data.dict$human_friendly_name[157:nrow(data.dict)]
################################################################################
# correlation between pgs and taking meds in spark
spark.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv") %>%
  rename_all(.funs = function(x) sub("corrected_", "", x)) %>%
  select(IID, "ADHD-Demontis", contains("cog")&contains("UKB")) %>%
  rename_all(.funs = function(x) str_replace_all(x, "-UKB-2020", ""))
pgs.rm <- inner_join(spark.rm.filt %>% select(IID = ParticipantID, ends_with("yn"), ends_with("effect")),
                     spark.pgs)
corr.table(pgs.rm %>% select(colnames(spark.pgs)[-1],ends_with("yn"), ends_with("effect")),
           method = "spearman") %>%
  filter(V1 %in% colnames(spark.pgs)[-1], !V2 %in% colnames(spark.pgs)[-1]) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  mutate(V2 = sub("ch_med_", "", V2)) %>%
  filter(!grepl("PM", V2)) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "***", ifelse(pval<0.01, "**", ifelse(pval<0.05, "*", ""))))) +
  geom_tile()+
  geom_text(size = 3) +
  redblu.col.gradient()+my.guides+null_labs +
  labs(caption = paste0("n(samples): ", nrow(pgs.rm), "\n",
                        "* pval < 0.05", "\n", 
                        "** pval < 0.01", "\n", 
                        "*** FDR < 0.05"))

pgs.rm %>% pivot_longer(cols = c(ends_with("yn"), ends_with("effect"))) %>% 
  filter(!is.na(value)) %>% 
  mutate(value=ifelse(value==1,"yes","no"),name = sub("ch_med_","",name),
         type=case_when(grepl("_yn",name)~"taking the drug", 
                        grepl("_effect",name)~"drug effective", T~""), 
         name = sub("_yn|_effect","",name))%>%
  ggplot(aes(as.factor(value),`ADHD-Demontis`,fill=value)) +
  geom_violin(show.legend=F)+geom_boxplot(width=0.2,fill="white")+
  ggpubr::stat_compare_means(label.y.npc = 0.9,color="red")+
  scale_fill_manual(values=palette.1)+
  ggh4x::facet_nested(rows=vars(type,name),scales="free",nest_line=element_line(linewidth=0.6))+
  bw.theme+labs(x="",y="ADHD PGS (Demontis)")
ggsave2("../2026/figs/spark-ADHD-PGS-to-meds-taking-and-effect.pdf",4,18)
################################################################################
################################################################################


################################################################################


################################################################################

