################################################################################
#                               pre-model lookups                              #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(Hmisc)
library(jtools)
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response"
setwd(project.dir)
################################################################################
spark.rm <- read_tsv("data/rawdata/spark-rm/raw/ChildSurvey20181031.csv")
data.dict <- readxl::read_excel("data/rawdata/spark-rm/data-dict-child_edited.xlsx", sheet = 1)
# keep basic ID questions, age, sex, and questions starting of q95
spark.rm.filt <- spark.rm[,c(1:4, 160:ncol(spark.rm))]
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
################################################################################
nas <- left_join(data.frame(question = colnames(spark.rm.filt), count = colSums(is.na(spark.rm.filt))),
                 data.dict%>%select(question = data_col, human_friendly_name)) %>%
  mutate(human_friendly_name = ifelse(is.na(human_friendly_name), question, human_friendly_name))
nas %>%
  drop_na() %>%
  ggplot(aes(x=count, y=human_friendly_name)) +
  geom_bar(stat = "identity")
################################################################################
# get samples and questions of interest
colnames(spark.rm.filt) <- nas$human_friendly_name
nas <- nas %>%
  mutate(int = grepl("ID", human_friendly_name)) %>%
  mutate(int = ifelse(int==F,grepl("age", human_friendly_name),int)) %>%
  mutate(int = ifelse(int==F,grepl("yn", human_friendly_name),int)) %>%
  mutate(int = ifelse(int==F,grepl("curr", human_friendly_name),int)) %>%
  mutate(int = ifelse(int==F,grepl("effect", human_friendly_name),int)) %>%
  mutate(int = ifelse(int==F,grepl("help", human_friendly_name),int)) %>%
  mutate(int = ifelse(question %in% c("q104_alias", "q105_alias", "q106_alias", "q107_alias", "q108_alias"), T, int))
sample <- spark.rm.filt %>% 
  select(nas$human_friendly_name[which(nas$int==T)]) 
# look at drug effectiveness (methylphenidate and amphetamines)
# effectiveness to age

# age not predictive of effectiveness
summary(lm(sample$ch_med_methyl_effect ~ sample$Dependent_age_at_eval_start_months))
# infromant age not predictive
summary(lm(sample$ch_med_methyl_effect ~ sample$Informant_age_at_eval_start_months))


############################ look at correlations between effectiveness and other questions
tt <- Hmisc::rcorr(sample%>%select(-c(contains("ID")))%>%as.matrix())
p.r <- tt$r %>%
  as.data.frame() %>%
  rownames_to_column("q") %>%
  pivot_longer(cols = colnames(tt$r), names_to = "q2", values_to = "r")
p.pval <- tt$P %>%
  as.data.frame() %>%
  rownames_to_column("q") %>%
  pivot_longer(cols = colnames(tt$P), names_to = "q2", values_to = "pval")
p.ready <- inner_join(p.r, p.pval)
p1 <- p.ready %>%
  mutate(effect = grepl("effect", q)) %>%
  mutate(effect = ifelse(effect==F, grepl("help", q), effect)) %>%
  # mutate(r=as.numeric(r)) %>%
  # drop_na(r) %>%
  filter(effect==T) %>%
  ggplot(aes(x=q, y=q2, fill = r, label = ifelse(pval<0.05, round(r,5), ""))) +
  geom_tile() +
  geom_text(size=2) +
  xlab("") +
  ylab("") +
  scale_fill_gradient(high = redblu.col[1], low = redblu.col[2]) +
  guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5))
p2 <- sample %>%
  mutate(ch_med_methyl_effect = as.factor(ch_med_methyl_effect)) %>%
  mutate(ch_med_methyl_curr = as.factor(ch_med_methyl_curr)) %>%
  drop_na(ch_med_methyl_curr, ch_med_methyl_effect) %>%
  ggplot(aes(x=ch_med_methyl_curr , fill = ch_med_methyl_effect)) +
  geom_bar()
patchwork::wrap_plots(p1,p2, nrow = 2, heights = c(4,1))
################################################################################
# look at spark mph-samples
mph.samples <- sample %>% 
  select(InformantID, Informant_age=Informant_age_at_eval_start_months, 
         ParticipantID, Dependent_age=Dependent_age_at_eval_start_months,
         mph_effect=ch_med_methyl_effect) %>%
  drop_na(mph_effect)
metadata <- read.delim("/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/master.tsv") %>%
  filter(IID %in% mph.samples$ParticipantID)
mph.samples <- inner_join(mph.samples%>%mutate(IID=ParticipantID), metadata)

summary(lm(mph.samples$mph_effect ~ mph.samples$Informant_age))
summary(lm(mph.samples$mph_effect ~ mph.samples$Dependent_age))
summary(lm(mph.samples$mph_effect ~ mph.samples$sex))

library(ggExtra)
ggMarginal(mph.samples %>% ggplot(aes(x=Dependent_age, y=mph_effect)) +geom_point()+geom_smooth())

write_tsv(mph.samples, "data/derivatives/spark-rm-mph-effect.tsv")
# writeLines(mph.samples$IID, "data/derivatives/mph-sample-IDs")
# writeLines(mph.samples$FID, "data/derivatives/mph-sample-FIDs")
################################################################################
# run GCTA to extract samples of interest 
bed.fam <- read_table("/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/genotypes/final.fam", col_names = F)
bed.fam <- bed.fam %>% filter(X2 %in% mph.samples$ParticipantID)
write_delim(bed.fam, "data/derivatives/mph-sample-IDs.fam", delim = "\t", col_names = F)
gcta_command <- paste(
  "/Dedicated/jmichaelson-wdata/trthomas/bin/gcta/gcta_1.93.2beta/gcta64", # GCTA executable
  "--bfile", "/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/genotypes/final", #PLINK files
  "--keep",  paste0(project.dir, "/data/derivatives/mph-sample-IDs.fam"), #list of SNPs
  "--thread-num", 8, # threads
  "--recode",
  "--out", paste0(project.dir, "/data/derivatives/spark-mph-samples-genotypes"),
  sep = " ")
gcta_command
################################################################################
p3 <- mph.samples %>% 
  mutate(sex = as.factor(sex)) %>%
  ggplot(aes(x=sex)) +
  geom_bar()
p4 <- mph.samples %>%
  ggplot(aes(x=Dependent_age)) +
  geom_histogram()
p5 <- mph.samples %>%
  mutate(mph_effect = as.factor(mph_effect)) %>%
  ggplot(aes(x=mph_effect, fill=sex)) +
  geom_bar() +
  scale_fill_manual(values = redblu.col)

patchwork::wrap_plots(p3,p4,p5,ncol = 3)
################################################################################
# get pgs for spark
mph.spark <- read_tsv("data/derivatives/spark-rm-mph-effect.tsv")
geno.pcs <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/PCA/all/PCs.tsv")
pgs.meta <- read_delim("/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/PGS/clean_PGS_names.tsv") %>%
  mutate(dup = ifelse(duplicated(PGS_shortname), 1, 0))
pgs <- read_delim("/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/PGS/all_PGS.tsv", name_repair = "minimal") %>%
  filter(PGS_name %in% pgs.meta$PGS_name) %>%
  dplyr::select(IID, PGS_name, PGS_raw) %>%
  pivot_wider(names_from = "PGS_name", values_from = "PGS_raw")
pgs.filt <- left_join(pgs, geno.pcs[,1:6])
x <- pgs.filt[,1:(ncol(pgs.filt)-5)]
y <- pgs.filt[,(ncol(pgs.filt)-5):ncol(pgs.filt)]
pgs.corrected <- cbind(x[,1], lapply(x[,2:ncol(x)], function(z) 
  scale(residuals(lm(as.numeric(z) ~ y$pc_01 + y$pc_02 + y$pc_03 + y$pc_04 + y$pc_05)), center = T, scale = T)))

colnames(pgs.corrected)[2:(ncol(pgs.corrected))] <- paste0("corrected_", colnames(pgs)[2:(ncol(pgs.corrected))])
pgs.filt <- left_join(pgs.corrected, pgs.filt) %>% filter(IID %in% mph.spark$IID)
write_tsv(pgs.filt, "data/derivatives/spark-rm-pgs.tsv")
# write_tsv(pgs.corrected, "data/derivatives/spark-abcd-corrected-pgs.tsv")
###############################################################################
# run GCTA to extract samples of interest of all genotypes, not just the imputed ones
# the file is so big, so include variant names of interest
tissue <- "Brain_Frontal_Cortex_BA9"
gcta_command <- paste(
  "/Dedicated/jmichaelson-wdata/trthomas/bin/gcta/gcta_1.93.2beta/gcta64", # GCTA executable
  "--bfile", "/Dedicated/jmichaelson-wdata/lcasten/spark/language_impairment_gwas/data/geno/merged/merged_hg19", #PLINK files
  "--keep",  paste0(project.dir, "/data/derivatives/mph-sample-IDs.fam"), 
  "--extract", paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/tx-imputation/UTMOST-GTEx-model-weights/tmp/rsid-ID02-UTMOST-for-",tissue), #list of SNPs
  "--thread-num", 8, # threads
  "--recode",
  "--out", paste0(project.dir, "/data/derivatives/spark-mph-samples-genotypes-LC-merged-", tissue),
  sep = " ")
gcta_command
################################################################################
# run GCTA to extract samples of interest of all genotypes, not just the imputed ones
# the file is so big, so include variant names of interest for celltype
gcta_command <- paste(
  "/Dedicated/jmichaelson-wdata/trthomas/bin/gcta/gcta_1.93.2beta/gcta64", # GCTA executable
  "--bfile", "/Dedicated/jmichaelson-wdata/lcasten/spark/language_impairment_gwas/data/geno/merged/merged_hg19", #PLINK files
  "--keep",  paste0(project.dir, "/data/derivatives/mph-sample-IDs.fam"), 
  "--extract", paste0("/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/Excitatory-ID_37-FDR-sig"), #list of SNPs
  "--thread-num", 8, # threads
  "--recode",
  "--out", paste0(project.dir, "/data/derivatives/spark-mph-samples-genotypes-LC-merged-Excitatory-FDR-sig"),
  sep = " ")
gcta_command

################################################################################
################################################################################
# check correlation between MPH effectiveness and ASR
spark.asr <- read_tsv("data/rawdata/ASR_combined_2e.tsv")
tmp <- cbind(spark.asr[,1:8], do.call(cbind, lapply(spark.asr[,9:ncol(spark.asr)], function(x) as.numeric(x)))) %>%
  as.data.frame() %>%
  select(where(not_all_na))
mph.asr <- inner_join(mph.spark %>% mutate(IID = InformantID) %>% select(-sex), tmp)
plot.t <- corr.table(x = mph.asr%>%select(mph_effect), y = mph.asr%>%select(starts_with("q"))) %>%
  filter(V1 == "mph_effect") %>%
  filter(V2 != "mph_effect") 
p1 <- plot.t %>% 
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, paste0("r: ", r, "\n", "p: ",pval), ""))) +
  geom_tile() +
  geom_text(size = 2) +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1]) +
  my.guides + 
  xlab("") + ylab("")
# subscales instead?
spark.asr.subscales <- read_tsv("data/rawdata/ASR_subscales.tsv")
mph.asr.sc <- inner_join(mph.spark %>% mutate(IID=InformantID), spark.asr.subscales)
p2<- corr.table(x = mph.asr.sc%>%select(mph_effect), y = mph.asr.sc%>%select(starts_with("ASR"))) %>%
  filter(V1 == "mph_effect") %>%
  filter(V2 != "mph_effect") %>% 
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, paste0("r: ", r, "\n", "p: ",pval), ""))) +
  geom_tile() +
  geom_text(size = 2) +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1]) +
  my.guides + 
  xlab("") + ylab("")
patchwork::wrap_plots(p1,p2, nrow = 1)
################################################################################
# check correlation between mph_effect vector and subject communication?
spark.scq <- read_csv("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12/scq_2022-12-12.csv")
mph.scq <- inner_join(mph.spark, spark.scq%>%rename(IID=subject_sp_id)%>%select(-sex))
corr.table(mph.scq%>%select(mph_effect), mph.scq%>%select(starts_with("q"), final_score)) %>%
  filter(V1 == "mph_effect") %>%
  filter(V2 != "mph_effect") %>%
  ggplot(aes(x=V1, y=V2, fill = r, label=ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 3)), ""))) +
  geom_tile()+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  geom_text(size=3) +
  xlab("")+ylab("")+my.guides
################################################################################
# is mph_effect reported by parents is detectable in CBCL?
# indep section
mph.spark <- read_tsv("data/derivatives/spark-rm-mph-effect.tsv")
cbcl.scales <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/cbcl-scales.rds")
spark.cbcl <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/cbcl/static_data/CBCL_syn-subscales_from-imputed.tsv") %>%
  filter(IID %in% mph.spark$IID)
meta <- inner_join(mph.spark%>%select(-c(sex, cohort)), spark.cbcl)
summary(lm(syn_attention ~ mph_effect, data = meta))
p1.bar <- meta %>%
  select(IID, starts_with("syn"), total) %>%
  pivot_longer(cols = c(starts_with("syn"), total), names_to = "var", values_to = "val") %>%
  mutate(var2 = sub("syn_","", var)) %>%
  ggplot(aes(x=val)) + 
  geom_histogram()+
  facet_wrap("var2", scales = "free", ncol = 3)
p1 <- corr.table(meta%>% select(mph_effect), meta%>%select(starts_with("syn"), total)) %>%
  filter(V1 == "mph_effect") %>%
  filter(V2 != "mph_effect") %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.1, paste0("r: ", round(r, 3), "\n","p: ", round(pval, 3)), ""))) +
  geom_tile() +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1]) +
  geom_text(size=3) +
  my.guides + xlab("") + ylab("") +
  labs(caption = paste0("n(samples): ", nrow(meta)))
patchwork::wrap_plots(p1, p1.bar, widths = c(1,2))

spark.cbcl <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/cbcl/static_data/CBCL_item-level_imputed.tsv") %>%
  filter(IID %in% mph.spark$IID)
meta <- inner_join(mph.spark%>%select(-c(sex, cohort)), spark.cbcl)

p2.bar <- meta %>%
  select(IID, starts_with("q")) %>%
  pivot_longer(cols = c(starts_with("q")), names_to = "var", values_to = "val") %>%
  mutate(val = as.factor(val)) %>%
  mutate(var2 = sub("_.*","", var)) %>%
  ggplot(aes(x=val)) + 
  geom_bar()+
  facet_wrap("var2", scales = "free_y", ncol = 6)
p2 <- left_join(corr.table(meta%>% select(mph_effect), meta%>%select(starts_with("q"), cbcl_total)) %>%
  filter(V1 == "mph_effect") %>%
  filter(V2 != "mph_effect"), cbcl.scales %>% rename(V2 = cbcl_item)) %>%
  drop_na(cbcl_syndrome) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 3)), ""))) +
  geom_tile() +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1]) +
  geom_text(size=3) +
  my.guides + xlab("") + ylab("") +
  labs(caption = paste0("n(samples): ", nrow(meta))) +
  facet_grid(rows  = vars(cbcl_syndrome), space = "free", scales = "free")  
patchwork::wrap_plots(p2, p2.bar, widths = c(1,3))
patchwork::wrap_plots(patchwork::wrap_plots(p1, p1.bar, widths = c(1,2)), patchwork::wrap_plots(p2, p2.bar, widths = c(1,3)), 
                      widths = c(1.5,3))
################################################################################
# check pgs correlations with mph_effect
# indep section
mph.spark <- read_tsv("data/derivatives/spark-rm-mph-effect.tsv")
spark.all <- left_join(mph.spark %>% select(-starts_with("corrected")), read_tsv("data/derivatives/spark-rm-pgs.tsv"))
corr.table(spark.all %>% select(mph_effect), spark.all %>% select(starts_with("corrected"))) %>%
  filter(V1 == "mph_effect", V2 != "mph_effect") %>%
  mutate(V2 = sub("corrected_", "", V2)) %>%
  ggplot(aes(x=V1, y = V2, fill = r, label=ifelse(pval<0.05, paste0("r: ", round(r,3), ",  p: ", round(pval,5)), ""))) +
  geom_tile() +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1]) +
  geom_text(size = 3)+
  my.guides + xlab("")+ylab("") +
  labs(caption = paste0("n(samples): ", nrow(spark.all)))
spark.all %>%
  mutate(mph_effect = as.factor(mph_effect)) %>%
  ggplot(aes(x = mph_effect, y = `corrected_cog_symbol_digit-UKB-2020`, fill = mph_effect)) +
  geom_boxplot(show.legend = F)+
  stat_compare_means(method = "t.test") +
  scale_fill_manual(values = boxplot.colors)
################################################################################

