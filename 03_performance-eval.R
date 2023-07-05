################################################################################
#               models performance for drug response prediction                #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(jtools)
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response"
setwd(project.dir)
################################################################################
# combine models output 
models <- data.frame(file = list.files("data/derivatives/m-outputs/spark", pattern = "model")) %>% 
  mutate(model_d = sub("model-", "", sub(".rds", "", file))) %>%
  mutate(model = paste0("model_", c(1:length(file)))) %>%
  mutate(weights_source = sub("-.*", "", model_d))
m.out <- list()
for (i in 1:nrow(models)) {
  # i=1
  m.o <- read_rds(paste0("data/derivatives/m-outputs/spark/", models$file[i])) %>%
    column_to_rownames("IID")
  colnames(m.o)[1] <- models$model[i]
  m.out[[i]] <- m.o
}
m.out <- do.call(cbind, m.out) %>%
  rownames_to_column("IID")
################################################################################
mph.spark <- read_tsv("data/derivatives/spark-rm-mph-effect.tsv")
meta <- left_join(mph.spark, m.out)
to.save <- list(meta = meta, models = models)
save(to.save, file = "data/derivatives/m-outputs/spark/combined-output.RData")
################################################################################
models.corr <- Hmisc::rcorr(meta%>%select(mph_effect,starts_with("model"))%>%as.matrix())
# what if we kept males only?
# models.corr <- Hmisc::rcorr(meta%>%filter(sex=="M")%>%select(mph_effect,starts_with("model"))%>%as.matrix())
# models.corr <- Hmisc::rcorr(meta%>%filter(sex=="F")%>%select(mph_effect,starts_with("model"))%>%as.matrix())
p.r <- models.corr$r %>%
  as.data.frame() %>%
  rownames_to_column("effect") %>%
  pivot_longer(cols = colnames(models.corr$r), names_to = "model", values_to = "r")
p.pval <- models.corr$P %>%
  as.data.frame() %>%
  rownames_to_column("effect") %>%
  pivot_longer(cols = colnames(models.corr$P), names_to = "model", values_to = "pval")
p.ready <- left_join(models, inner_join(p.r, p.pval) %>% filter(effect == "mph_effect"))
p <- p.ready %>%
  ggplot(aes(x=effect, y=model, fill = r, label = ifelse(pval<0.05, paste0("r=",round(r,5),"\n p=",round(pval,4)), ""))) +
  geom_tile() +
  geom_text(size=3.5) +
  xlab("") +
  ylab("") +
  scale_fill_gradient2(high = redblu.col[1], low = redblu.col[2], mid = "white") + 
  guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5)) +
  facet_grid(rows = "weights_source", space = "free", scales = "free") +
  labs(caption = paste0("model_d reference: F1-F2-F3-F4-F5",
                        "\n",
                        "F1 (weights source): [tissue] frontal cortex BA4 or [celltype] Excitatory",
                        "\n",
                        "F2 (genes used): [targets] MPH gene targets or [all] all genes",
                        "\n",
                        "F3 (imputed_tx corrected for gen-PCs): [T] True or [F] False",
                        "\n",
                        "F4 (z-scaled distributions): [T] True or [F] False",
                        "\n",
                        "F5 (method to compute response): [1] Pearson correlation or [2] Jaccard Index"))
library(ggridges)
library(gridExtra)
t <- gridExtra::tableGrob(models%>%select(-file), theme = ttheme_minimal(),
                          rows = NULL)

p2 <- left_join(meta %>% 
  select(IID, mph_effect, starts_with("model")) %>%
  pivot_longer(starts_with("model"), names_to = "model", values_to = "response"), models) %>%
  ggplot(aes(y = model, x = response)) +
  geom_density_ridges() +
  geom_vline(xintercept = 0, color = "red", linetype = 2) +
  facet_grid(rows = "weights_source", space = "free", scales = "free")
cowplot::plot_grid(patchwork::wrap_plots(p,p2, widths = c(1,3)),t, nrow = 1, rel_widths = c(3,1))
# cowplot::plot_grid(patchwork::wrap_plots(p,p2, widths = c(1,3)),t, nrow = 1, rel_widths = c(3,1), labels = "only females are kept for analysis", label_size = 8)
################################################################################
meta %>%
  ggplot(aes(x=model_4))+
  geom_histogram()
################################################################################
y <- meta$mph_effect

jtools::export_summs(lm(y~meta$model_1),lm(y~meta$model_2),lm(y~meta$model_3),lm(y~meta$model_4),
                     lm(y~meta$model_5),lm(y~meta$model_6),lm(y~meta$model_7),lm(y~meta$model_8),
                     lm(y~meta$model_9),lm(y~meta$model_10),lm(y~meta$model_11),lm(y~meta$model_12),
                     lm(y~meta$model_13),lm(y~meta$model_14),lm(y~meta$model_15),lm(y~meta$model_16),
                     lm(y~meta$model_17))
meta %>%
  mutate(mph_effect = as.factor(mph_effect)) %>%
  ggplot(aes(x=mph_effect, y=model_4, fill = mph_effect)) +
  geom_boxplot() 
################################################################################
load("data/derivatives/m-outputs/spark/combined-output.RData")
# meta <- to.save[[1]]

meta <- read_rds("data/derivatives/spark-mph-sample-meta-predicted-resp-pgs.rds")
models <- to.save[[2]]
meta.binarized <- meta %>%
  mutate_at(.vars = paste0("model_", c(1:17)), function(x) ifelse(x>0, 0, 1))

left_join(meta.binarized %>%
  dplyr::select(IID, sex, mph_effect, starts_with("model")) %>%
  pivot_longer(starts_with("model_"), names_to = "model", values_to = "predicted_response") %>%
  mutate(predicted_response = as.factor(predicted_response)), models) %>%
  mutate(model = factor(model, levels = models$model))%>%
  mutate(mph_effect = as.factor(mph_effect)) %>%
  # filter(sex == "F") %>%
  # ggplot(aes(x=predicted_response, fill = predicted_response)) +
  ggplot(aes(x=predicted_response, fill = mph_effect)) +
  geom_bar() +
  # labs(subtitle = "only females included") +
  facet_wrap("model", scales = "free", ncol = 8, dir = "h")
################################################################################
meta %>% 
  dplyr::select(IID, sex, mph_effect, starts_with("model")) %>%
  pivot_longer(starts_with("model_"), names_to = "model", values_to = "predicted_response") %>%
  mutate(model = factor(model, levels = models$model))%>%
  filter(model != "model_17") %>%
  ggplot(aes(x=predicted_response, color = model)) +
  geom_density(position = "identity")
################################################################################
models.corr <- Hmisc::rcorr(meta%>%select(mph_effect,starts_with("model"))%>%as.matrix())
# what if we kept males only?
# models.corr <- Hmisc::rcorr(meta%>%filter(sex=="M")%>%select(mph_effect,starts_with("model"))%>%as.matrix())
p.r <- models.corr$r %>%
  as.data.frame() %>%
  rownames_to_column("effect") %>%
  pivot_longer(cols = colnames(models.corr$r), names_to = "model", values_to = "r")
p.pval <- models.corr$P %>%
  as.data.frame() %>%
  rownames_to_column("effect") %>%
  pivot_longer(cols = colnames(models.corr$P), names_to = "model", values_to = "pval")
p.ready <- left_join(models, inner_join(p.r, p.pval) %>% filter(effect == "mph_effect"))
p <- p.ready %>%
  ggplot(aes(x=effect, y=model, fill = r, label = ifelse(pval<0.05, paste0("r=",round(r,5),"\n p=",round(pval,4)), ""))) +
  geom_tile() +
  geom_text(size=3.5) +
  xlab("") +
  ylab("") +
  scale_fill_gradient2(high = redblu.col[1], low = redblu.col[2], mid = "white") + 
  guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5)) +
  facet_grid(rows = "weights_source", space = "free", scales = "free") +
  labs(caption = paste0("model_d reference: F1-F2-F3-F4-F5",
                        "\n",
                        "F1 (weights source): [tissue] frontal cortex BA4 or [celltype] Excitatory",
                        "\n",
                        "F2 (genes used): [targets] MPH gene targets or [all] all genes",
                        "\n",
                        "F3 (imputed_tx corrected for gen-PCs): [T] True or [F] False",
                        "\n",
                        "F4 (z-scaled distributions): [T] True or [F] False",
                        "\n",
                        "F5 (method to compute response): [1] Pearson correlation or [2] Jaccard Index"))
p

################################################################################
meta %>% 
  mutate(mph_effect = as.factor(mph_effect)) %>%
  ggplot(aes(x=scale(model_4, center = F), fill = mph_effect)) +
  geom_histogram()

imputed.tx <- read_rds("data/derivatives/spark-mph-samples-LC-raw-imputed-tx-Excitatory.rds")
tmp <- left_join(meta%>%select(IID, mph_effect), imputed.tx%>%rownames_to_column("IID"))

registerDoMC(6)
weights <- foreach(i = 3:ncol(tmp), .combine = rbind) %dopar% {
  # i=4
  gene <- colnames(tmp)[i]
  return(summary(glm(mph_effect~ ., data = tmp%>%select(mph_effect,i)))$coefficients[-1,] %>% 
           t() %>% as.data.frame() %>% 
           mutate(gene = colnames(tmp)[i]))
}
sig.genes <- weights %>% filter(`Pr(>|t|)`<0.05)
summary(glm(mph_effect ~ ., data = tmp %>% select(mph_effect, sig.genes$gene)))
################################################################################
meta.eur <- meta %>% filter(iWES_gen_pop=="EUR")
jtools::export_summs(glm(meta.eur$mph_effect ~ meta.eur$model_1),glm(meta.eur$mph_effect ~ meta.eur$model_2),glm(meta.eur$mph_effect ~ meta.eur$model_3),glm(meta.eur$mph_effect ~ meta.eur$model_4),
                     glm(meta.eur$mph_effect ~ meta.eur$model_5),glm(meta.eur$mph_effect ~ meta.eur$model_6),glm(meta.eur$mph_effect ~ meta.eur$model_7),glm(meta.eur$mph_effect ~ meta.eur$model_8),
                     glm(meta.eur$mph_effect ~ meta.eur$model_9),glm(meta.eur$mph_effect ~ meta.eur$model_10),glm(meta.eur$mph_effect ~ meta.eur$model_11),glm(meta.eur$mph_effect ~ meta.eur$model_12),
                     glm(meta.eur$mph_effect ~ meta.eur$model_13),glm(meta.eur$mph_effect ~ meta.eur$model_14),glm(meta.eur$mph_effect ~ meta.eur$model_15),glm(meta.eur$mph_effect ~ meta.eur$model_16),
                     glm(meta.eur$mph_effect ~ meta.eur$model_17))
################################################################################
library(ggExtra)
ggMarginal(p = meta %>%
             mutate(mph_effect=as.factor(mph_effect)) %>%
             ggplot(aes(x=corrected_ADHD_Demontis, y=corrected_ADHD_PGC_2019, color = mph_effect)) +
             geom_point() +
             geom_smooth(), type = "density", groupColour = T)
ggMarginal(p = meta %>%
             mutate(mph_effect=as.factor(mph_effect)) %>%
             ggplot(aes(x=corrected_ADHD_PGC_2019, y=model_4, color = mph_effect)) +
             geom_point() +
             geom_smooth(method = "glm"), type = "density", groupColour = T)

################################################################################
# what if you corrected the predicted response for genetic PCs at the end?
geno.pcs <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/PCA/all/PCs.tsv") %>%
  filter(IID %in% meta$IID)
rownames(geno.pcs) = geno.pcs$IID
meta.2 <- left_join(meta, geno.pcs[,1:6])
meta.2$corrected_model_4 <- predict(glm(mph_effect~., data = meta.2%>%select(mph_effect, starts_with("pc_"))))
meta.2 %>%
  mutate(mph_effect=as.factor(mph_effect)) %>%
  ggplot(aes(x=mph_effect, y=corrected_model_4, fill = mph_effect)) +
  geom_violin()

################################################################################

################################################################################
tr.cbcl <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/cbcl/static_data/CBCL_syn-subscales_from-imputed.tsv")
inferred.cbcl <- read_tsv("data/derivatives/inferred-cbcl-dates-from-age-by-sleep-rm.tsv")
meta3 <- left_join(meta, left_join(tr.cbcl%>%select(IID, age_months, starts_with("syn")), inferred.cbcl))
meta3 %>% select(mph_effect, syn_attention, cbcl_age) %>%
  # mutate(mph_effect = as.factor(mph_effect)) %>%
  drop_na(syn_attention) %>%
  mutate(corr_attention = residuals(lm(syn_attention ~ cbcl_age))) %>%
  ggplot(aes(y=mph_effect, x=corr_attention)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor()
library(ggpubr)
summary(lm(meta3$syn_attention ~ meta3$sex))
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# May 22nd
# look at correlation between predicted response and spark-rm sleep data?
spark.rm <- read_tsv("data/rawdata/spark-rm/raw/ChildSurvey20181031.csv")
data.dict <- readxl::read_excel("data/rawdata/spark-rm/data-dict-child_edited.xlsx", sheet = 1) %>%
  mutate(notext = !grepl("text", ignore.case = T, x = values)) %>%
  filter(notext==T) %>%
  drop_na(values) %>%
  mutate(yn_flip = grepl("1 = No$", ignore.case = T, x = values)) %>%
  mutate(ynn_flip = grepl("2 = Not sure$", ignore.case = T, x = values)) %>%
  mutate(sub_1 = ifelse(data_col %in% c(paste0("q0", c(1:9), "_alias")), T, 
                        ifelse(data_col %in% c(paste0("q", c(20:29, 31:47, 62:69,72, 74:83,89,90:94, 104:108), "_alias")), T, F))) %>%
  mutate(flip_5 = ifelse(data_col %in% c(paste0("q", c(16:19, 60:61, 71,73, 88), "_alias")), T, F)) %>%
  mutate(flip_3 = ifelse(data_col %in% c(paste0("q", c(30), "_alias")), T, F)) %>%
  filter(!data_col %in% c("q59_sq01_alias", "q59_sq02_alias", "q87_sq01_alias", "q87_sq02_alias"))

# keep basic ID questions, age, sex, and other non text questions
spark.rm.filt <- spark.rm %>% select(c(1:4), data.dict$data_col)
################################################################################
# non-text y/n questions to modify their values
# make the flips of yn questions so that yes =1, no =0, NA
spark.rm.filt[spark.rm.filt == -666 | spark.rm.filt == -999] <- NA
for (i in 1:nrow(data.dict)){
  # i=1
  # print(i)
  q <- data.dict$data_col[i]
  index <- which(colnames(spark.rm.filt)==q)
  q_human <- data.dict$human_friendly_name[i]
  if (data.dict$yn_flip[i] == T) {
    spark.rm.filt <- spark.rm.filt %>% 
      rename(qq = which(colnames(spark.rm.filt)==q)) %>%
      mutate(qq = ifelse(is.na(qq), NA, ifelse(qq == 2, 0,1)))
  } else if (data.dict$ynn_flip[i] == T) {
    spark.rm.filt <- spark.rm.filt %>% 
      rename(qq = which(colnames(spark.rm.filt)==q)) %>%
      mutate(qq = ifelse(is.na(qq), NA, ifelse(qq %in% c(2,3), 0,1)))
  } else if (data.dict$sub_1[i] == T) {
    spark.rm.filt <- spark.rm.filt %>% 
      rename(qq = which(colnames(spark.rm.filt)==q)) %>%
      mutate(qq = ifelse(is.na(qq), NA, qq-1))
  } else if (data.dict$flip_5[i] == T) {
    spark.rm.filt <- spark.rm.filt %>% 
      rename(qq = which(colnames(spark.rm.filt)==q)) %>%
      mutate(qq = ifelse(is.na(qq), NA, ifelse(qq==1, 4, ifelse(qq==2, 3, ifelse(qq==4, 1, ifelse(qq==5,0,2))))))
  } else if (data.dict$flip_3[i] == T) {
    spark.rm.filt <- spark.rm.filt %>% 
      rename(qq = which(colnames(spark.rm.filt)==q)) %>%
      mutate(qq = ifelse(is.na(qq), NA, ifelse(qq==1, 2, ifelse(qq==2, 1, 0))))
  }
  # print(colnames(spark.rm.filt)[which(colnames(spark.rm.filt)=="qq")])
  colnames(spark.rm.filt)[index] <- q_human
  # print(colnames(spark.rm.filt)[which(colnames(spark.rm.filt)==q_human)])
}
write_tsv(spark.rm.filt, "data/rawdata/spark-rm/child_rmatch_ME.tsv")
################################################################################
temp <- inner_join(meta, spark.rm.filt)
# filter them to keep methyl_curr = 1 (currently taking MPH)
curr.mph <- temp %>% filter(ch_med_methyl_curr==1)
nas <- data.frame(q = colnames(curr.mph%>%select(starts_with("ch")))) %>%
  mutate(na = colSums(is.na(curr.mph[,q]))) %>% 
  filter(na < 450)
t <- corr.table(x = curr.mph %>% select(paste0("model_", c(1,3,5,7,9)), mph_effect) %>%
             mutate(model_1 = -model_1, 
                    model_3 = -model_3,
                    model_5 = -model_5,
                    model_7 = -model_7,
                    model_9 = -model_9), y = curr.mph%>%select(starts_with("ch"))%>%select(nas$q)) %>%
  filter(V1 %in% c(paste0("model_", c(1,3,5,7,9)), "mph_effect")) %>%
  filter(!V2 %in% c(paste0("model_", c(1:24)), "ch_med_methyl_effect", "mph_effect")) %>% 
  filter(!is.na(r))
sig.q <- t%>%filter(pval<0.05)%>%distinct(V2)%>%arrange(V2)
p10 <- t %>% filter(V2 %in% sig.q$V2) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label=ifelse(pval<0.05, paste0("r: ", round(r, 3), ", ", "p: ", round(pval, 4)), ""))) +
  geom_tile() +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1]) +
  geom_text(size=2.5, nudge_x = 0, nudge_y = 0) +
  my.guides +
  xlab("") + ylab("")
p11 <- curr.mph %>% 
  select(sig.q$V2) %>%
  pivot_longer(starts_with("ch"), names_to = "question", values_to = "answer") %>%
  drop_na(answer) %>%
  mutate(answer = as.factor(answer)) %>%
  ggplot(aes(x=answer)) +
  geom_bar() +
  facet_wrap("question", ncol = 4, scales = "free")
patchwork::wrap_plots(p10,p11, ncol = 2)
################################################################################
################################################################################
################################################################################
################################################################################
################# correlation between mph_effect and other data ################
################################################################################
################################################################################
################################################################################
spark.mph <- read_tsv("data/derivatives/spark-rm-mph-effect.tsv")
meta <- left_join(spark.mph, read_rds("data/derivatives/m-outputs/spark/all-spark-samples/model-celltype-all-TRUE-TRUE-1-TT.rds") %>%
  rename(model_3 = m))
################################################################################
# look at predicted response for MPH in spark and their scq
spark.scq <- read_csv("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12/scq_2022-12-12.csv")
mph.scq <- inner_join(meta, spark.scq%>%rename(IID=subject_sp_id)%>%select(-sex))

p1.bar <- mph.scq %>% select(IID, starts_with("q"), final_score) %>%
  pivot_longer(cols = c(starts_with("q"), final_score), names_to = "var", values_to = "val") %>%
  mutate(val = as.factor(val)) %>%
  ggplot(aes(x= val)) +
  geom_bar() +
  facet_wrap("var", scales = "free", ncol = 4)


p1 <- corr.table(mph.scq%>%select(mph_effect, model_3)%>%mutate(model_3 = -model_3), mph.scq%>%select(starts_with("q"), final_score)) %>%
  filter(V1 %in% c("mph_effect")) %>%
  filter(!V2 %in% c("mph_effect", "model_3")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label=ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 3)), ""))) +
  geom_tile()+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  geom_text(size=3) +
  labs(caption = paste0("n(samples): ", length(unique(mph.scq$IID))), title = "assessment: scq") +
  xlab("")+ylab("")+my.guides
################################################################################
# look at predicted response for MPH in spark and their iq
spark.iq <- read_csv("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12/iq_2022-12-12.csv")
mph.iq <- inner_join(meta, spark.iq%>%select(IID=subject_sp_id, test_year, age_test_date_months, ends_with("score"), ends_with("ratio")))

p2.bar <- mph.iq %>% select(IID, ends_with("score")) %>%
  pivot_longer(cols = c(ends_with("score")), names_to = "var", values_to = "val") %>%
  ggplot(aes(x= val)) +
  geom_histogram() +
  facet_wrap("var", scales = "free", nrow = 4)

p2 <- corr.table(mph.iq%>%select(mph_effect, model_3)%>%mutate(model_3 = -model_3), mph.iq%>%select(ends_with("score"))) %>%
  filter(V1 %in% c("mph_effect")) %>%
  filter(!V2 %in% c("mph_effect", "model_3")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label=ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 3)), ""))) +
  geom_tile()+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  geom_text(size=3) +
  labs(caption = paste0("n(samples): ", length(unique(mph.iq$IID))), title = "assessment: iq") +
  xlab("")+ylab("")+my.guides


# instead of the iq test use the predicted cog-impair
spark.cogimpair <- read_csv("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12/predicted_iq_experimental_2022-12-12.csv")
mph.cogimpair <- inner_join(meta, spark.cogimpair%>%select(IID=subject_sp_id, derived_cog_impair)%>%drop_na(derived_cog_impair))

p22.bar <- mph.cogimpair %>% select(IID, derived_cog_impair) %>%
  pivot_longer(cols = c(derived_cog_impair), names_to = "var", values_to = "val") %>%
  mutate(val = as.factor(val)) %>%
  ggplot(aes(x= val)) +
  geom_bar() +
  facet_wrap("var", scales = "free", nrow = 4)+xlab("")

ft <- fisher.test(x = mph.cogimpair$mph_effect, y = mph.cogimpair$derived_cog_impair)

p22 <- mph.cogimpair %>%
  mutate(mph_effect = as.factor(mph_effect)) %>%
  mutate(derived_cog_impair = as.factor(derived_cog_impair)) %>%
  ggplot(aes(x=mph_effect, fill=derived_cog_impair)) +
  # geom_point(position = "jitter")+
  # geom_smooth(method = "lm")
  geom_bar()+
  scale_fill_manual(values = boxplot.colors) +
  labs(caption = paste0("n(samples): ", nrow(mph.cogimpair), "\n",
                        "Fisher's Exact Test for Count Data", "\n",
                        "p-value: ", round(ft$p.value,4), "\n",
                        "confint: ", round(ft$conf.int[1],4), "    ", round(ft$conf.int[2],4), "\n",
                        "OR: ", round(ft$estimate,4)))

p22.bar+p22
################################################################################
# look at predicted response for MPH in spark and their rbs-r
spark.rbsr <- read_csv("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12/rbsr_2022-12-12.csv")
mph.rbsr <- inner_join(meta, spark.rbsr%>%select(IID=subject_sp_id, age_at_eval_months, ends_with("score")))

p3.bar <- mph.rbsr %>% select(IID, ends_with("score")) %>%
  pivot_longer(cols = c(ends_with("score")), names_to = "var", values_to = "val") %>%
  ggplot(aes(x= val)) +
  geom_bar() +
  facet_wrap("var", scales = "free", nrow = 4)


p3 <- corr.table(mph.rbsr%>%select(mph_effect, model_3)%>%mutate(model_3 = -model_3), mph.rbsr%>%select(ends_with("score"))) %>%
  filter(V1 %in% c("mph_effect")) %>%
  filter(!V2 %in% c("mph_effect", "model_3")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label=ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 3)), ""))) +
  geom_tile()+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  geom_text(size=3) +
  labs(caption = paste0("n(samples): ", length(unique(mph.rbsr$IID))), title = "assessment: rbsr")+
  xlab("")+ylab("")+my.guides
################################################################################
# look at predicted response for MPH in spark and their bg-history child
spark.bgchild <- read_csv("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12/background_history_child_2022-12-12.csv")
mph.bgchild <- inner_join(meta, spark.bgchild%>%select(IID=subject_sp_id, age_at_eval_months, 
                                                       mother_highest_education, father_highest_education,
                                                       language_age_level,sped_y_n, sped_aide, sped_speech, sped_ot, sped_pt, 
                                                       sped_behavior, sped_soc_skills, intervention_behavioral, 
                                                       intervention_speech_language, intervention_ot_sensory))

encode.education <- function(x) {
  temp <- data.frame(x = x, encoded = NA) %>%
    mutate(encoded = ifelse(x == "did_not_attend_high_school", 1,
                            ifelse(x == "some_high_school", 2,
                                   ifelse(x == "ged_diploma", 3,
                                          ifelse(x == "trade_school", 4,
                                                 ifelse(x == "high_school_graduate", 5,
                                                        ifelse(x == "some_college", 6,
                                                               ifelse(x == "associate_degree", 7,
                                                                      ifelse(x == "baccalaureate_degree", 8,
                                                                             ifelse(x == "graduate_or_professional_degree", 8, NA))))))))))
  return(temp$encoded)
}
encode.to.age <- function(x) {
  temp <- data.frame(x=x, encoded = NA) %>%
    mutate(encoded = ifelse(x == "signif_below_age", 1,
                            ifelse(x == "slight_below_age", 2,
                                   ifelse(x == "at_age", 3,
                                          ifelse(x == "above_age", 4, NA)))))
  return(temp$encoded)
}
mph.bgchild <- mph.bgchild %>%
  mutate_at(.vars = c("mother_highest_education", "father_highest_education"), .funs = function(x) encode.education(x)) %>%
  mutate_at(.vars = "language_age_level", .funs = function(x) encode.to.age(x)) %>%
  mutate_at(.vars = c("sped_y_n", "sped_aide", "sped_speech", "sped_ot", "sped_pt", 
                      "sped_behavior", "sped_soc_skills", "intervention_behavioral", 
                      "intervention_speech_language", "intervention_ot_sensory"), .funs = function(x) ifelse(is.na(x), 0, x))

p4.bar <- mph.bgchild %>% select(IID, mother_highest_education, father_highest_education,
                       language_age_level,sped_y_n, sped_aide, sped_speech, sped_ot, sped_pt, 
                       sped_behavior, sped_soc_skills, intervention_behavioral, 
                       intervention_speech_language, intervention_ot_sensory) %>%
  pivot_longer(cols = c(mother_highest_education, father_highest_education,
                        language_age_level,sped_y_n, sped_aide, sped_speech, sped_ot, sped_pt, 
                        sped_behavior, sped_soc_skills, intervention_behavioral, 
                        intervention_speech_language, intervention_ot_sensory), names_to = "var", values_to = "val") %>%
  mutate(val = as.factor(val)) %>%
  ggplot(aes(x= val)) +
  geom_bar() +
  facet_wrap("var", scales = "free", nrow = 4)


p4 <- corr.table(mph.bgchild%>%select(mph_effect, model_3)%>%mutate(model_3 = -model_3), mph.bgchild%>%select(ends_with("education"),
                                                                                                        ends_with("level"), 
                                                                                                        starts_with("sped"),
                                                                                                        starts_with("intervention"))) %>%
  filter(V1 %in% c("mph_effect")) %>%
  filter(!V2 %in% c("mph_effect", "model_3")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label=ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 3)), ""))) +
  geom_tile()+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  geom_text(size=3) +
  labs(caption = paste0("n(samples): ", length(unique(mph.bgchild$IID))), title = "assessment: background history child")+
  xlab("")+ylab("")+my.guides
################################################################################
# look at predicted response for MPH in spark and their bg-history adult
spark.bgadult <- read_csv("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12/background_history_adult_2022-12-12.csv")
mph.bgadult <- inner_join(meta %>% mutate(IID = InformantID) %>% select(-sex), spark.bgadult%>%select(IID=subject_sp_id, age_at_eval_months, 
                                                       highest_level_education,
                                                       starts_with("support"),
                                                       starts_with("adult")))
encode.rarely.to.most <- function(x) {
  temp <- data.frame(x=x, encoded = NA) %>%
    mutate(encoded = ifelse(x == "rarely", 1,
                            ifelse(x == "sometimes", 2,
                                   ifelse(x == "most", 3, NA))))
  return(temp$encoded)
}
mph.bgadult <- mph.bgadult %>%
  mutate_at(.vars = c("highest_level_education"), .funs = function(x) encode.education(x)) %>%
  mutate_at(.vars = "supports_how_much", .funs = function(x) encode.rarely.to.most(x)) %>%
  mutate_at(.vars = c("supports" , "supports_daily_living", "supports_employ", "supports_health", "supports_finances",              
                      "supports_community_transp","supports_life","supports_social","supports_counseling",
                      "supports_housing","supports_how_much","supports_daily_living_need","supports_employ_need",
                      "supports_health_need","supports_finances_need","supports_community_transp_need","supports_life_need",
                      "supports_social_need","supports_counseling_need","supports_housing_need","adult_needs_financial",
                      "adult_needs_employ","adult_needs_housing","adult_needs_ed","adult_needs_therapy_access",
                      "adult_needs_daily_living","adult_needs_med_tx","adult_needs_non_med_tx","adult_needs_social_access",
                      "adult_needs_parity","adult_needs_research_cause","adult_needs_acceptance"), .funs = function(x) ifelse(is.na(x), 0, x))
p5.bar <- mph.bgadult %>% select(IID, highest_level_education, starts_with("support"), starts_with("adult")) %>%
  pivot_longer(cols = c(highest_level_education, starts_with("support"), starts_with("adult")), names_to = "var", values_to = "val") %>%
  mutate(val = as.factor(val)) %>%
  ggplot(aes(x= val)) +
  geom_bar() +
  facet_wrap("var", scales = "free", ncol = 4)

p5 <- corr.table(mph.bgadult%>%select(mph_effect, model_3)%>%mutate(model_3 = -model_3), mph.bgadult%>%select(ends_with("education"),
                                                                                                        starts_with("support"),
                                                                                                        starts_with("adult"))) %>%
  filter(V1 %in% c("mph_effect")) %>%
  filter(!V2 %in% c("mph_effect", "model_3")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label=ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 3)), ""))) +
  geom_tile()+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  geom_text(size=3) +
  labs(caption = paste0("n(samples): ", length(unique(mph.bgadult$IID))), title = "assessment: background history adult")+
  xlab("")+ylab("")+my.guides
################################################################################
# look at predicted response for MPH in spark and their vineland
spark.vineland <- read_csv("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12/vineland-3_2022-12-12.csv", na = c("-9999.0"))
mph.vineland <- inner_join(meta, spark.vineland%>%select(IID=subject_sp_id, ends_with("est")))

p6.bar <- mph.vineland %>% select(IID, ends_with("est")) %>%
  pivot_longer(cols = c(ends_with("est")), names_to = "var", values_to = "val") %>%
  ggplot(aes(x= val)) +
  geom_histogram() +
  facet_wrap("var", scales = "free", ncol = 4)


p6 <- corr.table(mph.vineland%>%select(mph_effect, model_3)%>%mutate(model_3 = -model_3), mph.vineland%>%select(ends_with("est"))) %>%
  filter(V1 %in% c("mph_effect")) %>%
  filter(!V2 %in% c("mph_effect", "model_3")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label=ifelse(pval<0.3, paste0("r: ", round(r, 3), ",  p: ", round(pval, 3)), ""))) +
  geom_tile()+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  geom_text(size=3) +
  labs(caption = paste0("n(samples): ", length(unique(mph.vineland$IID))), title = "assessment: vineland") +
  xlab("")+ylab("")+my.guides
################################################################################

# combine plots
patchwork::wrap_plots(p1, p1.bar, widths = c(1,3))
patchwork::wrap_plots(patchwork::wrap_plots(p2, p2.bar, widths = c(1,2)), patchwork::wrap_plots(p22, p22.bar), ncol = 1, heights = c(2,1))
patchwork::wrap_plots(p3, p3.bar, widths = c(1,2))
patchwork::wrap_plots(p4, p4.bar, widths = c(1,2))
patchwork::wrap_plots(p5, p5.bar, widths = c(1,2))
patchwork::wrap_plots(p6, p6.bar, widths = c(1,2))
################################################################################
