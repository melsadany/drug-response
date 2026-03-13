################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/2026")
setwd(project.dir)
################################################################################
################################################################################
spark.dir <- correct_path("/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12")
################################################################################
################################################################################
# read the SCQ data, and correct the yes/no values
spark.scq.raw <- read_csv(paste0(spark.dir, "/scq_2022-12-12.csv"))
spark.scq <- spark.scq.raw %>% filter(missing_values == 0) %>%
  mutate_at(.vars = vars(starts_with("q")), .funs = function(x) ifelse(is.na(x), 0, x)) %>%
  rename(IID = subject_sp_id)
################################################################################
################################################################################
# read the child background data
spark.bg.child.raw <- read_csv(paste0(spark.dir, "/background_history_child_2022-12-12.csv"), na = c("888"))
spark.bg.child <- spark.bg.child.raw %>% select(IID = subject_sp_id, contains("age_level")) %>%
  drop_na()
################################################################################
################################################################################
# read the cognitive data 
spark.iq.raw <- read_csv(paste0(spark.dir, "/predicted_iq_experimental_2022-12-12.csv"))
spark.iq <- spark.iq.raw %>% select(IID = subject_sp_id, derived_cog_impair) %>%
  drop_na()
################################################################################
################################################################################
# read spark RM data from the eating/sleep survey
# keep the methylphenidate effectiveness question and flip the values
# the raw is coded as 1 for yeas, and 2 for no
rm.q.dict <- readxl::read_excel("../data/rawdata/spark-rm/data-dict-child_edited.xlsx", sheet = 1)
spark.rm <- read_tsv("../data/rawdata/spark-rm/raw/ChildSurvey20181031.csv", na = c("-666", "-999")) %>% 
  select(1:8, q102_sq03_alias,q102_sq04_alias) %>%
  rename(mph_effect = q102_sq03_alias) %>%
  filter(! is.na(mph_effect)) %>%
  mutate(mph_effect = ifelse(mph_effect == 2, 0, 1)) %>%
  rename(IID = ParticipantID)
################################################################################
################################################################################
# corr w scq 
# decide on what samples you want to keep. 
# you might keep all samples or only keep the ones with verbal communication with phrases==1

# all types of communicators
tmp <- inner_join(spark.rm%>%select(IID, mph_effect), 
                  spark.scq%>%select(IID, starts_with("corr"), starts_with("q")))

# # for verbal communicators only
# tmp <- inner_join(spark.rm%>%select(IID, mph_effect), 
#                   spark.scq%>%select(IID, starts_with("corr"), starts_with("q"))) %>% 
#   filter(q01_phrases==1) %>% 
#   select(-c(q01_phrases))
scq.to.mph.eff <- do.call(rbind,
                          lapply(tmp %>% select(starts_with("q")), function(x) {
                            df <- fisher_table(tmp$mph_effect, x)
                            return(df)})) %>% 
  as.data.frame() %>% 
  rownames_to_column("V2") %>%
  mutate(V1 = "mph_effect",
         sig = ifelse(pval<0.05, paste0("pval < 0.05"), paste0("pval >= 0.05")),
         V2 = sub("q.._","",V2), FDR = p.adjust(pval,"fdr"))
scq.to.mph.eff %>% mutate(CI = paste0("(",signif(conf_min,3)," - ", signif(conf_max,3), ")"),
                          OR=signif(OR,3), pval = signif(pval,3), FDR = signif(FDR,3)) %>% 
  select(SCQ_item = V2, OR, CI, pval, FDR) %>% arrange(by = pval) %>% 
  mutate(SCQ_item = str_replace_all(SCQ_item,"_"," ")) %>%
  write_tsv("data/derivatives/SPARK-scq-to-mph-eff.tsv")

library(showtext)
p1 <- scq.to.mph.eff %>%
  ggplot(aes(x=OR, y = reorder(V2,desc(OR)))) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.5),labels = scales::label_parse()) +
  scale_shape_manual(values = c(1, 2)) + 
  geom_errorbarh(aes(xmin = conf_min, xmax = conf_max, alpha = sig), 
                 size = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  bw.theme +
  labs(x = "OR (methylphenidate effectiveness and SCQ items)", y="", alpha="",
       caption = paste0("n(samples): ", nrow(tmp), "\n", "odds ratio from Fisher's Exact test"))
p1
ggsave(p1,filename = "../2026/figs/spark-scq-to-mph-effectiveness-RM.pdf",width = 6,height = 14)
p2 <- corr.table(tmp%>%select(mph_effect,starts_with("q"))) %>%
  filter(V1 == "mph_effect", V2 !=V1) %>%
  ggplot(aes(x=V1, y=V2, fill=r, label=ifelse(pval<0.05, paste0("r: ", round(r, 4), ",  p: ", round(pval, 4)),"")))+
  geom_tile() +
  geom_text(size=3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  my.guides + labs(x="", y="", caption = paste0("n(samples): ", nrow(tmp))) +
  theme(axis.text.x.bottom = element_text(hjust = 0.5, angle = 0))
p2
# conclusion: some sig correlations were found
################################################################################
################################################################################
# corr w cog impairment
tmp2 <- inner_join(spark.rm%>%select(IID, mph_effect), spark.iq)
t<- fisher.test(tmp2$mph_effect, tmp2$derived_cog_impair)
p3 <- tmp2 %>% mutate(mph_effect = ifelse(mph_effect==1,"yes","no"), derived_cog_impair=ifelse(derived_cog_impair==1,T,F)) %>%
  ggplot(aes(mph_effect, fill=derived_cog_impair))+
  geom_bar() +
  scale_fill_manual(values = palette.1)+
  labs(subtitle = paste0("Fisher's Exact test", "\n",
                         "OR: ", round(t$estimate[[1]],4), 
                        "\npvalue: ", round(t$p.value,4))) +
  bw.theme+labs(x="effective methylphenidate",fill="cognitive impairment")
p3
ggsave2("../2026/figs/spark-cognitive-impairment-to-mph-effectiveness-RM.pdf",4,5)
# conclusion: sig correlation
################################################################################
################################################################################
# corr w child background data
tmp3 <- inner_join(spark.rm%>%select(IID, mph_effect), spark.bg.child)
m <- lm(mph_effect ~ ., data = tmp3%>%select(-IID))
summary(m)$coefficients %>% as.data.frame() %>%
  rownames_to_column("var") %>% filter(var != "(Intercept)") %>%
  mutate(confint_min = confint(m)[-1,1], confint_max = confint(m)[-1,2]) %>%
  mutate(sig = ifelse(`Pr(>|t|)`<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y = var), color = var) +
           geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5) +
           geom_vline(xintercept = 0, linetype = "dashed", size = 0.2, color = "red") +
           scale_alpha_manual(values = c(1, 0.5)) +
           scale_shape_manual(values = c(1, 2)) + 
           geom_errorbarh(aes(xmin = confint_min, xmax = confint_max,alpha = sig), 
                          size = 0.8, height = 0, show.legend = F, 
                          position = position_dodge(width = 0.6)) +
           # scale_color_manual(values = six.colors[1:4]) +
           my.guides+labs(x ="Estimate", y="", 
                          title = "predicting MPH_effect by cog/function/language age level", 
                          caption = paste0("lm(mph_effect ~ cog_age_level + langcog_age_level + functioncog_age_level)")) 
# conclusion: no sig correlation for any
################################################################################
################################################################################
# corr w predicted mph response
spark.tt.pred <- read_rds("../data/derivatives/m-outputs/spark/mph-samples/model-celltype-all-FALSE-TRUE-1.rds")
tmp4 <- inner_join(spark.rm%>%select(IID, mph_effect), spark.tt.pred %>% mutate(m=scale(-m)))
tmp4 %>% mutate(mph_effect = case_when(mph_effect==1~"yes",mph_effect==0~"no"))%>%
  ggplot(aes(mph_effect,m,fill=mph_effect)) +
  geom_violin(show.legend = F)+geom_boxplot(fill="white",width=0.2)+
  ggpubr::stat_compare_means(color="red",label.x.npc = 0.25)+
  scale_fill_manual(values=palette.1)+bw.theme+
  labs(x="effective methylphenidate",y="predicted response to methylphenidate")
ggsave2("../2026/figs/spark-mph-effect-w-predicted.pdf",4,5)
cor.test(tmp4$mph_effect, tmp4$m)
# conclusion: no sig correlation
################################################################################
################################################################################
################################################################################
