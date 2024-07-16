################################################################################
#                                 pre-res lookups                              #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response"
setwd(project.dir)
################################################################################
spark.dir <- "/Dedicated/jmichaelson-sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v9_2022-12-12"
################################################################################
# read the SCQ data, and correct the yes/no values
spark.scq.raw <- read_csv(paste0(spark.dir, "/scq_2022-12-12.csv"))
# scq.flip.q <- c(paste0("q0", 9), paste0("q", 19:40))
# scq.q.dict <- data.frame(question_n = c("q09", paste0("q", c(10,19,26,27,28,29,30,31,32,33,36,37,39,40)), # reciprocal soc intera
#                                         paste0("q0", c(2,3,4,5,6)), paste0("q", c(20,21,22,23,24,25,34)), #communication
#                                         paste0("q0", c(7,8)), paste0("q", c(11,12,13,14,15,16)), # restrictive, repetitive stereotype behav
#                                         paste0("q", c(38,17,18)) # NA
#                                         ),
#                          subscale = c(rep("reciprocal_social_interaction", 15),
#                                       rep("communication", 12),
#                                       rep("restricted_repetitive_stereotyped_behav", 8),
#                                       rep("NOT", 3))) %>%
#   mutate(question_name = colnames(spark.scq.raw %>% select(starts_with("q")))[readr::parse_number(question_n)]) %>%
#   mutate(yn_flip = ifelse(question_n %in% scq.flip.q, T, F))
# spark.scq.flipped <- spark.scq.raw %>% filter(missing_values == 0) %>%
#   mutate_at(.vars = vars(starts_with("q")), .funs = function(x) ifelse(is.na(x), 0, x)) %>%
#   mutate_at(.vars = vars(scq.q.dict$question_name[scq.q.dict$yn_flip]), .funs = function(x) ifelse(x == 0, 1, 0)) %>%
#   rename(IID = subject_sp_id)
spark.scq <- spark.scq.raw %>% filter(missing_values == 0) %>%
  mutate_at(.vars = vars(starts_with("q")), .funs = function(x) ifelse(is.na(x), 0, x)) %>%
  rename(IID = subject_sp_id)
################################################################################
# read the child background data
spark.bg.child.raw <- read_csv(paste0(spark.dir, "/background_history_child_2022-12-12.csv"), na = c("888"))
spark.bg.child <- spark.bg.child.raw %>% select(IID = subject_sp_id, contains("age_level")) %>%
  drop_na()
################################################################################
# read the cognitive data 
spark.iq.raw <- read_csv(paste0(spark.dir, "/predicted_iq_experimental_2022-12-12.csv"))
spark.iq <- spark.iq.raw %>% select(IID = subject_sp_id, derived_cog_impair) %>%
  drop_na()
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
p1 <- do.call(rbind,
        lapply(tmp %>% select(starts_with("q")), function(x) {
          test <- fisher.test(tmp$mph_effect, x)
          df <- data.frame(pval = test$p.value,
                           OR = test$estimate[[1]],
                           confint_min = test$conf.int[1],
                           confint_max = test$conf.int[2])
          return(df)})) %>% 
  as.data.frame() %>% 
  rownames_to_column("V2") %>%
  mutate(V1 = "mph_effect") %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  mutate(V2 = as.factor(V2)) %>%
  ggplot(aes(x=OR, y = V2)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_shape_manual(values = c(1, 2)) + 
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 size = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = six.colors[1:4]) +
  labs(x = "", y="", caption = paste0("n(samples): ", nrow(tmp), "\n", 
                                      "odds ratio from Fisher's Exact test"))
p2 <- corr.table(tmp%>%select(mph_effect), tmp%>%select(starts_with("q"))) %>%
  filter(V1 == "mph_effect", V2 !=V1) %>%
  ggplot(aes(x=V1, y=V2, fill=r, label=ifelse(pval<0.05, paste0("r: ", round(r, 4), ",  p: ", round(pval, 4)),"")))+
  geom_tile() +
  geom_text(size=3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  my.guides + labs(x="", y="", caption = paste0("n(samples): ", nrow(tmp))) +
  theme(axis.text.x.bottom = element_text(hjust = 0.5, angle = 0))
# conclusion: some sig correlations were found
################################################################################
# corr w cog impairment
tmp2 <- inner_join(spark.rm%>%select(IID, mph_effect), spark.iq)
t<- fisher.test(tmp2$mph_effect, tmp2$derived_cog_impair)
p3 <- tmp2 %>% 
  mutate(mph_effect = as.factor(mph_effect), derived_cog_impair = as.factor(derived_cog_impair)) %>%
  ggplot(aes(mph_effect, fill=derived_cog_impair))+
  geom_bar() +
  scale_fill_manual(values = boxplot.colors)+
  labs(subtitle = paste0("Fisher's Exact test", "\n",
                         "OR: ", round(t$estimate[[1]],4), 
                        "\npvalue: ", round(t$p.value,4))) +
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5))
# conclusion: sig correlation
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
# corr w predicted mph response
spark.tt.pred <- read_rds("../data/derivatives/m-outputs/spark/mph-samples/model-celltype-all-FALSE-TRUE-1.rds")
tmp4 <- inner_join(spark.rm%>%select(IID, mph_effect), spark.tt.pred)
cor.test(tmp4$mph_effect, tmp4$m)
# conclusion: no sig correlation
################################################################################
# corr w pgs
spark.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv")
tmp5 <- inner_join(spark.rm%>%select(IID, mph_effect), spark.pgs)
p4 <- corr.table(tmp5%>%select(mph_effect), scale(tmp5%>%select(starts_with("corrected")), scale = T, center = T), method = "spearman") %>%
  filter(V1 == "mph_effect", V2 !=V1) %>%
  mutate(V2 = sub("corrected_", "", V2)) %>%
  # filter(V2 %in% c("ADHD-Demontis", "cog_memory-UKB-2020", "cog_reaction_time-UKB-2020",
  #                  "cog_gFactor-UKB-2020", "cog_matrix-UKB-2020", "cog_symbol_digit-UKB-2020",
  #                  "cog_tower_rearranging-UKB-2020", "cog_trail_making_testB-UKB-2020",
  #                  "cog_verbal_numerical_reasoning-UKB-2020", "EA-Lee", "CP-Lee")) %>%
  mutate(V2 = sub("cog_", "", sub("-UKB-2020", "", sub("-Lee", "", sub("-Demontis", "", V2))))) %>%
  rename(corr = r) %>%
  ggplot(aes(x=V1, y=V2, fill=corr, label=ifelse(pval<0.1, paste0("ρ: ", round(corr, 4), ",  p: ", round(pval, 4)),"")))+
  geom_tile() +
  geom_text(size=3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  my.guides+labs(x="", y="", fill = "ρ", caption = paste0("n(samples): ", nrow(tmp))) +
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5))
# conclusion: sig correlation was found for a few number of PGS
################################################################################
################################################################################
# combine
patchwork::wrap_plots(p1,patchwork::wrap_plots(p4,p3, nrow = 2))
ggsave("figs/0724_report/SPARK-mph-eff-results.png", bg = "white",
       width = 9, height = 12, units = "in", dpi = 360)
################################################################################
################################################################################
# supplementary figs
# SCQ dist
tmp %>% select(starts_with("q"), mph_effect) %>%
  pivot_longer(cols = c(starts_with("q"), mph_effect), names_to = "q", values_to = "v") %>%
  mutate(v = as.factor(v))%>%
  ggplot(aes(x=v, fill = v))+
  geom_bar(show.legend = F)+
  facet_wrap("q", ncol = 5)+
  scale_fill_manual(values = six.colors[c(1,3)])+
  labs(x="", caption = paste0("n(samples): ", nrow(tmp)))
# pgs dist
tmp5 %>%
  select(starts_with("corrected")) %>%
  pivot_longer(cols = starts_with("corrected"), names_to = "pgs", values_to = "v") %>%
  filter(grepl("ADHD", pgs) | grepl("cog", pgs) | grepl("EA_", pgs)) %>%
  mutate(pgs = sub("corrected_", "", pgs)) %>%
  ggplot(aes(x=v))+
  geom_histogram()+
  facet_wrap("pgs", scales = "free", ncol = 4) + labs(x="", caption = paste0("n(samples): ", nrow(tmp5)))
################################################################################