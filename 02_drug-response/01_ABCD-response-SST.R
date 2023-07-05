################################################################################
#                   predcit drug response for ABCD participants                #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response"
setwd(project.dir)
################################################################################
################################################################################
################################################################################
# dimensionality reduction for task features of SST
abcd.raw.dir <- "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_4_0"
sst.raw <- read_delim(paste0(abcd.raw.dir, "/abcd_sst02.txt"), delim = "\t")
question.dict <- data.frame(q0 = colnames(sst.raw), description = t(sst.raw)[,1], row.names = NULL)
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meds-matrix_CMAP-meds.rds") %>%
  as.data.frame() %>%
  select(c(1:6), methylphenidate) %>%
  rename(IID = ID)

#########################

# decided to keep count of answers in trial run 1 (r1) only
sst.r1 <- left_join(sst.raw[-1,] %>%
                      select(IID = subjectkey, interview_age, sex, eventname,
                             ends_with("_nt")&contains("r1"), 
                             tfmri_sst_beh_switchflag,tfmri_sst_beh_performflag,
                             -c(contains("_beh_go_nt"), contains("_beh_s_nt"))) %>%
                      filter(is.na(tfmri_sst_beh_switchflag) == F),
                    abcd.meds) %>% drop_na() %>%
  mutate_at(.vars = vars(ends_with("_nt"), ends_with("flag"), ends_with("age")), .funs = function(x) as.numeric(x))
hist(sst.r1)

#########################

# check if there's age correlation with these questions
corr.table(sst.r1%>%select(interview_age), sst.r1%>%select(ends_with("nt"))) %>%
  filter(V1 == "interview_age", V2 != V1) %>%
  ggplot(aes(x=V1, y = V2, fill = r, label = paste0("r: ", round(r, 3), ",  p: ", round(pval, 5)))) +
  geom_tile() + 
  geom_text()+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  my.guides+labs(x = "", y="", caption = paste0("n(samples): ", nrow(sst.r1)))
# age is found to be significantly affecting all, so correct for it

#########################

# do archetypal analysis for reducing number of dimensions after correcting for age effect
library(archetypes)
x.r1 <- apply(sst.r1%>%select(ends_with("_nt")), 
              MARGIN = 2, 
              FUN = function(x) residuals(glm(x ~ sst.r1$interview_age, family = quasipoisson())))
sst.a.r1 <- stepArchetypes(x.r1, 
                           k = 1:10)
#########################
# check the RSS and decide on what number of archetypes to keep
screeplot(sst.a.r1)
# kept 4 only
a4.r1 <- bestModel(sst.a.r1[[4]])
# check the archetypes loadings
a4.r1$archetypes %>% as.data.frame() %>%
  mutate(A = paste0("A", 1:4)) %>%
  pivot_longer(cols = -c("A"), names_to = "question", values_to = "loading") %>%
  ggplot(aes(x=A, y = question, fill = loading, label = ifelse(abs(loading)>2, round(loading, 2), "")))+
  geom_tile()+
  geom_text(size=3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1]) +
  my.guides+labs(y="", x="", caption = paste0("n(samples): ", nrow(sst.r1), "\n",
                                              "questions with missing values were dropped", "\n", 
                                              "question val = residuals(glm(x ~ sst.r1$interview_age, family = quasipoisson()))"))

# save the alphas 
sst.archetypes.r1 <- cbind(a4.r1$alphas, sst.r1%>%select(1:4, methylphenidate))
colnames(sst.archetypes.r1)[1:4] <- paste0("A", 1:4, "_r1_nt")
write_tsv(sst.archetypes.r1, "data/derivatives/abcd-sst-archetypes-r1-nt.tsv")
################################################################################
# correlation between abcd pgs and sst archetypes & abcd predicted mph and sst archetypes
# indep section
sst.a <- read_tsv("data/derivatives/abcd-sst-archetypes-r1-nt.tsv")
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv")
colnames(abcd.pgs) <- str_replace_all(pattern = "-",replacement =  "_",string =  colnames(abcd.pgs))
abcd.pred <- read_rds("../data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(predicted = m) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
abcd.c <- inner_join(sst.a, abcd.pgs)
abcd.c <- left_join(abcd.c, abcd.pred)

#########################
# sst archetypes w predicted response to MPH
library(betareg)
m1 <- betareg(y ~ methylphenidate + predicted + methylphenidate:predicted, data = abcd.c %>% mutate(y = A1_r1_nt+1e-12))
m2 <- betareg(y ~ methylphenidate + predicted + methylphenidate:predicted, data = abcd.c %>% mutate(y = A2_r1_nt+1e-12))
m3 <- betareg(y ~ methylphenidate + predicted + methylphenidate:predicted, data = abcd.c %>% mutate(y = A3_r1_nt+1e-12))
m4 <- betareg(y ~ methylphenidate + predicted + methylphenidate:predicted, data = abcd.c %>% mutate(y = A4_r1_nt+1e-12))
rbind(
  summary(m1)$coefficients %>% 
    as.data.frame() %>% rownames_to_column("var") %>% mutate(archetype = "A1_r1_nt") %>%
    mutate(confint_min = confint(m1)[1:4,1],confint_max = confint(m1)[1:4,2]),
  summary(m2)$coefficients %>% 
    as.data.frame() %>% rownames_to_column("var") %>% mutate(archetype = "A2_r1_nt") %>%
    mutate(confint_min = confint(m2)[1:4,1],confint_max = confint(m2)[1:4,2]),
  summary(m3)$coefficients %>% 
    as.data.frame() %>% rownames_to_column("var") %>% mutate(archetype = "A3_r1_nt") %>%
    mutate(confint_min = confint(m3)[1:4,1],confint_max = confint(m3)[1:4,2]),
  summary(m4)$coefficients %>% 
    as.data.frame() %>% rownames_to_column("var") %>% mutate(archetype = "A4_r1_nt") %>%
    mutate(confint_min = confint(m4)[1:4,1],confint_max = confint(m4)[1:4,2])
) %>%
  mutate(sig = ifelse(`mean.Pr...z..`<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  filter(var != "(Intercept)") %>%
  mutate(var = factor(var, levels = c("predicted", "methylphenidate", "methylphenidate:predicted"))) %>%
  mutate(archetype = as.factor(archetype)) %>%
  ggplot(aes(x=`mean.Estimate`, y = var), color = var) +
  geom_point(aes(color = var, alpha = sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_shape_manual(values = c(1, 2)) + 
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, color = var, alpha = sig), 
                 size = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = six.colors[1:4]) +
  facet_wrap("archetype", scales = "free_x", nrow = 1) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  my.guides+labs(x ="Estimate", y="", 
                 title = "predicting SST_r1_nt archetypess by MPH predicted response", 
                 caption = paste0("models: betareg(y ~ methylphenidate + predicted + methylphenidate:predicted, data = abcd.c %>% mutate(y = A[1:4]+1e-12))", "\n",
                                  "predicted values are scaled and the A[1:4] are the raw alphas from the archetypes function", "\n", 
                                  "n(samples):", nrow(abcd.c), "\n",
                                  "n(samples taking MPH): ", nrow(abcd.c%>%filter(methylphenidate==1)), "\n",
                                  "n(unique samples):", nrow(abcd.c%>%distinct(IID, .keep_all = T)), "\n",
                                  "n(unique samples taking MPH): ", nrow(abcd.c%>%filter(methylphenidate==1)%>%distinct(IID, .keep_all = T))))

#########################
# sst archetypes w pgs 
# corr.table(abcd.c%>% select(starts_with("A")), scale(abcd.c%>%select(starts_with("corrected"))), method = "spearman") %>%
#   filter(V1 %in% c(paste0("A", 1:4, "_r1_nt")), ! V2 %in% c(paste0("A", 1:4, "_r1_nt"))) %>%
#   mutate(V2 = sub("corrected_", "", V2)) %>%
#   filter(V2 %in% c("ADHD-Demontis", "cog_memory-UKB-2020", "cog_reaction_time-UKB-2020",
#                    "cog_gFactor-UKB-2020", "cog_matrix-UKB-2020", "cog_symbol_digit-UKB-2020",
#                    "cog_tower_rearranging-UKB-2020", "cog_trail_making_testB-UKB-2020",
#                    "cog_verbal_numerical_reasoning-UKB-2020", "EA-Lee")) %>%
#   mutate(V2 = sub("cog_", "", sub("-UKB-2020", "", sub("-Lee", "", sub("-Demontis", "", V2))))) %>%
#   rename(corr = r) %>%
#   ggplot(aes(x=V1, y=V2, fill=corr, label=ifelse(pval<0.05, paste0("corr: ", round(corr, 3), ",  p: ", round(pval, 4)), "")))+
#   geom_tile() +
#   scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
#   geom_text(size=3)+
#   my.guides + labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c), "\n",
#                                                 "Archetypes weights are the raw alphas from the archetypes function here", "\n", 
#                                                 "the pgs were z-scaled for the selected sample", "\n", 
#                                                 "shown corr values are spearman correlation"))
# maybe you should get these as a betareg model with including taking mph and interaction
library(betareg)
A1.pgs <- do.call(rbind, lapply(abcd.c%>%select(starts_with("corrected")), function(x) {
  m = betareg(y ~ methylphenidate + pgs + methylphenidate:pgs, 
          data = abcd.c %>% mutate(y = A1_r1_nt+1e-12, pgs = scale(x, scale = T, center = T)))
  summary(m)$coefficients %>% 
    as.data.frame() %>% rownames_to_column("var") %>% mutate(archetype = "A1_r1_nt") %>%
    mutate(confint_min = confint(m)[1:4,1],confint_max = confint(m)[1:4,2])
})) %>% as.data.frame() %>%
  rownames_to_column("V2") %>%
  mutate(V2 = as.factor(sub("\\.[0-9]", "", sub("corrected_", "", V2))))
A2.pgs <- do.call(rbind, lapply(abcd.c%>%select(starts_with("corrected")), function(x) {
  m = betareg(y ~ methylphenidate + pgs + methylphenidate:pgs, 
              data = abcd.c %>% mutate(y = A2_r1_nt+1e-12, pgs = scale(x, scale = T, center = T)))
  summary(m)$coefficients %>% 
    as.data.frame() %>% rownames_to_column("var") %>% mutate(archetype = "A2_r1_nt") %>%
    mutate(confint_min = confint(m)[1:4,1],confint_max = confint(m)[1:4,2])
})) %>% as.data.frame() %>%
  rownames_to_column("V2") %>%
  mutate(V2 = as.factor(sub("\\.[0-9]", "", sub("corrected_", "", V2))))
A3.pgs <- do.call(rbind, lapply(abcd.c%>%select(starts_with("corrected")), function(x) {
  m = betareg(y ~ methylphenidate + pgs + methylphenidate:pgs, 
              data = abcd.c %>% mutate(y = A3_r1_nt+1e-12, pgs = scale(x, scale = T, center = T)))
  summary(m)$coefficients %>% 
    as.data.frame() %>% rownames_to_column("var") %>% mutate(archetype = "A3_r1_nt") %>%
    mutate(confint_min = confint(m)[1:4,1],confint_max = confint(m)[1:4,2])
})) %>% as.data.frame() %>%
  rownames_to_column("V2") %>%
  mutate(V2 = as.factor(sub("\\.[0-9]", "", sub("corrected_", "", V2))))
A4.pgs <- do.call(rbind, lapply(abcd.c%>%select(starts_with("corrected")), function(x) {
  m = betareg(y ~ methylphenidate + pgs + methylphenidate:pgs, 
              data = abcd.c %>% mutate(y = A4_r1_nt+1e-12, pgs = scale(x, scale = T, center = T)))
  summary(m)$coefficients %>% 
    as.data.frame() %>% rownames_to_column("var") %>% mutate(archetype = "A4_r1_nt") %>%
    mutate(confint_min = confint(m)[1:4,1],confint_max = confint(m)[1:4,2])
})) %>% as.data.frame() %>%
  rownames_to_column("V2") %>%
  mutate(V2 = as.factor(sub("\\.[0-9]", "", sub("corrected_", "", V2))))
rbind(A1.pgs, A2.pgs, A3.pgs, A4.pgs) %>% 
  mutate(sig = ifelse(`mean.Pr...z..`<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  filter(var != "(Intercept)") %>%
  filter(V2 %in% c("ADHD_Demontis", "cog_memory_UKB_2020", "cog_reaction_time_UKB_2020",
                   "cog_gFactor-UKB_2020", "cog_matrix_UKB_2020", "cog_symbol_digit_UKB_2020",
                   "cog_tower_rearranging_UKB_2020", "cog_trail_making_testB_UKB_2020",
                   "cog_verbal_numerical_reasoning_UKB_2020", "EA_Lee")) %>%
  mutate(V2 = sub("cog_", "", sub("_UKB_2020", "", sub("_Lee", "", sub("_Demontis", "", V2))))) %>%
  mutate(var = factor(var, levels = c("pgs", "methylphenidate", "methylphenidate:pgs"))) %>%
  mutate(archetype = as.factor(archetype)) %>%
  ggplot(aes(x=`mean.Estimate`, y = V2), color = var) +
  geom_point(aes(color = var, alpha = sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_shape_manual(values = c(1, 2)) + 
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, color = var, alpha = sig), 
                 size = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = six.colors[1:4]) +
  facet_wrap("archetype", scales = "free_x", nrow = 1) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  my.guides+labs(x ="Estimate", y="", 
                 title = "predicting SST_r1_nt archetypess by MPH predicted response", 
                 caption = paste0("models: betareg(y ~ methylphenidate + pgs + methylphenidate:pgs, data = abcd.c %>% mutate(y = A[1:4]+1e-12))", "\n",
                                  "pgs are z-scaled for the selected samples and the A[1:4] are the raw alphas from the archetypes function", "\n", 
                                  "n(samples):", nrow(abcd.c), "\n",
                                  "n(samples taking MPH): ", nrow(abcd.c%>%filter(methylphenidate==1)), "\n",
                                  "n(unique samples):", nrow(abcd.c%>%distinct(IID, .keep_all = T)), "\n",
                                  "n(unique samples taking MPH): ", nrow(abcd.c%>%filter(methylphenidate==1)%>%distinct(IID, .keep_all = T))))

################################################################################
# sst archetypes w nih toolbox metrics
tbx.raw <- read_delim(paste0(abcd.raw.dir, "/abcd_tbss01.txt"), delim = "\t")
tbx.filt <- tbx.raw[-1,] %>%
            select(IID = subjectkey, interview_age, sex, eventname,
                   ends_with("_agecorrected")) %>% 
  # drop_na() %>%
  mutate_at(.vars = vars(ends_with("_agecorrected"), ends_with("age")), .funs = function(x) as.numeric(x)) %>%
  mutate(NA_tot = rowSums(is.na(tbx.raw[-1,]%>%select(ends_with("_agecorrected"))))) %>%
  filter(NA_tot != 10)
tmp <- inner_join(sst.a, tbx.filt)
corr.table(tmp%>%select(starts_with("A")), tmp%>%select(ends_with("_agecorrected")), method = "spearman") %>%
  filter(V1 %in% paste0("A", 1:4, "_r1_nt"), !V2 %in% paste0("A", 1:4, "_r1_nt")) %>%
  mutate(V2 = sub("nihtbx_", "", V2)) %>%
  mutate(V2 = sub("_agecorrected", "", V2)) %>%
  rename(corr = r) %>%
  ggplot(aes(x=V1, y=V2, fill = corr, label = ifelse(pval<0.05, paste0("ρ: ", round(corr, 3), "\n", "pval: ", round(pval, 10)), "")))+
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(high = redblu.col[1], low = redblu.col[2])+
  my.guides+labs(x="", y="", fill = "ρ", caption = paste0("n(samples): ", nrow(tmp), "\n", 
                                              "the NIH toolbox values are all age-corrected and were retrieved from the raw file", "\n",
                                              "the reported values are Spearman rank-based correlation coefficients (ρ)"))

#########################
# NIH tbx w predicted response to MPH
tmp2 <- inner_join(abcd.pred, tmp)
do.call(rbind, apply(tmp2 %>% select(ends_with("corrected")), MARGIN = 2, 
      function(x) {m = glm(y ~ predicted + methylphenidate + predicted:methylphenidate, 
                           data = tmp2 %>% mutate(y = x), family = quasipoisson())
      summary(m)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("var") %>%
        mutate(confint_min = confint(m)[,1], confint_max = confint(m)[,2])})) %>%
  rownames_to_column("V2") %>%
  mutate(V2 = sub("nihtbx_", "", V2)) %>%
  mutate(V2 = sub("_agecorr.*", "", V2)) %>%
  filter(var != "(Intercept)") %>%
  mutate(sig = ifelse(`Pr(>|t|)`<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  mutate(var = factor(var, levels = c("predicted", "methylphenidate", "predicted:methylphenidate"))) %>%
  mutate(V2 = as.factor(V2)) %>%
  ggplot(aes(x=Estimate, y = var), color = var) +
  geom_point(aes(color = var, alpha = sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_shape_manual(values = c(1, 2)) + 
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, color = var, alpha = sig), 
                 size = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = six.colors[1:3]) +
  facet_wrap("V2", scales = "free_x", nrow = 2) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  my.guides+labs(x ="Estimate", y="", 
                 title = "predicting NIH tbx data by MPH predicted response", 
                 caption = paste0("models: glm(y ~ predicted + methylphenidate + predicted:methylphenidate, family = quasipoisson())", "\n",
                                  "predicted values are z-scaled and the y values are the raw age-corrected values from the NIH TBX", "\n", 
                                  "n(samples):", nrow(tmp2), "\n",
                                  "n(samples taking MPH): ", nrow(tmp2%>%filter(methylphenidate==1)), "\n",
                                  "n(unique samples):", nrow(tmp2%>%distinct(IID, .keep_all = T)), "\n",
                                  "n(unique samples taking MPH): ", nrow(tmp2%>%filter(methylphenidate==1)%>%distinct(IID, .keep_all = T))))

#########################
# cbcl w predicted response to MPH
abcd.cbcl <- read_delim(paste0(abcd.raw.dir, "/abcd_cbcl01.txt"), delim = "\t")
cbcl.scales <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/cbcl-scales.rds")
abcd.cbcl.imp <- abcd.cbcl[-1,] %>%
  mutate_at(.vars = vars(starts_with("cbcl_q"), ends_with("age")), function(x) as.numeric(x)) %>%
  mutate_at(.vars = vars(starts_with("cbcl_q")), function(x) ifelse(is.na(x), median(x, na.rm = T), x))
abcd.cbcl.imp$syn_attention <- rowSums(abcd.cbcl.imp%>%select(cbcl.scales$q0[cbcl.scales$cbcl_syndrome=="attention"]))
abcd.cbcl.imp$syn_rulebreaking <- rowSums(abcd.cbcl.imp%>%select(cbcl.scales$q0[cbcl.scales$cbcl_syndrome=="rulebreaking"]))
abcd.cbcl.imp$syn_aggressive <- rowSums(abcd.cbcl.imp%>%select(cbcl.scales$q0[cbcl.scales$cbcl_syndrome=="aggressive"]))
abcd.cbcl.imp$syn_withdrawn <- rowSums(abcd.cbcl.imp%>%select(cbcl.scales$q0[cbcl.scales$cbcl_syndrome=="withdrawn"]))
abcd.cbcl.imp$syn_other <- rowSums(abcd.cbcl.imp%>%select(cbcl.scales$q0[cbcl.scales$cbcl_syndrome=="other"]))
abcd.cbcl.imp$syn_thought <- rowSums(abcd.cbcl.imp%>%select(cbcl.scales$q0[cbcl.scales$cbcl_syndrome=="thought"]))
abcd.cbcl.imp$syn_social <- rowSums(abcd.cbcl.imp%>%select(cbcl.scales$q0[cbcl.scales$cbcl_syndrome=="social"]))
abcd.cbcl.imp$syn_somatic <- rowSums(abcd.cbcl.imp%>%select(cbcl.scales$q0[cbcl.scales$cbcl_syndrome=="somatic"]))

####
# check if there's age correlation between cbcl syndromes and age
corr.table(abcd.cbcl.imp %>% select(interview_age), abcd.cbcl.imp %>% select(starts_with("syn"))) %>%
  filter(V1 == "interview_age", V2 != V1) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 3)), "")))+
  geom_tile()+
  geom_text(size=3)+
  scale_fill_gradient2(high = redblu.col[1], low = redblu.col[2])+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.cbcl.imp)))
# conc: age is sig affecting the tot syndrome score

# check if there's sex correlation between cbcl syndromes and age
corr.table(abcd.cbcl.imp %>% select(sex) %>% mutate(sex = ifelse(sex == "M", 0, 1)), abcd.cbcl.imp %>% select(starts_with("syn"))) %>%
  filter(V1 == "sex", V2 != V1) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 3)), "")))+
  geom_tile()+
  geom_text(size=3)+
  scale_fill_gradient2(high = redblu.col[1], low = redblu.col[2])+
  my.guides+labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.cbcl.imp)))
# conc: sex is sig affecting the tot syndrome score

hist(abcd.cbcl.imp %>% select(starts_with("syn"))) # quasi-poisson dist
# you either decide to correct for age, sex, and their interactions first or do that in the model at the end
cbcl.corrected <- cbind(abcd.cbcl.imp %>% select(IID = subjectkey, starts_with("interview"), sex, eventname), apply(abcd.cbcl.imp%>%select(starts_with("syn")), MARGIN = 2, 
                        function(x) residuals(glm(y ~ interview_age + sex + interview_age:sex, data = abcd.cbcl.imp %>% mutate(y = x), family = quasipoisson()))))

tmp3 <- inner_join(abcd.pred, inner_join(abcd.meds%>%select(IID, eventname, methylphenidate), abcd.cbcl.imp %>% select(IID = subjectkey, starts_with("interview"), sex, eventname, starts_with("syn"))))
hist(tmp3)
coeff <- do.call(rbind, apply(tmp3 %>% select(starts_with("syn")), MARGIN = 2, function(x) {
  m <- glm(y ~ methylphenidate + predicted + predicted:methylphenidate + interview_age + sex + interview_age:sex, data = tmp3 %>% mutate(y = x, interview_age = scale(interview_age)[,1]), family = quasipoisson())
  summary(m)$coefficients %>% as.data.frame() %>%
    rownames_to_column("var") %>%
    mutate(confint_min = confint(m)[,1], confint_max = confint(m)[,2])
  }))
coeff %>% 
  as.data.frame() %>%
  rownames_to_column("V2") %>%
  mutate(V2 = as.factor(sub("\\.[1-9]*", "", V2))) %>%
  filter(var != "(Intercept)") %>%
  mutate(sig = ifelse(`Pr(>|t|)`<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  mutate(var = factor(var, levels = c("predicted", "methylphenidate", "methylphenidate:predicted", "interview_age", "sexM", "interview_age:sexM"))) %>%
  ggplot(aes(x=Estimate, y = var), color = var) +
  geom_point(aes(color = var, alpha = sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_shape_manual(values = c(1, 2)) + 
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, color = var, alpha = sig), 
                 size = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = six.colors[1:6]) +
  facet_wrap("V2", scales = "free_x", nrow = 3) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  my.guides+labs(x ="Estimate", y="", 
                 title = "predicting CBCL syndromes by MPH predicted response", 
                 caption = paste0("models: glm(y ~ methylphenidate + predicted + predicted:methylphenidate + interview_age + sex + interview_age:sex, family = quasipoisson())", "\n",
                                  "predicted values are z-scaled and the y values are the raw syndrome scores after imputing missing values by the median per question", "\n",
                                  "the interview_age variable was z-scaled before being added to the model", "\n", 
                                  "n(samples):", nrow(tmp3), "\n",
                                  "n(samples taking MPH): ", nrow(tmp3%>%filter(methylphenidate==1)), "\n",
                                  "n(unique samples):", nrow(tmp3%>%distinct(IID, .keep_all = T)), "\n",
                                  "n(unique samples taking MPH): ", nrow(tmp3%>%filter(methylphenidate==1)%>%distinct(IID, .keep_all = T))))
################################################################################
################################################################################
################################################################################
# supplementary figs
# distribution of SST archetypes, pgs, age, predicted
p1 <- abcd.c %>% select(1:4) %>%
  pivot_longer(cols = starts_with("A"), names_to = "archetype", values_to = "alpha") %>%
  ggplot(aes(x=alpha))+
  geom_histogram()+
  facet_wrap("archetype", scales = "free_y", nrow = 1)
p2 <- abcd.c %>%
  ggplot(aes(x=interview_age))+
  geom_histogram()
p3 <- abcd.c %>%
  ggplot(aes(x=predicted))+
  geom_histogram()
p4 <- abcd.c %>%
  select(starts_with("corrected")) %>%
  pivot_longer(cols = starts_with("corrected"), names_to = "pgs", values_to = "v") %>%
  filter(grepl("ADHD", pgs) | grepl("cog", pgs) | grepl("EA_", pgs)) %>%
  mutate(pgs = sub("corrected_", "", pgs)) %>%
  ggplot(aes(x=v))+
  geom_histogram()+
  facet_wrap("pgs", scales = "free", ncol = 4) + labs(x="", caption = paste0("n(samples): ", nrow(abcd.c)))
pa <- patchwork::wrap_plots(p1,patchwork::wrap_plots(p2,p3, nrow = 1), nrow = 2)
patchwork::wrap_plots(pa, p4, ncol = 1, heights = c(3,4))

# distribution of nihtbx data, age, 
p5 <- tmp2 %>% pivot_longer(cols = starts_with("nih"), names_to = "task", values_to = "v") %>%
  ggplot(aes(x=v))+
  geom_histogram()+
  facet_wrap("task", scales = "free")+
  labs(x="", caption = paste0("n(samples): ", nrow(tmp2)))
p6 <- tmp2 %>%
  ggplot(aes(x=interview_age))+
  geom_histogram()
p7 <- tmp2 %>%
  ggplot(aes(x=predicted))+
  geom_histogram()
patchwork::wrap_plots(patchwork::wrap_plots(p6,p7, nrow = 1),p5, nrow = 2, heights = c(1,3))
################################################################################
