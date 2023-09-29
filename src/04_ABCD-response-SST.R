################################################################################
#                   predcit drug response for ABCD participants                #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response"
setwd(project.dir)
################################################################################
raw.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meds-matrix_CMAP-meds.rds") %>%
  as.data.frame() %>%
  filter(methylphenidate == 1) %>%
  rename(IID = ID)
tmp <- raw.meds[,7:ncol(raw.meds)]
raw.meds <- cbind(raw.meds[,1:6], tmp[,colSums(tmp)>0]) %>%
  mutate(tot_meds = rowSums(tmp))
rm(tmp)
gc()
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meds-matrix_CMAP-meds.rds") %>%
  as.data.frame() %>%
  select(c(1:6), methylphenidate) %>%
  # filter(methylphenidate == 1) %>%
  rename(IID = ID)
################################################################################
################################################################################
# combine models output 
models <- data.frame(file = list.files("data/derivatives/m-outputs/abcd/", pattern = "model")) %>% 
  mutate(model_d = sub("model-", "", sub(".rds", "", file))) %>%
  mutate(model = paste0("model_", c(1:length(file)))) %>%
  mutate(weights_source = sub("-.*", "", model_d))
m.out <- list()
for (i in 1:nrow(models)) {
  # i=1
  m.o <- read_rds(paste0("data/derivatives/m-outputs/abcd/", models$file[i])) %>%
    column_to_rownames("IID")
  colnames(m.o)[1] <- models$model[i]
  m.out[[i]] <- m.o
}
m.out <- do.call(cbind, m.out) %>%
  rownames_to_column("IID")

to.save <- list(m.out = m.out, models = models)
save(to.save, file = "data/derivatives/m-outputs/abcd/combined-output.RData")
################################################################################
################################################################################
################################################################################
################################################################################
######################## correlate response to SST fMRI ########################
################################################################################
################################################################################
################################################################################
abcd.raw.dir <- "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_4_0"
sst.raw <- read_delim(paste0(abcd.raw.dir, "/abcd_sst02.txt"), delim = "\t")
question.dict <- data.frame(q0 = colnames(sst.raw), description = t(sst.raw)[,1], row.names = NULL)

filtered.q.of.int <- c("tfmri_sst_all_beh_nrgo_rt","tfmri_sst_r1_beh_nrgo_rt","tfmri_sst_r2_beh_nrgo_rt",
                       "tfmri_sst_all_beh_crgo_rt","tfmri_sst_r1_beh_crgo_rt","tfmri_sst_r2_beh_crgo_rt",
                       "tfmri_sst_all_beh_crlg_rt","tfmri_sst_r1_beh_crlg_rt","tfmri_sst_r2_beh_crlg_rt",
                       "tfmri_sst_all_beh_crlg_mrt","tfmri_sst_r1_beh_crlg_mrt","tfmri_sst_r2_beh_crlg_mrt",
                       "tfmri_sst_all_beh_incrlg_rt","tfmri_sst_r1_beh_incrlg_rt","tfmri_sst_r2_beh_incrlg_rt",
                       "tfmri_sst_beh_performflag", "tfmri_sst_beh_switchflag")


sst.filt <- sst.raw %>%
  filter(subjectkey %in% m.out$IID) %>% 
  select(IID = subjectkey, interview_age, sex, eventname,
         # c(10,11,23:25, 26,28,30, 35:37, 38,40,42, 47:49, 50,52,54, 59:61, 62,64,66, 71:73,
         #   77:79, 83:85, 86,88,90, 95:97, 98:101)
         filtered.q.of.int) %>%
  filter(is.na(tfmri_sst_beh_switchflag) == F)
sst.mph.filt <- left_join(right_join(m.out,right_join(abcd.meds, sst.filt) %>% 
                                       # filter(methylphenidate %in% c(0,1))
                                       mutate(methylphenidate = ifelse(is.na(methylphenidate),0,methylphenidate))), 
                          raw.meds %>% select(IID, eventname, tot_meds)) %>% 
  mutate(tot_meds = ifelse(is.na(tot_meds), 0, tot_meds))

# write_rds(sst.mph.filt, "data/derivatives/abcd-mph-sst.rds", compress = "gz")

# look at correlation between age, number of meds taken and task performance in general
age.effect <- do.call(rbind, lapply(sst.mph.filt[,filtered.q.of.int], function(x) summary(lm(x ~ as.numeric(sst.mph.filt$interview_age)))$coefficients%>%
                       as.data.frame()%>%rownames_to_column("var"))) %>% rownames_to_column("question") %>%
  mutate(question = sub("\\.[1-9]", "", question)) %>%
  filter(var != "(Intercept)")%>%
  mutate(age_sig = ifelse(`Pr(>|t|)`<0.05, T,F))

meds.effect <- do.call(rbind, lapply(sst.mph.filt[,filtered.q.of.int], function(x) summary(lm(x ~ as.numeric(sst.mph.filt$tot_meds)))$coefficients%>%
                                      as.data.frame()%>%rownames_to_column("var"))) %>% rownames_to_column("question") %>%
  mutate(question = sub("\\.[1-9]", "", question)) %>%
  filter(var != "(Intercept)")%>%
  mutate(meds_sig = ifelse(`Pr(>|t|)`<0.05, T,F))

# include these questions affected by any of these covariates in a table
covar <- left_join(meds.effect %>% select(question, meds_sig), age.effect %>% select(question, age_sig))
# correct for these covariates if they showed significance
# either correct for age only, tot_meds only, both, or keep it as it is if there's no significant relations
plots.ls <- list()
for (i in 1:nrow(covar)) {
  # i=1
  que <- covar$question[i]
  if (covar$meds_sig[i] == T & covar$age_sig[i] == T) {
    df <- sst.mph.filt %>% select(IID, interview_age, tot_meds, which(colnames(sst.mph.filt)==que), methylphenidate, eventname) %>%
      drop_na() %>%
      rename(q = 4) %>%
      mutate(corrected_q = residuals(lm(q ~ as.numeric(interview_age) + tot_meds))) %>% 
      mutate(question = que) %>% mutate(age_sig = T, meds_sig = T) 
    plots.ls[[i]] <- df
  } else if (covar$meds_sig[i] == T) {
    df <- sst.mph.filt %>% select(IID, interview_age, tot_meds, which(colnames(sst.mph.filt)==que), methylphenidate, eventname) %>%
      drop_na() %>%
      rename(q = 4) %>%
      mutate(corrected_q = residuals(lm(q ~ tot_meds))) %>% 
      mutate(question = que) %>% mutate(age_sig = F, meds_sig = T) 
    plots.ls[[i]] <- df
  } else if (covar$age_sig[i] == T) {
    df <- sst.mph.filt %>% select(IID, interview_age, tot_meds, which(colnames(sst.mph.filt)==que), methylphenidate, eventname) %>%
      drop_na() %>%
      rename(q = 4) %>%
      mutate(corrected_q = residuals(lm(q ~ as.numeric(interview_age)))) %>% 
      mutate(question = que) %>% mutate(age_sig = T, meds_sig = F) 
    plots.ls[[i]] <- df
  }else {
    df <- sst.mph.filt %>% select(IID, interview_age, tot_meds, which(colnames(sst.mph.filt)==que), methylphenidate, eventname) %>%
      drop_na() %>%
      rename(q = 4) %>%
      mutate(corrected_q = q) %>% 
      mutate(question = que) %>% mutate(age_sig = F, meds_sig = F) 
    plots.ls[[i]] <- df
  }
}

plots.ls <- do.call(rbind, plots.ls)
library(ggpubr)
# plot the corrected question values against the predicted mph_response by each model
# make sure to only look at samples that reported taking  MPH on the same day
pdf("figs/11_abcd-sst-mph-taking-models-performance.pdf", width = 14, height = 12)
for (i in c(1,3,5,7,9,11,13,15)) {
  # i=1
  plots.ls.2[[length(plots.ls.2)+1]] <- left_join(left_join(plots.ls, m.out%>%select(IID, paste0("model_", i))), question.dict %>% rename(question = q0)) %>% 
    filter(methylphenidate == 1) %>%
    pivot_longer(cols = starts_with("model"), names_to = "model", values_to = "predicted") %>%
    mutate(predicted = -predicted) %>%
    ggplot(aes(x = as.numeric(predicted), y = as.numeric(corrected_q))) +
    geom_point() +
    geom_smooth(method = "lm") + 
    facet_wrap("description", scales = "free", nrow = 3) +
    labs(title = paste0("predicted response by model: ", i)) + 
    xlab("predicted_response") +
    ylab("corrected question response")+ 
    stat_cor(color = "red")
  print(plots.ls.2[[length(plots.ls.2)]])
}
dev.off()
tmp <- left_join(left_join(plots.ls, m.out%>%select(IID, paste0("model_", c(1,3,5,7,9,11,13,15,17)))), question.dict %>% rename(question = q0))
write_rds(tmp, "data/derivatives/abcd-mph-corrected-sst-w-predicted.rds")
rm(tmp)
gc()
################################################################################
# look at correlation of taking methylphenidate while considering the predicted response in model
library(ggstatsplot)
inner_join(plots.ls %>% 
  filter(question %in% filtered.q.of.int), m.out) %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  mutate(response = as.factor(ifelse(model_3 < 0, "responder", "non-responder"))) %>%
  mutate(corrected_q = as.numeric(corrected_q)) %>%
  ggplot(aes(fill = methylphenidate, y = corrected_q, x = response)) +
  geom_boxplot() +
  facet_wrap("question", scales = "free_y", ncol = 3) +
  stat_compare_means(size = 3, method = "t.test")
p2 <- inner_join(plots.ls %>% 
                   filter(question %in% filtered.q.of.int), m.out) %>%
  mutate(methylphenidate = as.factor(methylphenidate)) %>%
  mutate(response = as.factor(ifelse(model_3 < 0, "responder", "non-responder"))) %>%
  mutate(corrected_q = as.numeric(corrected_q)) %>%
  ggplot(aes(x = methylphenidate, y = corrected_q)) +
  geom_boxplot() +
  facet_wrap("question", scales = "free", ncol = 3) +
  stat_compare_means(size = 3, method = "t.test")

################################################################################
################################################################################
################################################################################
################################################################################
# dimensionality reduction for task features?
# indep section
abcd.raw.dir <- "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_4_0"
sst.raw <- read_delim(paste0(abcd.raw.dir, "/abcd_sst02.txt"), delim = "\t")
question.dict <- data.frame(q0 = colnames(sst.raw), description = t(sst.raw)[,1], row.names = NULL)
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meds-matrix_CMAP-meds.rds") %>%
  as.data.frame() %>%
  select(c(1:6), methylphenidate) %>%
  rename(IID = ID)

sst.filt <- left_join(sst.raw[-1,] %>%
  select(IID = subjectkey, interview_age, sex, eventname,
         ends_with("_rt"), 
         ends_with("_mrt"),
         tfmri_sst_beh_switchflag,tfmri_sst_beh_performflag) %>%
  filter(is.na(tfmri_sst_beh_switchflag) == F),
  abcd.meds)

x0 <- sst.filt %>% 
  select(c(1:4, 48),
         ends_with("_rt"),
         ends_with("_mrt"),
         ends_with("flag")) %>%
  mutate_at(.vars = vars(ends_with("rt"), ends_with("flag"), "interview_age"), .funs = function(x) as.numeric(x)) %>%
  mutate_at(.vars = vars(ends_with("rt")), .funs = function(x) ifelse(is.na(x), mean(x, na.rm = T), x)) %>%
  # mutate_at(.vars = vars(ends_with("_rt")), .funs = function(x) qnorm(x)) %>%
  # mutate_at(.vars = vars(ends_with("_rt")), .funs = function(x) ifelse(x=="-Inf", NA, x)) %>%
  drop_na()
x1 <- x0 %>%
  mutate_at(.vars = vars(ends_with("rt")), .funs = function(x) scale(x, center = T, scale = T)) %>%
  mutate_at(.vars = vars(ends_with("rt")), .funs = function(x) residuals(lm(x ~ scale(as.numeric(x0$interview_age), scale = T, center = T)))) %>%
  drop_na()

# AA
library(archetypes)
sst.aa <- stepArchetypes(x1[,-c(1:5)], k = 1:10)

screeplot(sst.aa)
a4 <- bestModel(sst.aa[[4]])
a4$archetypes %>% as.data.frame() %>%
  mutate(A = paste0("A", 1:4)) %>%
  pivot_longer(cols = -c("A"), names_to = "question", values_to = "loading") %>%
  ggplot(aes(x=A, y = question, fill = loading, label = ifelse(abs(loading)>2, round(loading, 2), "")))+
  geom_tile()+
  geom_text(size=3)+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1]) +
  my.guides+labs(y="", x="", caption = paste0("n(samples): ", nrow(x1), "\n",
                                              "questions with missing values were replaced by taking the mean of question values", "\n", 
                                              "question val = residuals(scale(raw_q) ~ scale(interview_age))"))
sst.archetypes <- cbind(a4$alphas, sst.filt[,1:4])
colnames(sst.archetypes)[1:4] <- paste0("A", 1:4)
write_tsv(sst.archetypes, "data/derivatives/abcd-sst-data-archetypes.tsv")
################################################################################
# correlation between abcd pgs and sst archetypes
# indep section
sst.a <- read_tsv("data/derivatives/abcd-sst-data-archetypes.tsv") %>%
  mutate_at(.vars = vars(starts_with("A")), .funs = function(x) qnorm(x)) %>%
  mutate_at(.vars = vars(starts_with("A")), .funs = function(x) ifelse(x == "-Inf", NA, x))

abcd.pgs <- read_tsv("data/derivatives/spark-abcd-corrected-pgs.tsv")
pgs.meta <- read_delim("/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/PGS/clean_PGS_names.tsv") %>%
  mutate(dup = ifelse(duplicated(PGS_shortname), 1, 0)) %>%
  mutate(pgs = str_replace_all(pattern = "-",replacement =  "_",string =  PGS_name))
colnames(abcd.pgs) <- str_replace_all(pattern = "-",replacement =  "_",string =  colnames(abcd.pgs))
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meds-matrix_CMAP-meds.rds") %>%
  as.data.frame() %>%
  select(c(1:6), methylphenidate) %>%
  rename(IID = ID)
abcd.c <- inner_join(abcd.meds%>%mutate(interview_age = as.numeric(interview_age)), sst.a %>% select(-starts_with("tfmri")))

abcd.c <- left_join(abcd.c, abcd.pgs) %>% drop_na()
corr.table(abcd.c%>% select(paste0("A", 1:4)), abcd.c%>%select(starts_with("corrected"))) %>%
  filter(V1 %in% c(paste0("A", 1:4)), ! V2 %in% c(paste0("A", 1:4))) %>%
  mutate(V2 = sub("corrected_", "", V2)) %>%
  ggplot(aes(x=V1, y=V2, fill=r, label=ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 4)), "")))+
  geom_tile() +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  geom_text(size=3)+
  my.guides + labs(x="", y="", caption = paste0("n(samples): ", nrow(abcd.c), "\n",
                                                "Archetypes weights were normalized using qnorm, and -Inf samples were dropped here"))
################################################################################
################################################################################
# check correlation between sst archetypes, methylphenidate, and predicted response
# indep section
sst.a <- read_tsv("data/derivatives/abcd-sst-data-archetypes.tsv")
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meds-matrix_CMAP-meds.rds") %>%
  as.data.frame() %>%
  select(c(1:6), methylphenidate) %>%
  rename(IID = ID)
abcd.c <- inner_join(abcd.meds%>%mutate(interview_age = as.numeric(interview_age)), sst.a %>% select(-starts_with("tfmri")))
abcd.pred <- read_rds("data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(model_3 = m)
table(abcd.pred$IID %in% abcd.c$IID)
table(abcd.c$IID %in% abcd.pred$IID)
table(abcd.c$methylphenidate)
meta <- inner_join(abcd.c, abcd.pred)

corr.table(meta %>% select(paste0("A", 1:4)), meta %>% select(methylphenidate, model_3)%>%mutate(model_3=-model_3)) %>%
  filter(V1 %in% c(paste0("A", 1:4)), ! V2 %in% c(paste0("A", 1:4))) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.1, paste0("r: ", round(r, 3), ",  p: ", round(pval, 4)), ""))) +
  geom_tile() +
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1], mid = "white") +
  geom_text(size = 3)+
  my.guides+xlab("")+ylab("")+
  labs(caption = paste0("n(samples): ", nrow(meta), "\n",
                        "n(unique samples): ", length(unique(meta$IID))))

m1 <- glm(A1 ~ methylphenidate+ model_3+(model_3*methylphenidate), data = meta %>%
               mutate(model_3 = scale(-model_3, scale=T, center=T), 
                      A1 = qnorm(A1))%>%filter(A1!="-Inf")
                      # A1 = sqrt(asin(A1)))
          , family = "gaussian")
m2 <- glm(A2 ~ methylphenidate+ model_3+(model_3*methylphenidate), data = meta %>%
               mutate(model_3 = scale(-model_3, scale=T, center=T), 
                      A2 = qnorm(A1))%>%filter(A2!="-Inf")
                      # A2 = sqrt(asin(A2)))
          , family = "gaussian")
m3 <- glm(A3 ~ methylphenidate+ model_3+(model_3*methylphenidate), data = meta %>%
               mutate(model_3 = scale(-model_3, scale=T, center=T), 
                      A3 = qnorm(A3))%>%filter(A3!="-Inf")
                      # A3 = sqrt(asin(A3)))
          , family = "gaussian")
m4 <- glm(A4 ~ methylphenidate+ model_3+(model_3*methylphenidate), data = meta %>%
               mutate(model_3 = scale(-model_3, scale=T, center=T), 
                      A4 = qnorm(A4))%>%filter(A4!="-Inf")
                      # A4 = sqrt(asin(A4)))
          , family = "gaussian")




rbind(summary(m1)$coefficients %>% as.data.frame() %>%
        rownames_to_column("var") %>%
        mutate(predicted = "A1", confint_min = confint(m1)[,1], confint_max = confint(m1)[,2]),
      summary(m2)$coefficients %>% as.data.frame() %>%
        rownames_to_column("var") %>%
        mutate(predicted = "A2", confint_min = confint(m2)[,1], confint_max = confint(m2)[,2]),
      summary(m3)$coefficients %>% as.data.frame() %>%
        rownames_to_column("var") %>%
        mutate(predicted = "A3", confint_min = confint(m3)[,1], confint_max = confint(m3)[,2]),
      summary(m4)$coefficients %>% as.data.frame() %>%
        rownames_to_column("var") %>%
        mutate(predicted = "A4", confint_min = confint(m4)[,1], confint_max = confint(m4)[,2])) %>%
  mutate(sig = ifelse(`Pr(>|t|)`<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  filter(var != "(Intercept)") %>%
  mutate(var = sub("scale", "", var)) %>%
  mutate(var = factor(var, levels = c("model_3", "methylphenidate", "methylphenidate:model_3"))) %>%
  mutate(predicted = as.factor(predicted)) %>%
  ggplot(aes(x=Estimate, y = var), color = var) +
  geom_point(aes(color = predicted, alpha = sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_shape_manual(values = c(1, 2)) + 
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, color = predicted, alpha = sig), 
                 size = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = six.colors[1:4]) +
  my.guides+labs(x ="Estimate", y="", 
                 title = "predicting SST archetypess by MPH predicted response", 
                 caption = paste0("models: glm(qnorm(A[1:4]) ~ methylphenidate + scale(model_3) + model_3*methylphenidate)", "\n",
                                  "n(samples):", "\n", 
                                  "A1 ", length(m1$residuals), "\n", 
                                  "A2 ", length(m2$residuals), "\n", 
                                  "A3 ", length(m3$residuals), "\n", 
                                  "A4 ", length(m4$residuals), "\n",
                                  "n(samples taking MPH): ", nrow(meta%>%filter(methylphenidate==1))))

# only the ones taking mph
meta %>%
  mutate(model_3 = -model_3) %>% filter(methylphenidate==1) %>%
  pivot_longer(cols = starts_with("A"), names_to = "archetype", values_to = "v") %>%
  ggplot(aes(x=model_3, y =v))+
  geom_point()+
  facet_wrap(facets = "archetype", scales = "free")+
  geom_smooth(method = "lm")+
  stat_cor()
################################################################################
################################################################################
################################################################################
# based on spark, the cog-symbol-digit pgs is predicitve of mph response
# look at the correlation between 

summary(lm(abcd.c$A2 ~ abcd.c$interview_age + abcd.c$methylphenidate + abcd.c$corrected_cog_symbol_digit_UKB_2020 + (abcd.c$methylphenidate*abcd.c$corrected_cog_symbol_digit_UKB_2020)))
summary(lm(abcd.c$A3 ~ abcd.c$interview_age + abcd.c$methylphenidate + abcd.c$corrected_cog_symbol_digit_UKB_2020 + (abcd.c$methylphenidate*abcd.c$corrected_cog_symbol_digit_UKB_2020)))
summary(lm(abcd.c$A4 ~ abcd.c$interview_age + abcd.c$methylphenidate + abcd.c$corrected_cog_symbol_digit_UKB_2020 + (abcd.c$methylphenidate*abcd.c$corrected_cog_symbol_digit_UKB_2020)))  
