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
  filter(methylphenidate == 1) %>%
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
######################## correlate response to MID fMRI ########################
################################################################################
################################################################################
################################################################################
abcd.raw.dir <- "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_4_0"
mid.raw <- read_delim(paste0(abcd.raw.dir, "/abcd_mid02.txt"), delim = "\t")
question.dict <- data.frame(q0 = colnames(mid.raw), description = t(mid.raw)[,1], row.names = NULL)

filtered.q.of.int <- c( #smal reward
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_srw_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_srwpfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_srwpfb_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_srwnfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_srwnfb_mrt"),
                       # large reward
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_lrw_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_lrwpfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_lrwpfb_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_lrwnfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_lrwnfb_mrt"),
                       # small loss
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_sl_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_slpfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_slpfb_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_slnfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_slnfb_mrt"),
                       # large loss
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_ll_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_llpfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_llpfb_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_llnfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_llnfb_mrt"),
                       # neutral trials
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_nt_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_ntpfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_ntpfb_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_ntnfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_ntnfb_mrt"),
                       # large and small rewards
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hrw_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hrwpfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hrwpfb_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hrwnfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hrwnfb_mrt"),
                       # large and small losses
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hl_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hlpfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hlpfb_mrt"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hlnfb_rate"),
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hlnfb_mrt"),
                       # total earnings
                       paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_t_earnings"))


mid.filt <- mid.raw %>%
  filter(subjectkey %in% m.out$IID) %>% 
  select(IID = subjectkey, interview_age, sex, eventname,
         # c(10,11,23:25, 26,28,30, 35:37, 38,40,42, 47:49, 50,52,54, 59:61, 62,64,66, 71:73,
         #   77:79, 83:85, 86,88,90, 95:97, 98:101)
         filtered.q.of.int) %>%
  filter(is.na("tfmri_mid_beh_switchflag") == F)
mid.mph.filt <- left_join(right_join(m.out,right_join(abcd.meds, mid.filt) %>% 
                                       # filter(methylphenidate %in% c(0,1))
                                       mutate(methylphenidate = ifelse(is.na(methylphenidate),0,methylphenidate))), 
                          raw.meds %>% select(IID, eventname, tot_meds)) %>% 
  mutate(tot_meds = ifelse(is.na(tot_meds), 0, tot_meds))

# write_rds(mid.mph.filt, "data/derivatives/abcd-mph-mid.rds", compress = "gz")

# look at correlation between age, number of meds taken and task performance in general
age.effect <- do.call(rbind, lapply(mid.mph.filt[,filtered.q.of.int], function(x) summary(lm(x ~ as.numeric(mid.mph.filt$interview_age)))$coefficients%>%
                       as.data.frame()%>%rownames_to_column("var"))) %>% rownames_to_column("question") %>%
  mutate(question = sub("\\.[1-9]", "", question)) %>%
  filter(var != "(Intercept)")%>%
  mutate(age_sig = ifelse(`Pr(>|t|)`<0.05, T,F))

meds.effect <- do.call(rbind, lapply(mid.mph.filt[,filtered.q.of.int], function(x) summary(lm(x ~ as.numeric(mid.mph.filt$tot_meds)))$coefficients%>%
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
    df <- mid.mph.filt %>% select(IID, interview_age, tot_meds, which(colnames(mid.mph.filt)==que), methylphenidate, eventname) %>%
      drop_na() %>%
      rename(q = 4) %>%
      mutate(corrected_q = residuals(lm(q ~ as.numeric(interview_age) + tot_meds))) %>% 
      mutate(question = que) %>% mutate(age_sig = T, meds_sig = T) 
    plots.ls[[i]] <- df
  } else if (covar$meds_sig[i] == T) {
    df <- mid.mph.filt %>% select(IID, interview_age, tot_meds, which(colnames(mid.mph.filt)==que), methylphenidate, eventname) %>%
      drop_na() %>%
      rename(q = 4) %>%
      mutate(corrected_q = residuals(lm(q ~ tot_meds))) %>% 
      mutate(question = que) %>% mutate(age_sig = F, meds_sig = T) 
    plots.ls[[i]] <- df
  } else if (covar$age_sig[i] == T) {
    df <- mid.mph.filt %>% select(IID, interview_age, tot_meds, which(colnames(mid.mph.filt)==que), methylphenidate, eventname) %>%
      drop_na() %>%
      rename(q = 4) %>%
      mutate(corrected_q = residuals(lm(q ~ as.numeric(interview_age)))) %>% 
      mutate(question = que) %>% mutate(age_sig = T, meds_sig = F) 
    plots.ls[[i]] <- df
  }else {
    df <- mid.mph.filt %>% select(IID, interview_age, tot_meds, which(colnames(mid.mph.filt)==que), methylphenidate, eventname) %>%
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
plots.ls.2 <- list()
pdf("figs/11_abcd-mid-mph-taking-model-3-performance.pdf", width = 6, height = 6)
for (i in c(3)) {
  # i=1
  for (j in 1:length(filtered.q.of.int)) {
    # j=1
    q <- filtered.q.of.int[j]
    pp <- left_join(left_join(plots.ls, m.out%>%select(IID, paste0("model_", i))), question.dict %>% rename(question = q0)) %>% 
      filter(question == filtered.q.of.int[j]) %>%
      filter(methylphenidate == 1) %>%
      pivot_longer(cols = starts_with("model"), names_to = "model", values_to = "predicted") %>%
      mutate(predicted = -predicted) %>%
      ggplot(aes(x = as.numeric(predicted), y = as.numeric(corrected_q))) +
      geom_point() +
      geom_smooth(method = "lm") + 
      # facet_wrap("description", scales = "free", nrow = 3) +
      labs(title = paste0("predicted response by model: ", i), subtitle = q) + 
      xlab("predicted_response") +
      ylab("corrected question response")+ 
      stat_cor(color = "red")
    print(pp)
  }
}
dev.off()
tmp <- left_join(left_join(plots.ls, m.out%>%select(IID, paste0("model_", c(1,3,5,7,9,11,13,15,17)))), question.dict %>% rename(question = q0))
write_rds(tmp, "data/derivatives/abcd-mph-corrected-mid-w-predicted.rds")
rm(tmp)
gc()
################################



####################
################################################################################


################################################################################

################################################################################