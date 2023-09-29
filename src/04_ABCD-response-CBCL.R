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
abcd.meds <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meph-meta_by-eventname.rds") %>%
  filter(base_meph == 1 | year1_meph == 1 | year2_meph == 1 | year3_meph == 1)
################################################################################
load("data/derivatives/m-outputs/abcd/combined-output.RData")
m.output <- to.save[[1]]
rm(to.save)
gc()
################################################################################
################################################################################
################################################################################
######################## correlate response to CBCL att ########################
################################################################################
################################################################################
################################################################################
# get correlation between taking MPH and att_problems total in ABCD
cbcl.sclaes <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/cbcl-scales.rds") %>%
  filter(scale == "att")
abcd.cbcl <- read_tsv("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_4_0/abcd_cbcl01.txt")[-1,] %>%
  dplyr::select(IID = subjectkey, starts_with("interview"), sex, eventname, cbcl.sclaes$q0) %>%
  drop_na() %>%
  as.data.frame() %>%
  filter(IID %in% m.output$IID)
colnames(abcd.cbcl)[6:ncol(abcd.cbcl)] <- cbcl.sclaes$cbcl_item
abcd.cbcl <- cbind(abcd.cbcl[c(1,2,4,5)], abcd.cbcl[c(3,6:ncol(abcd.cbcl))]%>%mutate_all(function(x) as.numeric(x)))
abcd.cbcl$att_tot <- rowSums(abcd.cbcl[,6:ncol(abcd.cbcl)])
abcd.cbcl.filt <- abcd.cbcl %>% 
  select(-starts_with("q")) %>%
  mutate(sex = as.factor(sex))
jtools::export_summs(lm(abcd.cbcl.filt$att_tot ~ abcd.cbcl.filt$interview_age), lm(abcd.cbcl.filt$att_tot ~ abcd.cbcl.filt$sex))

# check correlation between tot number of meds and att_tot
meds.cbcl.meta <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meds-matrix_CMAP-meds.rds") %>%
  as.data.frame() %>%
  filter(methylphenidate == 1)
meds.cbcl.meta.2 <- meds.cbcl.meta%>%
  select(1:6,colnames(meds.cbcl.meta[,7:ncol(meds.cbcl.meta)])[colSums(meds.cbcl.meta[,7:ncol(meds.cbcl.meta)])>0]) %>%
  mutate(tot_meds = rowSums(meds.cbcl.meta[,7:ncol(meds.cbcl.meta)]))

abcd.cbcl.filt.2 <- left_join(abcd.cbcl.filt, meds.cbcl.meta.2 %>%select(IID=ID, eventname, tot_meds)) %>%
  mutate(tot_meds = ifelse(is.na(tot_meds), 0, tot_meds))
jtools::export_summs(lm(abcd.cbcl.filt.2$att_tot ~ abcd.cbcl.filt.2$interview_age), 
                     lm(abcd.cbcl.filt.2$att_tot ~ abcd.cbcl.filt.2$sex), 
                     lm(abcd.cbcl.filt.2$att_tot ~ abcd.cbcl.filt.2$tot_meds))
# turns out the age, sex, and tot_meds affect att_tot
abcd.cbcl.filt.2$corrected_att_tot = residuals(lm(abcd.cbcl.filt.2$att_tot ~ abcd.cbcl.filt.2$sex + abcd.cbcl.filt.2$interview_age + abcd.cbcl.filt.2$tot_meds))

################################################################################
# check correlation between taking mph and cbcl_tot
correlations <- data.frame(ID = unique(abcd.cbcl.filt.2$IID), 
                           att_adjr2=NA,
                           age_sex_meds_corr_att_adjr2=NA,
                           att_r=NA,
                           age_sex_meds_corr_att_r=NA,
                           mph_tot = NA,
                           event_tot=NA)
meds.cbcl.meta <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/meds/meds-matrix_CMAP-meds.rds") %>%
  as.data.frame() %>%
  filter(methylphenidate == 1) %>%
  select(1:6, methylphenidate)

library(doMC)
registerDoMC(cores = 8)
for (i in 1:nrow(correlations)){
# correlations12 <- foreach (i = 1:nrow(correlations), .combine = rbind) %dopar% {
  # i=1
  subject <- correlations$ID[i]
  meph.df <- meds.cbcl.meta %>%
    filter(ID == subject) %>%
    select(ID, eventname, interview_date, methylphenidate)
  df.sub <- left_join(abcd.cbcl.filt.2 %>%
                        filter(IID==subject),
                      meph.df %>% rename(IID=ID, meph_stat=methylphenidate) %>% select(-interview_date)) %>%
    mutate(interview_age=as.numeric(interview_age)) %>%
    mutate_at(.vars = "meph_stat",funs(replace_na(.,0)))
  correlations$att_adjr2[i] <- summary(lm(att_tot ~ meph_stat, data = df.sub))$adj.r.squared
  correlations$age_sex_meds_corr_att_adjr2[i] <- summary(lm(corrected_att_tot ~ meph_stat , data = df.sub))$adj.r.squared
  correlations$att_r[i] <- cor(df.sub$att_tot, df.sub$meph_stat)
  correlations$age_sex_meds_corr_att_r[i] <- cor(df.sub$corrected_att_tot, df.sub$meph_stat)
  correlations$mph_tot[i] <- sum(df.sub$meph_stat)
  correlations$event_tot[i] <- nrow(df.sub)
  print(i)
  # return(correlations[i,])
}
################################################################################
# compare predicted response with r2
meta <- left_join(m.output %>% 
                    mutate(model_3 = -model_3) %>%
                    rename(cont_predicted_response = model_3),
                  correlations %>% rename(IID = ID)) %>%
  mutate(mph_prop = mph_tot/event_tot) %>%
  filter(mph_prop!=1) %>%
  # filter(mph_tot >=2) %>%
  # filter(mph_prop == 0.25) %>%
  mutate(mph_prop=as.factor(round(mph_prop,3)))

meta$age_sex_meds_corr_att_r_corr_prop = residuals(lm(meta$age_sex_meds_corr_att_r ~ meta$mph_prop))

summary(lm(meta$att_r ~ meta$cont_predicted_response))
summary(lm(meta$age_sex_meds_corr_att_r ~ meta$cont_predicted_response))
summary(lm(meta$att_adjr2 ~ meta$cont_predicted_response))
summary(lm(meta$age_sex_meds_corr_att_adjr2 ~ meta$cont_predicted_response))

p1 <- meta %>% ggplot(aes(x=cont_predicted_response, y=att_r, color = mph_prop)) +
  geom_point() +
  geom_smooth(method = "lm", show.legend = F) +
  facet_wrap("mph_prop")+
  labs(subtitle = "cor(att_tot, meph_stat, data = df.sub)")
p2 <- meta %>% ggplot(aes(x=cont_predicted_response, y=age_sex_meds_corr_att_r, color = mph_prop)) +
  geom_point() +
  geom_smooth(method = "lm", show.legend = F) +
  facet_wrap("mph_prop")+
  labs(subtitle = paste0("cor(corrected_att_tot, meph_stat, data = df.sub)", "\n",
                         "df.sub$corrected_att_tot = residuals(lm(att_tot ~ sex + interview_age + tot_meds, data = abcd.cbcl.filt.2))"))
p3 <- meta %>% ggplot(aes(x=cont_predicted_response, y=att_adjr2, color = mph_prop)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap("mph_prop") +
  labs(subtitle = "summary(lm(att_tot ~ meph_stat, data = df.sub))$adj.r.squared",
    caption = paste0("nsamples(mph_prop (0.25)): ", table(meta$mph_prop == "0.25")[[2]], "\n",
                        "nsamples(mph_prop (0.333)): ", table(meta$mph_prop == "0.333")[[2]], "\n",
                        "nsamples(mph_prop (0.5)): ", table(meta$mph_prop == "0.5")[[2]], "\n",
                        "nsamples(mph_prop (0.667)): ", table(meta$mph_prop == "0.667")[[2]], "\n",
                        "nsamples(mph_prop (0.75)): ", table(meta$mph_prop == "0.75")[[2]]))
p4 <- meta %>% ggplot(aes(x=cont_predicted_response, y=age_sex_meds_corr_att_adjr2, color = mph_prop)) +
  geom_point() +
  geom_smooth(method = "lm", show.legend = F) +
  facet_wrap("mph_prop") +
  labs(subtitle = paste0("summary(lm(corrected_att_tot ~ meph_stat , data = df.sub))$adj.r.squared", "\n",
                         "df.sub$corrected_att_tot = residuals(lm(att_tot ~ sex + interview_age + tot_meds, data = abcd.cbcl.filt.2))"))
p5 <- meta %>% ggplot(aes(x=cont_predicted_response, y=age_sex_meds_corr_att_r_corr_prop, color = mph_prop)) +
  geom_point() +
  geom_smooth(method = "lm", show.legend = F) +
  facet_wrap("mph_prop") +
  labs(subtitle = "age_sex_meds_corr_att_r_corr_prop = residuals(lm(meta$age_sex_meds_corr_att_r ~ meta$mph_prop))")

patchwork::wrap_plots(p1+stat_cor(),p2+stat_cor(),p3+stat_cor(),p4+stat_cor(), nrow = 2)
patchwork::wrap_plots(p1,p2,p5,p3,p4, nrow = 2)
patchwork::wrap_plots(p1+stat_cor(color = "black"),p2+stat_cor(color = "black"),p5+stat_cor(color = "black"),p3+stat_cor(color = "black"),p4+stat_cor(color = "black"), nrow = 2)

################################################################################


################################################################################

################################################################################