################################################################################
#              check correlation between CBCL and fMRI performance             #
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
################################################################################
load("data/derivatives/m-outputs/abcd/combined-output.RData")
m.output <- to.save[[1]]
rm(to.save)
gc()
################################################################################
# get the MID fMRI task data
mid.tmp <- read_rds("data/derivatives/abcd-mph-corrected-mid-w-predicted.rds") %>%
  select(IID, question, eventname, corrected_q) %>%
  pivot_wider(names_from = "question", values_from = "corrected_q")

abcd.mid.mph <- inner_join(right_join(m.output%>%select(IID, paste0("model_", c(1,3,5,7,9,11,13,15))),right_join(raw.meds, mid.tmp) %>%
                                       mutate(methylphenidate = ifelse(is.na(methylphenidate),0,methylphenidate))), 
                          raw.meds %>% select(IID, eventname, tot_meds)) %>% 
  mutate(tot_meds = ifelse(is.na(tot_meds), 0, tot_meds))

mid.q.of.int <- c( #smal reward
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_srw_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_srwpfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_srwpfb_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_srwnfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_srwnfb_mrt"),
  # large reward
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_lrw_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_lrwpfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_lrwpfb_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_lrwnfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_lrwnfb_mrt"),
  # small loss
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_sl_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_slpfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_slpfb_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_slnfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_slnfb_mrt"),
  # large loss
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_ll_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_llpfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_llpfb_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_llnfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_llnfb_mrt"),
  # neutral trials
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_nt_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_ntpfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_ntpfb_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_ntnfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_ntnfb_mrt"),
  # large and small rewards
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hrw_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hrwpfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hrwpfb_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hrwnfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hrwnfb_mrt"),
  # large and small losses
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hl_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hlpfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hlpfb_mrt"),paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hlnfb_rate"),
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_hlnfb_mrt"),
  # total earnings
  paste0("tfmri_mid_", c("all", "r1", "r2"),"_beh_t_earnings"))

################################################################################
################################################################################
# get the SST fMRI task data
sst.tmp <- read_rds("data/derivatives/abcd-mph-corrected-sst-w-predicted.rds") %>%
  select(IID, question, eventname, corrected_q) %>%
  pivot_wider(names_from = "question", values_from = "corrected_q")

abcd.sst.mph <- inner_join(right_join(m.output%>%select(IID, paste0("model_", c(1,3,5,7,9,11,13,15))),right_join(raw.meds, sst.tmp) %>% 
                                       mutate(methylphenidate = ifelse(is.na(methylphenidate),0,methylphenidate))), 
                          raw.meds %>% select(IID, eventname, tot_meds)) %>% 
  mutate(tot_meds = ifelse(is.na(tot_meds), 0, tot_meds))

sst.q.of.int <- c("tfmri_sst_all_beh_nrgo_rt","tfmri_sst_r1_beh_nrgo_rt","tfmri_sst_r2_beh_nrgo_rt",
                       "tfmri_sst_all_beh_crgo_rt","tfmri_sst_r1_beh_crgo_rt","tfmri_sst_r2_beh_crgo_rt",
                       "tfmri_sst_all_beh_crlg_rt","tfmri_sst_r1_beh_crlg_rt","tfmri_sst_r2_beh_crlg_rt",
                       "tfmri_sst_all_beh_crlg_mrt","tfmri_sst_r1_beh_crlg_mrt","tfmri_sst_r2_beh_crlg_mrt",
                       "tfmri_sst_all_beh_incrlg_rt","tfmri_sst_r1_beh_incrlg_rt","tfmri_sst_r2_beh_incrlg_rt",
                       "tfmri_sst_beh_performflag", "tfmri_sst_beh_switchflag")
################################################################################
################################################################################
# get cbcl
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
summary(lm(abcd.cbcl.filt$att_tot ~ abcd.cbcl.filt$sex))
summary(lm(abcd.cbcl.filt$att_tot ~ abcd.cbcl.filt$interview_age))
abcd.cbcl.filt$corrected_att_tot <- residuals(lm(abcd.cbcl.filt$att_tot ~ abcd.cbcl.filt$interview_age + abcd.cbcl.filt$sex))
################################################################################
################################################################################
# check correlation between corrected attention_tot and MID
mid.cbcl <- inner_join(abcd.cbcl.filt, abcd.mid.mph%>%select(-c(interview_age,sex, interview_date))) %>%
  filter(methylphenidate==1)
p1 <- corr.table(mid.cbcl%>%select(att_tot, corrected_att_tot, "model_3")%>%mutate(model_3=-model_3), mid.cbcl%>%select(mid.q.of.int)) %>%
  filter(V1%in%c("att_tot", "corrected_att_tot", "model_3")) %>%
  filter(!V2 %in% c("mph_effect", "att_tot", "corrected_att_tot", "model_3")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label=ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 4)), ""))) +
  geom_tile()+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  geom_text(size=3) +
  xlab("")+ylab("")+my.guides +
  labs(caption = paste0("n(samples): ", nrow(mid.cbcl), "\n", "n(unique subjects): ", length(unique(mid.cbcl$IID))))
################################################################################
################################################################################
# check correlation between attention_tot and SST
sst.cbcl <- inner_join(abcd.cbcl.filt, abcd.sst.mph%>%select(-c(interview_age,sex, interview_date))) %>%
  filter(methylphenidate==1)
p2 <- corr.table(sst.cbcl%>%select(att_tot, corrected_att_tot, "model_3")%>%mutate(model_3=-model_3), sst.cbcl%>%select(sst.q.of.int)) %>%
  filter(V1%in%c("att_tot", "corrected_att_tot", "model_3")) %>%
  filter(!V2 %in% c("mph_effect", "att_tot", "corrected_att_tot", "model_3")) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label=ifelse(pval<0.05, paste0("r: ", round(r, 3), ",  p: ", round(pval, 4)), ""))) +
  geom_tile()+
  scale_fill_gradient2(low = redblu.col[2], high = redblu.col[1])+
  geom_text(size=3) +
  xlab("")+ylab("")+my.guides +
  labs(caption = paste0("n(samples): ", nrow(sst.cbcl), "\n", "n(unique subjects): ", length(unique(sst.cbcl$IID))))
################################################################################
################################################################################
patchwork::wrap_plots(p1,p2)
