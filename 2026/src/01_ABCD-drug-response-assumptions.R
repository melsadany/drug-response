################################################################################
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),"/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/2026")
setwd(project.dir)
################################################################################
################################################################################
################################################################################
################################################################################
## read medication data for ABCD
## define a person as drug responder if they take it two consec years
## and not responder if they stop it
abcd.raw.dir <- correct_path("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core")
abcd.deriv.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/data/ABCD/abcd5")
abcd.pgs <- read_tsv("../data/derivatives/spark-abcd-corrected-pgs.tsv")
# colnames(abcd.pgs) <- sub("corrected_", "", colnames(abcd.pgs))
demo <- read_csv(paste0(abcd.deriv.dir, "/age-sex-by-eventname.csv"))

abcd.meds.r <- read_rds(paste0(abcd.deriv.dir, "/../meds/abcd5/abcd5-meds-matrix.rds"))
drop.sub <- read_csv(paste0(abcd.deriv.dir, "/../meds/abcd5/subjects-missing-med-name.csv"))
abcd.meds <- anti_join(abcd.meds.r, drop.sub)
abcd.pred <- read_rds("../data/derivatives/m-outputs/abcd/all-samples/model-celltype-all-FALSE-TRUE-1.rds") %>%
  rename(predicted = m) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])
med.count <- abcd.meds %>% select(-c(2,3)) %>% group_by(IID) %>% summarise_all(.funs=sum) %>% ungroup()


#### methylphenidate
mph.equv <- c("methylphenidate","ritalin","concerta", "methylin", "aptensio", 
              "daytrana","quillivant", "quillichew", "jornay", "adhansia", 
              "metadate","relexxii") # https://www.drugs.com/medical-answers/brands-methylphenidate-3510739/
abcd.mph <- abcd.meds %>% 
  rowwise() %>% mutate(methylphenidate = methylphenidate+ritalin+concerta+methylin+aptensio+
                         daytrana+quillivant+quillichew+jornay+metadate+relexxii) %>% ungroup() %>%
  select(1,2, methylphenidate) %>% 
  mutate(other = rowSums(abcd.meds[,-c(1:3)])-methylphenidate, MPH_only = (methylphenidate>0 & other==0))
med.count.mph <- abcd.mph %>% select(-2) %>% group_by(IID,MPH_only) %>% summarise_all(.funs=sum) %>% ungroup() %>%
  filter(methylphenidate !=0, other<1)
#### clonidine
abcd.cd <- abcd.meds %>% select(1,2, clonidine) %>% 
  mutate(other = rowSums(abcd.meds[,-c(1:3)])-clonidine, CL_only = (clonidine>0 & other==0))
med.count.cd <- abcd.cd %>% select(-2) %>% group_by(IID,CL_only) %>% summarise_all(.funs=sum) %>% ungroup() %>%
  filter(clonidine !=0, other<1)
################################################################################
################################################################################
inner_join(med.count.mph,abcd.pred) %>% 
  mutate(MPH_c = ifelse(methylphenidate>1, "two or more", "one"),
         MPH_c = factor(MPH_c,c("one","two or more"))) %>%
  ggplot(aes(MPH_c, predicted, fill = MPH_c))+
  geom_violin(show.legend = F) + geom_boxplot(width=0.2,fill="white") + 
  ggpubr::stat_compare_means(color="red",label.y.npc = 1,label.x.npc = 0.13) +
  scale_fill_manual(values = palette.1)+bw.theme+
  labs(x="number of years on MPH", y="predicted MPH responsiveness")
ggsave2("figs/abcd-predicted-MPH-response-vs-taking-mph-once-or-more.png",3,4)

inner_join(abcd.pgs %>% select(IID, contains("ADHD")),abcd.pred) %>% 
  filter(IID %in% unique(abcd.mph$IID[abcd.mph$methylphenidate!=0])) %>%
  pivot_longer(cols = -c(IID,predicted)) %>%
  ggplot(aes(value, predicted,color=name))+
  geom_point(shape=1)+geom_smooth(method="lm",se=F) + ci_ribbon.multi+
  ggpubr::stat_cor() + 
  bw.theme+labs(x="PGS", y="predicted MPH responsiveness")

################################################################################
################################################################################
## identify participants of interest
abcd.mph.of.int <- abcd.mph %>% 
  mutate(mph_broad = methylphenidate>0) %>% group_by(IID) %>%
  filter(any(mph_broad)) %>% ungroup()
abcd.mph.on.all <- abcd.mph.of.int %>% filter(mph_broad) # some participants have been on MPH all data points
abcd.mph.off <- abcd.mph.of.int %>% filter(!mph_broad)
table(abcd.mph.on.all$IID %in% abcd.mph.off$IID) # here
table(abcd.mph.off$IID %in% abcd.mph.on.all$IID)
abcd.mph.on <- abcd.mph.on.all %>% filter(IID %in% abcd.mph.off$IID)
table(abcd.mph.on$IID %in% abcd.mph.off$IID)

## timepoints on MPH
abcd.mph.long <- abcd.mph %>% 
  mutate(eventname=case_when(eventname=="baseline_year_1_arm_1"~"E0",
                             eventname=="1_year_follow_up_y_arm_1"~"E1",
                             eventname=="2_year_follow_up_y_arm_1"~"E2",
                             eventname=="3_year_follow_up_y_arm_1"~"E3",
                             eventname=="4_year_follow_up_y_arm_1"~"E4", T~"EE"),
         mph_broad = as.numeric(methylphenidate >0)) %>%
  pivot_wider(names_from = eventname, values_from = mph_broad, id_cols = IID) %>%
  select(IID, paste0("E",c(0:4))) %>% mutate_all(.funs = function(x) replace_na(x,0)) %>%
  mutate(any_MPH = (E0+E1+E2+E3+E4)>0) %>% filter(any_MPH) %>%
  mutate(consec_01 = E0+E1==2, consec_12 = E1+E2==2, consec_23 = E2+E3==2, consec_34 = E3+E4==2) %>%
  mutate(responder_2 = (E0+E1+E2+E3+E4)>1, responder_3 = (E0+E1+E2+E3+E4)>2, responder_4 = (E0+E1+E2+E3+E4)>3)
################################################################################
## get CBCL deltas

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
spark.df.raw <- read_delim("/sdata/Simons/SPARK/ResearchMatch/sleep_eating_gastro_2018/RM0018Michaelson_AdultSurvey20181031_FinalDataSetTabDelimited_BV_Age.csv")
spark.pred <- read_rds("../data/derivatives/m-outputs/spark/all-spark-samples/model-celltype-all-TRUE-TRUE-1-TT.rds") %>%
  rename(predicted = m) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])

spark.mph <- spark.df.raw %>% 
  mutate(form_completed=form_completed==1) %>% filter(form_completed) %>%
  select(ParticipantID, tried_MPH = q119_sq01_alias,
         which_MPH = q119_sq02_alias, currently_on_MPH = q119_sq03_alias,
         effective_MPH = q119_sq04_alias, MPH_side_effects=q119_sq05_alias) %>% 
  mutate_at(.vars = vars(c(tried_MPH, currently_on_MPH, effective_MPH)), .funs=function(x)ifelse(x==0, T, ifelse(x==1,F,NA)))
inner_join(spark.mph%>%rename(var=effective_MPH),spark.pred%>%rename(ParticipantID=IID)) %>% 
  drop_na(var) %>%
  # mutate(MPH_c = ifelse(methylphenidate>1, "two or more", "one"),
  #        MPH_c = factor(MPH_c,c("one","two or more"))) %>%
  ggplot(aes(var, predicted, fill = var))+
  geom_violin(show.legend = F) + geom_boxplot(width=0.2,fill="white") + 
  ggpubr::stat_compare_means(color="red",label.y.npc = 1,label.x.npc = 0.13) +
  scale_fill_manual(values = palette.1)+bw.theme+
  labs(x="number of years on MPH", y="predicted MPH responsiveness")
ggsave2("figs/abcd-predicted-MPH-response-vs-taking-mph-once-or-more.png",3,4)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
