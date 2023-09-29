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
ukb.dir <- "/Dedicated/jmichaelson-sdata/UK_Biobank/phenotypes"
################################################################################
pheno <- read_csv(paste0(ukb.dir, "/demographic_health_and_brain_data_dump.csv"))
exc.samp <- read_csv(paste0(ukb.dir, "/../samples-to-exclude.csv"), col_names = NULL)
pheno.dict <- data.frame(q_raw = colnames(pheno))
pheno.filt <- pheno %>%
  filter(!participant_id %in% exc.samp$X1) %>%
  select(participant_id, age_at_recruitment, sex, treatment_medication_code, 
         # mean_time_to_correctly_identify_matches, duration_to_complete_numeric_path_trail_number_1,
         # maximum_digits_remembered_correctly_2,
         # number_of_symbol_digit_matches_made_correctly_2,
         # duration_to_complete_numeric_path_trail_number_1_2,
         pheno.dict$q_raw[26:46], -c(pheno.dict$q_raw[c(37,38)])) %>%
  mutate(clonidine = ifelse(grepl("clonidine", treatment_medication_code, ignore.case = T, fixed = F), 1, 0))
demo <- read_csv(paste0(ukb.dir, "/demographic_info.csv"))
# ################################################################################
# # done once
# # create a fam file for sample of interest
# uk.plink.dir <- paste0(ukb.dir, "/../genotypes/imputed_v3/plink")
# fam <- read_delim(paste0(uk.plink.dir, "/ukb_imp_chr1_v3.fam"), col_names = NULL)
# sample <- uk.c %>%
#   filter(clonidine==1 | !is.na(trailMaking_duration_to_complete_numeric_path))
# fam.2 <- fam %>%
#   filter(X1 %in% sample$participant_id)
# write_delim(fam.2, "data/derivatives/uk-sample.fam", col_names = F)
# ################################################################################
# check correlation between predicted response and participants performance
drug.resp <- read_rds("data/derivatives/uk-clonidine-model-3.rds")
uk.c <- inner_join(pheno.filt %>% mutate(IID=as.character(participant_id)) %>%
                          select(-c(treatment_medication_code, participant_id)), 
                      drug.resp%>%rename(predicted=m)) %>%
  mutate(predicted = scale(-predicted, scale = T, center = T)[,1])

# uk.meta <- uk.c[rowSums(is.na(uk.meta.b))<16,]
uk.meta <- uk.c

missing.data.stats <- data.frame(var = colnames(uk.meta), 
                                 found = apply(uk.meta, MARGIN = 2, FUN = function(x) table(is.na(x))[[1]]),
                                 missing = apply(uk.meta, MARGIN = 2, FUN = function(x) ifelse(any(is.na(x)),table(is.na(x))[[2]], 0)), 
                                 clonidine_taking_found = apply(uk.meta, MARGIN = 2, FUN = function(x) sum(uk.meta%>%filter(!is.na(x))%>%select(clonidine))),
                                 clonidine_taking_missing = apply(uk.meta, MARGIN = 2, FUN = function(x) sum(uk.meta%>%filter(is.na(x))%>%select(clonidine))))
missing.data.stats %>% 
  pivot_longer(cols = c(ends_with("found"), ends_with("missing")), names_to = "miss", values_to = "count") %>%
  filter(!var%in%c("clonidine", "sex", "participant_id", "age_at_recruitment")) %>%
  ggplot(aes(x=miss, y=count))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = boxplot.colors)+
  geom_text(aes(label = count), size = 3) +
  facet_wrap("var", scales = "free_y", ncol = 5)


uk.meta %>% pivot_longer(cols = colnames(uk.meta)[3:21], names_to = "task", values_to = "score") %>%
  ggplot(aes(x=score))+
  geom_histogram()+
  facet_wrap("task", scales = "free")

# coeffs <- do.call(rbind, apply(uk.meta%>%
#                                  select(fluid_intelligence_score, mean_time_to_correctly_identify_matches),
#                                  # select(3:21), 
#                                MARGIN = 2, 
#       FUN = function(x) {m= glm(y ~ clonidine + predicted + clonidine:predicted + age + sex + age:sex, 
#                                 data = uk.meta %>% 
#                                   mutate(y=scale(x, scale=T, center = T)[,1], 
#                                          age = age_at_recruitment) %>% 
#                                   filter(!is.na(y)))
#       data.frame(summary(m)$coefficients, confint_min = confint(m)[,1], confint_max = confint(m)[,2]) %>%
#         rownames_to_column("var")
#       })) %>%
#   rownames_to_column("V2")

m1 <- glm(y ~ clonidine + predicted + clonidine:predicted + age + sex + age:sex, 
          data = uk.meta %>% 
            mutate(y=fluid_intelligence_score, 
                   age = age_at_recruitment) %>% 
            filter(!is.na(y)), family = poisson())
m2 <- glm(y ~ clonidine + predicted + clonidine:predicted + age + sex + age:sex, 
          data = uk.meta %>% 
            mutate(y=scale(mean_time_to_correctly_identify_matches, scale = T, center = T), 
                   age = age_at_recruitment) %>% 
            filter(!is.na(y)), family = gaussian())
coeffs <- rbind(data.frame(summary(m1)$coefficients, 
                           confint_min = confint(m1)[,1], 
                           confint_max = confint(m1)[,2]) %>%
                  rownames_to_column("var") %>%
                  rename(statistic = 4, pvalue = 5) %>%
                  mutate(V2 = "fluid_intelligence_score"), 
                data.frame(summary(m2)$coefficients, 
                           confint_min = confint(m2)[,1], 
                           confint_max = confint(m2)[,2]) %>%
                  rownames_to_column("var") %>%
                  rename(statistic = 4, pvalue = 5) %>%
                  mutate(V2 = "mean_time_to_correctly_identify_matches"))


# coeffs <- do.call(rbind, apply(uk.meta%>%
#                                  select(fluid_intelligence_score, mean_time_to_correctly_identify_matches),
#                                # select(3:21), 
#                                MARGIN = 2, 
#                                FUN = function(x) {m= glm(y ~ clonidine + predicted + clonidine:predicted + age + sex + age:sex, 
#                                                          data = uk.meta %>% 
#                                                            mutate(y=scale(x, scale=T, center = T)[,1], 
#                                                                   age = age_at_recruitment) %>% 
#                                                            filter(!is.na(y)))
#                                data.frame(summary(m)$coefficients, confint_min = confint(m)[,1], confint_max = confint(m)[,2]) %>%
#                                  rownames_to_column("var")
#                                })) %>%
#   rownames_to_column("V2")

coeffs %>%
  mutate(V2 = sub("\\.[1-7]", "", V2)) %>%
  filter(var!="(Intercept)")%>%
  # mutate(sig=ifelse(Pr...t..<0.05, "pval < 0.05", "pval \u2265 0.05"))%>%
  mutate(sig=ifelse(pvalue<0.05, "pval < 0.05", "pval \u2265 0.05"))%>%
  mutate(var = factor(var, levels = c("clonidine", "predicted", "clonidine:predicted",
                                      "age", "sexMale", "age:sexMale"))) %>%
  mutate(V2 = as.factor(V2)) %>%
  ggplot(aes(x=Estimate, y = var), color = var) +
  geom_point(aes(color = var, alpha = sig),  position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.2, color = "red") +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_shape_manual(values = c(1, 2)) + 
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, color = var, alpha = sig), 
                 size = 0.8, height = 0, show.legend = F, 
                 position = position_dodge(width = 0.6)) +
  scale_color_manual(values = six.colors) +
  facet_wrap("V2", scales = "free_x") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  my.guides+labs(x ="Estimate", y="", 
                 title = "predicting cog phenotypes by clonidine predicted response", 
                 caption = paste0("model used for getting estimates: ",
                                  "glm(y ~ clonidine + predicted + clonidine:predicted + age + sex + age:sex)", "\n", 
                                  "The y variable is fluid_intelligence_score/mean_time_to_correctly_identify_matches per model. ", "\n", 
                                  "The 'family' paramter in glm() was changed to be poisson for fluid_Intelligence_score", "\n",
                                  "n(samples for fluid_intelligence_score): ", nrow(uk.meta%>%filter(!is.na(fluid_intelligence_score))), 
                                  ", n(taking clonidine): ", nrow(uk.meta%>%filter(!is.na(fluid_intelligence_score),clonidine==1)), "\n",
                                  "n(samples for mean_time_to_correctly_identify_matches): ", nrow(uk.meta%>%filter(!is.na(mean_time_to_correctly_identify_matches))),
                                  ", n(taking clonidine): ", nrow(uk.meta%>%filter(!is.na(mean_time_to_correctly_identify_matches),clonidine==1))))
################################################################################
# supplementary figs
# cog tasks
p1 <- uk.meta %>% pivot_longer(cols = c(fluid_intelligence_score, mean_time_to_correctly_identify_matches), 
                               names_to = "var", values_to = "value") %>%
  ggplot(aes(x=value))+
  geom_histogram()+
  facet_wrap("var", scales = "free")+
  labs(x="")
p2 <- uk.meta %>% pivot_longer(cols = c(fluid_intelligence_score, mean_time_to_correctly_identify_matches), 
                               names_to = "var", values_to = "value") %>%
  select(var, value, age = age_at_recruitment, sex, IID) %>%
  drop_na() %>%
  ggplot(aes(x=age))+
  geom_histogram()+
  facet_wrap("var", scales = "free")+
  labs(x="", caption = paste0("n(samples) in fluid_intelligence_score: ", nrow(uk.meta %>% filter(!is.na(fluid_intelligence_score))), "\n",
                              "n(samples) in mean_time_to_correctly_identify_matches: ", nrow(uk.meta %>% filter(!is.na(mean_time_to_correctly_identify_matches)))))
patchwork::wrap_plots(p1,p2, nrow = 2)  
################################################################################