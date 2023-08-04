################################################################################
#                correlation between mph transcriptomics and PGS genes         #
################################################################################
# packages setup ----------------------------------------------------------
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(ggpubr);library(ggh4x)
####
# project dir setup -------------------------------------------------------
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response"
setwd(project.dir)
####
# load the CMAP signatures and gFactor FUMA  ------------------------------
cmap <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int_v2.rds")[c("methylphenidate",
                                                                                           "atomoxetine",
                                                                                           "guanfacine",
                                                                                           "clonidine",
                                                                                           "ibuprofen"),] %>% 
  t() %>% as.data.frame() %>%
  mutate_all(.funs = function(x) scale(x, scale = T, center = T)[,1]) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -c(gene), names_to = "drug", values_to = "exp")
  
gf <- read_table("data/ukbb-gFactor-FUMA/genes.txt")
# combine both together
cmap.gf <- full_join(cmap %>% 
                        mutate(gf = ifelse(gene %in% gf$symbol, T, F)), 
                      gf %>% 
                        rename(gene = symbol)%>%
                        select(-c(3:9))) %>%
  drop_na(exp, pLI)
# load ADHD FUMA
ad <- read_table("data/Demontis-ADHD-FUMA/genes.txt")
# combine both together
cmap.ad <- full_join(cmap %>% 
                       mutate(ad = ifelse(gene %in% ad$symbol, T, F)), 
                     ad %>% 
                       rename(gene = symbol)%>%
                       select(-c(3:9))) %>%
  drop_na(exp, pLI)
####
# correlation between gene exp in CMAP and FUMA measures ------------------
p1 <- cmap.gf %>%
  # pivot_longer(cols = c(pLI, ncRVIS, starts_with("posMap")), names_to = "fuma", values_to = "val") %>%
  pivot_longer(cols = c(pLI), names_to = "fuma", values_to = "val") %>%
  ggplot(aes(x=exp, y=val))+
  geom_point(size=0.3, na.rm = T)+
  geom_smooth(method = "glm", na.rm = T) +
  stat_cor(method = "spearman", na.rm = T) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_grid2(rows = vars(fuma), cols = vars(drug), scales = "free", independent = T)+
  labs(x="gene expression in CMAP drug signature",
       caption = paste0("n(genes): ", nrow(cmap.gf%>%drop_na(exp,pLI)%>%distinct(gene)), "\n",
                        "FUMA variables description: \n",
                        "\tpLI : pLI score from ExAC database. The probability of being loss-of-function intolerant. The higher the score is, the more intolerant to loss-of-function mutations the gene is. \n",
                        "\tncRVIS : Non-coding residual variation intolerance score. The higher the score is, the more intolerant to noncoding variants the gene is.\n",
                        "\tposMapSNPs (posMap): The number of SNPs mapped to gene based on positional mapping (after functional filtering if parameters are given).\n",
                        "\tposMapMaxCADD (posMap): The maximum CADD score of mapped SNPs by positional mapping."))
p2 <- corr.table(cmap.gf%>%distinct(gene,pLI)%>%select(pLI), 
                 cmap.gf %>% pivot_wider(names_from = "drug", values_from = "exp", id_cols = "gene")%>%select(-1),
                 method = "spearman") %>%
  filter(V1 == "pLI", V2 != V1) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = paste0(rho, ": ", round(r,3), "\np: ",round(pval, 3))))+
  geom_tile()+geom_text(size=3)+
  my.guides+redblu.col.gradient+null_labs
p3 <- cmap.ad %>%
  # pivot_longer(cols = c(pLI, ncRVIS, starts_with("posMap")), names_to = "fuma", values_to = "val") %>%
  pivot_longer(cols = c(pLI), names_to = "fuma", values_to = "val") %>%
  ggplot(aes(x=exp, y=val))+
  geom_point(size=0.3, na.rm = T)+
  geom_smooth(method = "glm", na.rm = T) +
  stat_cor(method = "spearman", na.rm = T) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_grid2(rows = vars(fuma), cols = vars(drug), scales = "free", independent = T)+
  labs(x="gene expression in CMAP drug signature",
       caption = paste0("n(genes): ", nrow(cmap.ad%>%drop_na(exp,pLI)%>%distinct(gene))))
p4 <- corr.table(cmap.ad%>%distinct(gene,pLI)%>%select(pLI), 
                 cmap.ad %>% pivot_wider(names_from = "drug", values_from = "exp", id_cols = "gene")%>%select(-1),
                 method = "spearman") %>%
  filter(V1 == "pLI", V2 != V1) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = paste0(rho, ": ", round(r,3), "\np: ",round(pval, 3))))+
  geom_tile()+geom_text(size=3)+
  my.guides+redblu.col.gradient+null_labs
patchwork::wrap_plots(p3+labs(title = "ADHD"),p4,p1+labs(title = "gFactor"),p2, widths = c(3,0.6))
####