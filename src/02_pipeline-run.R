# run the pipeline with different variations
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response"
setwd(project.dir)
source("src/01_drug_response-functions.R")


dataset <- "spark"
# dataset <- "abcd"
dry.run = F

for (weights_source in c("tissue", "celltype")) {
  if (weights_source == "tissue") {
    if (dataset == "spark") {
      genotypes.path <- "data/derivatives/spark-mph-samples-genotypes-LC-merged-Brain_Frontal_Cortex_BA9.xmat.gz"
    }else {
      genotypes.path <- "data/derivatives/abcd-mph-samples-genotypes-TT-Brain_Frontal_Cortex_BA9.xmat.gz"
    }
  }else {
    if (dataset == "spark") {
      genotypes.path <- "data/derivatives/spark-mph-samples-genotypes-LC-merged-Excitatory-FDR-sig.xmat.gz"
    }else {
      genotypes.path <- "data/derivatives/abcd-mph-samples-genotypes-TT-Excitatory-FDR-sig.xmat.gz"
    }
  }
  for (genes in c("all", "targets")) {
    # print(genes)
    if (genes == "targets" && weights_source == "celltype") {
      break
    }
    for (genetic_correct in c(T,F)) {
      for (scaling in c(T,F)) {
        for (response_method in c(1,2)) {
          if (dry.run == T) {
            print(paste0("data/derivatives/m-outputs/", dataset, "/model-", weights_source,
                   "-", genes,"-", genetic_correct, "-", scaling, "-", response_method, 
                   ".rds"))
          }else {
            m.output <- Go.BP(from = weights_source, genes_set = genes, correct = genetic_correct, 
                              cores = 18, scale = scaling, method = 1, 
                              genotypes_path = genotypes.path)
            write_rds(m.output, paste0("data/derivatives/m-outputs/", dataset, "/model-", weights_source,
                                       "-", genes,"-", genetic_correct, "-", scaling, "-", response_method, 
                                       ".rds"), compress = "gz")
          }
        }
      }
    }
  }
}


# if you only want the imputed tx:
genotypes_path <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/data/derivatives/genotypes/genotypes-TT-Brain_Frontal_Cortex_BA9.xmat.gz"
# genotypes_path <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/data/derivatives/genotypes/genotypes-TT-Excitatory-FDR-sig.xmat.gz"
from == "tissue"
# from == "celltype"
  
if (from == "tissue") {
  tissue <- "Brain_Frontal_Cortex_BA9"
  genotypes <- fread(file = genotypes_path, header = T, nThread = cores)
  genotypes <- genotypes[-1,-1]
  gc()
  print(paste0("Done with: ", "reading genotypes file for tissue"))
  # get tissue weights
  tissue.weights <- read.table(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/tx-imputation/UTMOST-GTEx-model-weights/tmp/rsid-for-", tissue), row.names = 1)
  ready.weights <- tissue.weights %>% dplyr::select(variant = ID_02_UTMOST, gene, weight) %>%
    distinct(variant, gene, .keep_all = T) %>%
    filter(variant %in% colnames(genotypes))
  gc()
  print(paste0("Done with: ", "reading weights file for tissue"))
} else if (from == "celltype") {
  genotypes <- fread(file = genotypes_path, header = T, nThread = cores)
  genotypes <- genotypes[-1,-1]
  gc()
  print(paste0("Done with: ", "reading genotypes file for celltype"))
  # get celltype weights
  celltype <- "Excitatory"
  celltype.weights <- read_rds(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/", celltype,"-weights-fdr-sig.rds"))
  ready.weights <- celltype.weights %>%
    filter(FDR<0.05) %>%
    dplyr::select(variant=ID_37, gene, weight=beta) %>%
    distinct(variant, gene, .keep_all = T) %>%
    filter(variant %in% colnames(genotypes))
  rm(celltype.weights)
  gc()
  print(paste0("Done with: ", "reading weights file for celltype"))
}
# impute tx for all genes
imputed.tx <- impute.tx(genotypes = genotypes, 
                        weights = ready.weights, threads = cores) %>% as.data.frame()
gc()

################################################################################
# supp figures
celltype.weights <- read_rds(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/", celltype,"-weights-fdr-sig.rds")) %>%
  filter(FDR<0.05) %>%
  dplyr::select(variant=ID_37, gene, weight=beta) %>%
  distinct(variant, gene, .keep_all = T) %>%
  mutate(chr = as.numeric(sub("chr", "", sub(":.*", "", variant)))) %>%
  drop_na()
p1<- celltype.weights %>%
  ggplot(aes(x=weight))+
  geom_histogram()
p2<- celltype.weights %>%
  mutate(chr = as.factor(chr))%>%
  ggplot(aes(x=chr))+
  geom_bar()+
  labs(caption = paste0("number of significant (FDR<0.05) eQTL-gene association per chromosome"))
p3<- celltype.weights %>%
  mutate(chr = as.factor(chr))%>%
  group_by(chr, gene) %>%
  select(chr, gene) %>%
  mutate(count = n()) %>%
  distinct(chr, gene, .keep_all = T)%>%
  ggplot(aes(x=chr, y=count))+
  geom_boxplot()+
  labs(caption = paste0("number of significant (FDR<0.05) eQTL association per gene per chromosome"))
p4<- celltype.weights %>%
  mutate(chr = as.factor(chr))%>%
  distinct(chr, gene)%>%
  ggplot(aes(x=chr))+
  geom_bar()+
  labs(caption = paste0("number of genes with eQTL weights per chromosome"))
patchwork::wrap_plots(p1,p2,p3,p4, ncol = 2)
################################################################################