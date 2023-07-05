################################################################################
#                   genotypes effect on drug effectiveness                     #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(jtools)
library(biomaRt)
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response"
setwd(project.dir)
################################################################################
mph.spark <- read_tsv("data/derivatives/spark-rm-mph-effect.tsv")
genotypes <- fread(file = "data/derivatives/spark-mph-samples-genotypes.xmat.gz", header = T, nThread = 4)
genotypes <- genotypes[-1,-1]
# genotypes.5 <- fread(file = "data/derivatives/spark-mph-samples-genotypes-LC-merged-chr5.xmat.gz", header = T, nThread = 4)
# genotypes.5 <- genotypes.5[-1,-1]
# genotypes.16 <- fread(file = "data/derivatives/spark-mph-samples-genotypes-LC-merged-chr16.xmat.gz", header = T, nThread = 4)
# genotypes.16 <- genotypes.16[-1,-1]

gc()
geno.pcs <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/PCA/all/PCs.tsv") %>%
  filter(IID %in% mph.spark$IID)
rownames(geno.pcs) = geno.pcs$IID

mph_genes <- c("SLC6A3","SLC6A2","HTR1A")
gene.info <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/genomics/bioMart-genes-lengths-ENSG-symbol.rds") %>%
  filter(hgnc_symbol %in% mph_genes)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl <- getBM(attributes=c('ensembl_gene_id','chromosome_name', 'start_position', 'end_position', 'strand'),
      filters=c('ensembl_gene_id'),
      values=list(gene.info$ensembl_gene_id),
      mart=ensembl)
gene.info <- left_join(gene.info, ensembl)

genotyped.variants <- data.frame(variant = colnames(genotypes)) %>%
  filter(variant != "IID") %>%
  mutate(pos = sub("\\chr[0-9]+:", "", variant)) %>%
  mutate(pos = as.integer(sub("\\:.*", "", pos))) %>%
  mutate(chr = as.integer(sub("chr", "", sub("\\:.*", "", variant)))) %>%
  filter(chr %in% c(5,16)) %>%
  mutate(gene1= ifelse(chr==gene.info$chromosome_name[1] & pos>gene.info$start_position[1] & gene.info$end_position[1]>pos, T, F)) %>%
  mutate(gene2= ifelse(chr==gene.info$chromosome_name[2] & pos>gene.info$start_position[2] & gene.info$end_position[2]>pos, T, F)) %>%
  mutate(gene3= ifelse(chr==gene.info$chromosome_name[3] & pos>gene.info$start_position[3] & gene.info$end_position[3]>pos, T, F)) %>%
  filter(gene1 == T | gene2 ==T | gene3==T)
genes.genotypes <- genotypes %>% as.data.frame() %>% dplyr::select(c("IID", genotyped.variants$variant))
gc()
###############################################################################
# look at genotypes of these genes and their correlation with drug effectiveness
x <- genes.genotypes %>% column_to_rownames("IID")%>%mutate_all(function(x) as.numeric(x))
y <- mph.spark %>% dplyr::select(IID, mph_effect)%>% column_to_rownames("IID")%>% as.matrix()
x <- x[rownames(y),]%>% as.matrix()
all(rownames(x)==rownames(y))
dim(y)
dim(x)

jtools::export_summs(lm(y~x))
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################

################################################################################
