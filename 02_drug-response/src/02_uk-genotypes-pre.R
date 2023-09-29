################################################################################
#                                 genomics prework                             #
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
# samples of interest
# I decided to only keep the ones with clonidine status OR have some value in trailMaking_duration_to_complete_numeric_path

################################################################################
# section to be done once
# keep variants of interest for drug response
celltype.weights <- read_rds(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/Inhibitory-weights.rds"))
celltype.weights.2 <- celltype.weights %>%
  mutate(FDR = p.adjust(nominal_pval, method = "fdr")) %>%
  filter(FDR<0.05) %>%
  distinct(ID_37, gene, .keep_all = T) %>%
  mutate(CHR = as.numeric(sub("chr", "", sub(":.*", "", ID_37)))) %>%
  mutate(POS = sub(".*:(\\d+):.*", "\\1", ID_37)) %>%
  mutate(REF = sub(".*:(.):.*", "\\1", ID_37)) %>%
  mutate(ALT = sub(".*:([^:]+)$", "\\1",ID_37)) %>%
  mutate(uk_37 = paste0(CHR, ":", POS, "_", REF, "_", ALT)) %>%
  drop_na()
write_lines(celltype.weights.2$uk_37, "/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/Inhibitory-UK_37-FDR-sig")

uk.plink.dir <- paste0(ukb.dir, "/../genotypes/imputed_v3/plink")

registerDoMC(cores = 10)
foreach (i = c(12:16)) %dopar% {
  print(i)
  gcta_command <- paste0(
    "/Dedicated/jmichaelson-wdata/trthomas/bin/gcta/gcta_1.93.2beta/gcta64", # GCTA executable
    " --bfile ", "ukb_imp_chr", i, "_v3", 
    " --keep ", "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response/data/derivatives/uk-sample.fam", 
    " --extract ", paste0("/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/Inhibitory-UK_37-FDR-sig"), #list of SNPs
    " --thread-num ", 10, # threads
    " --recode ",
    " --out ", paste0(project.dir, "/data/derivatives/uk-genotypes-sample/uk-Inhibitory-FDR-sig-chr", i))
  print(gcta_command)
  system(gcta_command)
}
# ################################################################################
# # combining genotypes matrices for all chromosomes
# base.chr.plink <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response/data/derivatives/uk-genotypes-sample/uk-Inhibitory-FDR-sig-chr"
# all.m <- read_delim("data/derivatives/uk-sample.fam", col_names = NULL) %>%
#   select(IID="X1")
# for (i in 1:22) {
#   # i=22
#   # all.m[[i]] <- fread(paste0(base.chr.plink, 22 ,".xmat.gz"), nThread = 8)[-1,-1]
#   chr.m <- fread(paste0(base.chr.plink, 22 ,".xmat.gz"), nThread = 8)[-1,-1]
#   all.m <- full_join(all.m, chr.m %>% mutate(IID = as.numeric(IID)))
#   rm(chr.m)
#   gc()
# }
# 
# ################################################################################
# run your pipeline
source("../src/01_drug_response-functions.R")
cores = 8
genes_set = "all"
approach = 1
part <- read_delim("/Dedicated/jmichaelson-wdata/msmuhammad/projects/drug-response/02_drug-response/data/derivatives/uk-sample.fam", col_names = NULL) %>%
  select(X2) %>% mutate(X2=as.character(X2))

registerDoMC(cores = 8)
imputed.tx.all <- foreach (i = 1:22, .combine = cbind) %dopar% {
  # i=22
  genotypes.path <- paste0("data/derivatives/uk-genotypes-sample/uk-Inhibitory-FDR-sig-chr", i, ".xmat.gz")
  
  # section 1
  genotypes <- fread(file = genotypes.path, header = T, nThread = cores)
  genotypes <- genotypes[-1,-1] %>% 
    mutate_at(.vars = vars(colnames(genotypes)[-c(1,2)]), .funs = function(x) ifelse(is.na(x), median(as.numeric(x), na.rm = T), as.numeric(x)))
  gc()
  print(paste0("Done with: ", "reading genotypes file for celltype"))
  # get celltype weights
  celltype <- "Inhibitory"
  celltype.weights <- read_rds(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/", celltype,"-weights-fdr-sig.rds"))
  ready.weights <- celltype.weights %>%
    filter(FDR<0.05) %>%
    dplyr::select(variant=uk_37, gene, weight=beta) %>%
    distinct(variant, gene, .keep_all = T) %>%
    filter(variant %in% colnames(genotypes))
  rm(celltype.weights)
  gc()
  print(paste0("Done with: ", "reading weights file for celltype"))
  
  # section 2
  imputed.tx <- impute.tx(genotypes = genotypes, 
                          weights = ready.weights, threads = cores) %>% as.data.frame()
  gc()
  print(paste0("Done with: ", "imputing tx"))
  if (is.null(imputed.tx)) {
    print("exit")
    return(NULL)
  }
  imputed.tx <- imputed.tx[part$X2,]
  return(imputed.tx)
}
write_rds(imputed.tx.all, "data/derivatives/uk-clonidine-sample-inhibit-imputed-tx.rds", compress = "gz")
##########
# assume you have a vector of imputed tx for all genes
imputed.tx.all <- read_rds("data/derivatives/uk-clonidine-sample-inhibit-imputed-tx.rds")
# section 3
cmap.drug.sig <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/LINCS/cmap.of.int_v2.rds")
drug.sig <- cmap.drug.sig["clonidine",]%>%as.data.frame()%>%rename(clonidine=1)
print(paste0("Done with: ", "reading clonidine signature"))

m <- predict.response(samples.tx = if (ncol(imputed.tx.all)>1) {t(scale(imputed.tx.all))} else {t(imputed.tx.all)}, 
                      drug.sig = scale(drug.sig), approach = approach, threads = cores, set = genes_set) %>%
  as.data.frame() %>%
  rownames_to_column("IID") %>%
  rename(m = 2)
print(paste0("Done with: ", "predicting drug response with scaling option and approach 1"))


write_rds(m, paste0("data/derivatives/uk-clonidine-model-3.rds"), compress = "gz")
################################################################################