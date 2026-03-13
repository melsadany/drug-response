celltype="Excitatory"
celltype.weights <- read_rds(correct_path(paste0("/Dedicated/jmichaelson-wdata/msmuhammad/data/celltypes-cis-eQTLs/data/derivatives/", celltype,"-weights-fdr-sig.rds"))) %>%
  dplyr::select(variant=ID_37, gene, weight=beta) %>%
  distinct(variant, gene, .keep_all = T) %>%
  mutate(chr = as.numeric(sub("chr", "", sub(":.*", "", variant)))) %>%
  drop_na()
p1<- celltype.weights %>%
  ggplot(aes(x=weight))+
  geom_histogram(fill=palette.1[2]) + bw.theme + 
  # scale_y_log10()+
  labs(x=paste0("eQTL beta"))
p2<- celltype.weights %>%
  mutate(chr = as.factor(chr))%>%
  ggplot(aes(x=chr,fill=chr))+
  geom_bar(show.legend = F)+bw.theme+scale_fill_manual(values = c(palette.1,palette.2,six.colors))+
  scale_y_log10()+
  labs(x="Chromosome",y="Count of significant associations")
p3<- celltype.weights %>%
  mutate(chr = as.factor(chr))%>%
  group_by(chr, gene) %>%
  select(chr, gene) %>%
  mutate(count = n()) %>%
  distinct(chr, gene, .keep_all = T)%>%
  ggplot(aes(x=chr, y=count, fill=chr))+
  geom_boxplot(show.legend = F)+bw.theme + scale_fill_manual(values = c(palette.1,palette.2,six.colors))+
  scale_y_log10()+
  labs(x="Chromosome", y = "Count of significant cis-eQTLs\n(per gene)")
p3
p4<- celltype.weights %>%
  mutate(chr = as.factor(chr))%>%
  distinct(chr, gene)%>%
  ggplot(aes(x=chr,fill=chr))+
  geom_bar(show.legend = F)+bw.theme+scale_fill_manual(values = c(palette.1,palette.2,six.colors))+
  labs(y = ("Count of genes with significant eQTLs"),x="Chromosome")
p4
pdf("figs/excitatory-eQTL-supp.pdf",6,12)
patchwork::wrap_plots(p1,p2,p3,p4, ncol = 1)
dev.off()
ggsave2("figs/excitatory-eQTL-supp.pdf",6,12)
