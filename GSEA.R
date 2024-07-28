rm(list = ls());gc()
library(dplyr)
library(readr)
library(stringr)
library(limma)
library(clusterProfiler)
library(tidyverse)
library(enrichplot)
library(patchwork)
library(openxlsx)

# load data -----------
read_gct <- function(file) {
  lines <- readLines(file)
  data <- read.delim(text = lines[-(1:2)], header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  return(data)
}
exp <- read_gct("Desktop/BxPC_expression_matrix.gct") %>% select(-'DESCRIPTION') %>% mutate(across(where(is.numeric), ~ round(.x, 5))) # [1] 17941     8
exp <- log2(exp+1)
group <- str_split_fixed(colnames(exp),'_',3)[,2]
gs <- qusage::read.gmt('Desktop/PDAC_subtypes.gmt') %>% 
  enframe(name = "name", value = "string") %>%
  unnest(string)


# DE analysis -----------
myCompare <- c("Sph_Par")
pd <- data.frame(group)
pd$contrast <- as.factor(pd$group)
design <- model.matrix(~ 0 + contrast , data = pd)
fit <- lmFit(exp, design)
contrast  <- makeContrasts(
  Sph_Par = contrastSph - contrastPar,
  levels = design)
fits <- contrasts.fit(fit, contrast)
ebFit <- eBayes(fits)

if(T){
limma.res <- topTable(ebFit, coef = 1, adjust.method = 'fdr', number = Inf)
limma.res <- limma.res %>% filter(!is.na(adj.P.Val)) %>%
    mutate( logP = -log10(P.Value) ) %>%
    mutate( group = "not sig") %>%
    mutate( tag = myCompare ) %>%
    mutate( Gene = rownames(limma.res))
limma.res$group[which( (limma.res$P.Value < 0.05) & (limma.res$logFC > 0) )] = "up"
limma.res$group[which( (limma.res$P.Value < 0.05) & (limma.res$logFC < 0) )] = "down"
caseName <- str_split_i(myCompare,"_",1)
ctrlName <- str_split_i(myCompare,"_",2)
id_case <- which(pd$group==caseName)
id_control <- which(pd$group==ctrlName)
AVG_all <-  rowMeans(exp)
AVG_case <- rowMeans(exp[,id_case])
AVG_control <- rowMeans(exp[,id_control])
AVG <- cbind(AVG_all,AVG_case,AVG_control)
AVG <- AVG %>% as.data.frame() %>% mutate(Gene=rownames(AVG))
limma.res <- limma.res %>% left_join(AVG,by=c("Gene"="Gene"))
limma.res <- limma.res %>% dplyr::select(Gene,AVG_all,AVG_case,AVG_control,everything()) %>%
  mutate(cutoff="p<0.05") %>% 
  filter(!group=='not sig') %>% 
  arrange(desc(logFC))
}



# DE analysis -----------
gene_df <- limma.res %>% mutate(rank = rank(logFC,  ties.method = "random")) %>%
  arrange(desc(rank))

genelist <- gene_df$logFC
names(genelist) <- gene_df$Gene 

gsea_res <- GSEA(genelist, 
                 TERM2GENE = gs,
                 verbose = T
)

gsea_dat <- as.data.frame(gsea_res)
write.xlsx(gsea_dat,"Desktop/GSEA.xlsx",overwrite = T,colNames=T)

for (i in 1:3) {
  # i <- 4
  pdf(paste0("Desktop/gsea_plot_", i, ".pdf"), width = 15, height = 15)
  gseaplot2(gsea_res,
            geneSetID = rownames(gsea_dat)[i],
            title = "",
            color = "green",
            base_size = 12,
            rel_heights = c(1.5, 0.5, 1),
            subplots = 1:3,
            pvalue_table = T,
            ES_geom = "line")
  dev.off()
}


# sc

a <- readRDS('../Case18 PanIN-sc/Result/02-PDAC_subtype/4-Final_celltype/sce.all.int.rds')
FeatureStatPlot(
  srt = a, bg.by = "celltype", group.by = "celltype",
  stat.by = 'GPR110', add_box = TRUE
)
ggsave(filename="FeatureStatPlot.pdf",width = 18,height = 10)

P1 <- FeatureDimPlot(
  srt = a, features = c("GPR110"),
  reduction = "TSNE", theme_use = "theme_blank"
)

P2 <- CellDimPlot(srt = a, group.by = c('celltype'),reduction = "tsne",theme_use = "theme_blank")

P <- P2+P1
ggsave(filename="final_umap.pdf",width = 12,height = 6)

