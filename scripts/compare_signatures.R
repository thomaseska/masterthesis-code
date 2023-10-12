# compare DWLS signature to NSForest result
library(data.table)
library(dplyr)
library(tidyr)
# load
dwls <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_finalEven_tpm_ct0_rep0/signature.rds")
nsf <- fread("/nfs/data/tcell_deconvolution/data/bassez/ns_forest_result.csv")
nsf$V1 <- NULL

nsf$NSForest_markers <- nsf[, gsub("\\[|\'|\\]", "", NSForest_markers)]
nsf <- nsf %>% separate_rows(NSForest_markers, sep = ",\\s*") %>% as.data.table()

dwls <- as.data.table(dwls, keep.rownames = "Gene")
dwls <- melt(dwls, id.vars = "Gene", variable.name = "cell_type", value.name = "value")

nsf_with_scores <- merge(nsf, dwls, 
                         by.x = c("clusterName", "NSForest_markers"),
                         by.y = c("cell_type", "Gene"),
                         all.x = TRUE,
                         all.y = FALSE)
ggplot(dwls[ value <100], aes(value))+geom_density()+geom_density(data = nsf_with_scores, aes(value), color ="red")

nsf_with_scores_ext <- merge(nsf, dwls, 
                             by.x = c("NSForest_markers"),
                             by.y = c( "Gene"),
                             all.x = TRUE,
                             all.y = FALSE)
nsf_with_scores_ext <- nsf_with_scores_ext[, c(1, 2, 12, 13)]
nsf_with_scores_ext[,original_cluster:= (clusterName==cell_type),]
nsf_with_scores_ext <- nsf_with_scores_ext[order(clusterName, cell_type, NSForest_markers)]

nsf_with_scores_ext$cell_type <- nsf_with_scores_ext$cell_type %>% as.character()

nsf_with_scores_ext <- merge(nsf_with_scores_ext, pretty_names, by.x="cell_type", by.y = "data_name")
nsf_with_scores_ext <- merge(nsf_with_scores_ext, pretty_names, by.x="clusterName", by.y = "data_name")

setorder(nsf_with_scores_ext, plot_name.y, plot_name.x, NSForest_markers, na.last = TRUE)

geneorder <- nsf_with_scores_ext$NSForest_markers %>% unique()
typeorder <- nsf_with_scores_ext$clusterName %>% unique()



ggplot(nsf_with_scores_ext, aes(factor(cell_type, levels = typeorder), factor(NSForest_markers, levels = geneorder), fill = value))+
  geom_tile()+
  geom_tile(data=nsf_with_scores_ext, aes(color = original_cluster, size = original_cluster), alpha=0)+
  scale_colour_manual("original", values = c("white", "red", "black"))+ 
  scale_size_manual("original", values = c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

nsf_with_scores_ext[, zscores:=zscore(value), by=NSForest_markers]
ggplot(nsf_with_scores_ext, aes(plot_name.x, factor(NSForest_markers, levels = geneorder), fill = zscores))+
  geom_tile()+
  geom_tile(data=nsf_with_scores_ext, aes(color = original_cluster, size = original_cluster), alpha=0)+
  scale_colour_manual("NS-Forest \nMarker", values = c("white", "red", "black"))+ 
  scale_size_manual("NS-Forest \nMarker", values = c(0, 0.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab("Cell Type")+
  ylab("Marker Gene")+
  scale_fill_continuous("Z-score")
