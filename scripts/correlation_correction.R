# deconvolution results show systematic under- or overprediction per cell type 
# try to correct this through:
# correct for linear regression intercept per cell type
# get noise score per cell type through lin reg residuals

library(ggplot2)
library(ggpubr)


# load a run
deconv <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_finalEven_tpm_ct0_rep0/deconvolution.rds")
metrics <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_finalEven_tpm_ct0_rep0/results_metric.rds")


ground_truth <- metrics$facs_groud_truth %>% as.data.table(keep.rownames = "cell_type")
ground_truth <- data.table::melt(ground_truth, id.vars = "cell_type",
                                 variable.name = "sample",
                                 value.name = "fraction")
deconv <- as.data.table(deconv, keep.rownames = "sample" )
predicted <- data.table::melt(deconv, id.vars = "sample",
                              variable.name = "cell_type",
                              value.name = "fraction")

fractions <- merge(ground_truth, predicted,
                   by = c("cell_type", "sample"), all = TRUE,
                   suffixes = c(".true", ".predicted"))


ggplot(fractions, aes(fraction.true, fraction.predicted))+
  geom_point(aes(color=cell_type))+
  stat_cor()+
  geom_abline(slope = 1)+
  #coord_equal(ratio = 1)+
  labs(title = "Predicted Fractions VS. True Fraction per Cell Type
       Bassez, DWLS, all cell types, 10000 cell/sample")

# perform linear regression for one celltype, get slope and intercept

# T Cell CD4 EX
linreg <- lm(fraction.predicted~fraction.true, data = fractions[cell_type == "TCELL_CD4_EX_Proliferating"])
coeff <- coefficients(linreg)
slope <- coeff[2]
intercept <- coeff[1]



ggplot(fractions[cell_type == "TCELL_CD4_EX_Proliferating"], aes(fraction.true, fraction.predicted))+
  geom_point(aes(color=cell_type))+
  stat_cor()+
  geom_abline(slope = 1)+
  stat_smooth(method = "lm", geom = "smooth")+
  #coord_equal(ratio = 1)+
  xlim(0, NA)+
  ylim(0, NA)
  labs(title = "Predicted Fractions VS. True Fraction per Cell Type
       Bassez, DWLS, all cell types, 10000 cell/sample")

adj_frac <- fractions
adj_frac
adj_frac[cell_type == "TCELL_CD4_EX_Proliferating", fraction.corrected := (fraction.predicted - intercept),]

fractions$diff <- fractions$fraction.true - fractions$fraction.predicted
ggplot(fractions, aes(sample, cell_type, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "True Fractions - \nPredicted Fractions")


celltypes <- fractions$cell_type %>% unique()
corrected_fractions <- lapply(celltypes, function(ctype){
  c_frac <- fractions[cell_type == ctype]
  lr <- lm(fraction.predicted ~ fraction.true, data = c_frac)
  coeff <- coefficients(lr)
  
  c_frac$slope <- coeff[2]
  c_frac$intercept <- coeff[1]
  mn_true <- c_frac$fraction.true %>% mean()
  mn_pred <- c_frac$fraction.predicted %>% mean()
  if (((coeff[1]>0) & ((mn_true - mn_pred) < 0)) | ((coeff[1]<0) & ((mn_true - mn_pred) > 0))) {
    c_frac[, fraction.corrected := fraction.predicted - abs(intercept),]
  } else {
    c_frac[, fraction.corrected := fraction.predicted + abs(intercept),]
  }

  return(c_frac)
})
adj_frac <- rbindlist(corrected_fractions)

adj_frac <- merge(adj_frac, pretty_names, by.x="cell_type", by.y="data_name")

adj_frac$cor_diff <- adj_frac$fraction.true - adj_frac$fraction.corrected
ggplot(adj_frac, aes(sample, cell_type, fill=cor_diff))+geom_tile()+
  geom_text(aes(label= round(cor_diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "True Fractions - \nCorrected Fractions")

# boxplot shows improvement per celltype
box_frac <- melt(adj_frac,
                 measure.vars = c("diff", "cor_diff"),
                 variable.name = "adjust",
                 value.name = "diff")
ggplot(box_frac, aes(plot_name, diff, color = adjust))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  geom_abline(slope=0)+
  xlab("Cell Type")+
  ylab("D")+
  scale_color_manual("Value", values = c("#619CFF", "#F8766D"), labels = c("Original", "Adjusted"))


# scatter plot shows predicted and corrected
sca_frac <- melt(adj_frac,
                 measure.vars = c("fraction.predicted", "fraction.corrected"),
                 variable.name = "adjusted",
                 value.name = "fraction.predicted")
ggplot(sca_frac[cell_type=="TCELL_CD4_N"], aes(fraction.true, fraction.predicted, color=adjusted))+
  geom_point()+
  stat_cor()+
  geom_abline(slope = 1)+
  stat_smooth(method = "lm", geom = "smooth")+
  xlim(0, NA)+
  ylim(0, NA)

# improvement

adj_frac[, improvement := abs(diff) - abs(cor_diff)]
ggplot(adj_frac, aes(fraction.true, improvement))+
  geom_point()+
  stat_cor()+
  geom_smooth()


smoo_frac <- adj_frac[, mean(improvement), by=cell_type]
colnames(smoo_frac) <- c("cell_type", "mean_improvement")
ggplot(adj_frac, aes(slope, improvement))+
  geom_smooth()+
  geom_point()

ggplot(smoo_frac, aes(cell_type, mean_improvement))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  