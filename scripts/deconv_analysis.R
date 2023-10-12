# cleaned up script for deconvolution result analysis

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

zscore <- function(vec) {
  return((vec - mean(vec))/sd(vec))
}

dscore <- function(diff, truth) {
  return(diff / (truth + 1e-4))
}

sqdval <- function(d) {
  if (d >= 0) {
    return(sqrt(d))
  } else {
    return(-sqrt(abs(d)))
  }
}

sqdscore <- function(dvec) {
  sqd <- sapply(dvec, sqdval)
  return(sqd)
}

revsqdval <- function(sqd) {
  if (sqd >= 0) {
    return(sqd^2)
  } else {
    return(-(sqd^2))
  }
}
revsqdscore <- function(sqdvec) {
  dvec <- sapply(sqdvec, revsqdval)
  return(dvec)
}

create_frac_tbl <- function(deconv, metrics, missing = "none", scenario = "even") {
  
  ground_truth <- metrics$facs_groud_truth %>% as.data.table(keep.rownames = "cell_type")
  ground_truth <- data.table::melt(ground_truth, id.vars = "cell_type",
                                   variable.name = "sample",
                                   value.name = "fraction")
  if (missing != "none") {
    to_add <- data.table("cell_type"=missing,
                         "sample"=paste0(scenario, "_sample", 1:10),
                         "fraction"=0)
    ground_truth <- rbind(ground_truth, to_add)
  }
  
  deconv <- as.data.table(deconv, keep.rownames = "sample" )
  predicted <- data.table::melt(deconv, id.vars = "sample",
                                variable.name = "cell_type",
                                value.name = "fraction")
  
  fractions <- merge(ground_truth, predicted,
                     by = c("cell_type", "sample"), all = TRUE,
                     suffixes = c(".true", ".predicted"))
  
  fractions$diff <- fractions$fraction.true - fractions$fraction.predicted
  
  fractions[, diff.d:=dscore(diff, fraction.true), by= cell_type]
  
  fractions <- merge(fractions, pretty_names, by.x="cell_type", by.y="data_name", all.x=TRUE)
  
  return(fractions)
}

# use for bayesprism
deconv <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/bayesprism_finalBayes_counts_finalEven_counts_ct0_rep0/deconvolution.rds")
metrics <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/bayesprism_finalBayes_counts_finalEven_counts_ct0_rep0/results_metric.rds")

deconv_m <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/bayesprism_finalBayes_counts_finalMirror_counts_ct0_rep0/deconvolution.rds")
metrics_m <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/bayesprism_finalBayes_counts_finalMirror_counts_ct0_rep0/results_metric.rds")


# use for dwls
deconv <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_finalEven_tpm_ct0_rep0/deconvolution.rds")
metrics <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_finalEven_tpm_ct0_rep0/results_metric.rds")

deconv_m <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_finalMirror_tpm_ct0_rep0/deconvolution.rds")
metrics_m <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_finalMirror_tpm_ct0_rep0/results_metric.rds")


frac_dt <- create_frac_tbl(deconv, metrics)
frac_m_dt <- create_frac_tbl(deconv_m, metrics_m)


frac_dt$scenario <- "Even Fractions"
frac_m_dt$scenario <- "Mirrored Fractions"

frac_dt <- rbind(frac_dt, frac_m_dt)



frac_dt <- merge(frac_dt, sample_lookup, by= "sample")

# basic scatterplot
ggplot(frac_m_dt, aes(fraction.true, fraction.predicted))+
  geom_point()+
  stat_cor(method = "spearman")+
  geom_abline(slope = 1)

p1 <- ggplot(frac_dt[scenario == "Even Fractions"], aes(fraction.true, fraction.predicted))+
  geom_point()+
  stat_cor(method = "spearman")+
  geom_abline(slope=1)+
  # facet_wrap(vars(scenario), scales = "free")+
  xlab("True Fraction")+
  ylab("Predicted Fraction")

#scatterplot with regressionlines per cell type
ggplot(frac_dt[scenario=="Mirrored Fractions"], aes(fraction.true, fraction.predicted))+
  geom_point(aes(color=plot_name))+
  stat_cor(method = "spearman")+
  geom_abline(slope = 1, linewidth=1.5)+
  stat_smooth(method = "lm", geom = "smooth", aes(color = plot_name))+
  xlab("True Fraction")+
  ylab("Predicted Fraction")+
  scale_color_discrete(name = "Cell Type")+
  facet_wrap(vars(type), scales ="free")

# scatterplot true frac x D 
p1 <- ggplot(frac_dt[scenario=="Mirrored Fractions"], aes(fraction.true, abs(diff)))+
  geom_point()+
  stat_cor(method = "spearman")+
  xlab("True Fraction")+
  ylab("Absolute Difference")+
  labs(title = "Difference Score")+
  theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(frac_dt[scenario=="Mirrored Fractions"], aes(fraction.true, abs(diff.d)))+
  geom_point()+
  stat_cor(method = "spearman")+
  xlab("True Fraction")+
  ylab("Absolute Normalized Difference")+
  labs(title = "Normalized Difference Score")+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(p1, p2, labels = "AUTO")

# heatmap
ggplot(frac_dt[scenario == "Even Fractions"], aes(sample_new, plot_name, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "D")+
  xlab("Sample")+
  ylab("Cell Type")

# heatmap mirror
ggplot(frac_dt[scenario == "Mirrored Fractions"], aes(sample_new, plot_name, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "D")+
  xlab("Sample")+
  ylab("Cell Type")

# heatmap norm
scale_breaks <- c(0, -sqrt(1), -sqrt(10), -sqrt(25))

ggplot(frac_dt[scenario == "Mirrored Fractions"], aes(sample_new, plot_name, fill=sqdscore(diff.d)))+geom_tile()+
  geom_text(aes(label= round(sqdscore(diff.d), 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "ND")+
  xlab("Sample")+
  ylab("Cell Type")

# boxplot
ggplot(frac_dt[scenario == "Even Fractions"], aes(plot_name, diff))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab("Cell Type")+
  ylab("D")+
  geom_abline(slope = 0)

# boxplot norm 
ggplot(frac_dt[scenario == "Mirrored Fractions"], aes(plot_name, sqdscore(diff.d)))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab("Cell Type")+
  ylab("ND")+
  scale_y_continuous(breaks = scale_breaks, labels = revsqdscore(scale_breaks))+
  geom_abline(slope = 0)


p1 <- ggplot(frac_dt[scenario == "Even Fractions"], aes(fraction.true, fraction.predicted))+
  geom_point()+
  stat_cor(method = "spearman")+
  geom_abline(slope=1)+
  # facet_wrap(vars(scenario), scales = "free")+
  xlab("True Fraction")+
  ylab("Predicted Fraction")


frac_e_corrs <- frac_dt[scenario=="Even Fractions",
                        SpearmanRho(x = fraction.true, y = fraction.predicted, conf.level = 0.95),
                        by = plot_name]
frac_e_corrs$conf.bound <- rep(c("rho", "conf.lower", "conf.upper"), 23)
frac_e_corrs <- dcast(frac_e_corrs, formula = ... ~ conf.bound, value.var = "V1")

p2 <- ggplot(frac_e_corrs, aes(plot_name, rho))+geom_col()+
  geom_errorbar(aes(ymin = conf.lower, ymax=conf.upper))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab("Cell Type")+
  ylab("Spearman Correlation")




ggarrange(p1, p2, labels = "AUTO")


frac_m_corrs <- frac_dt[scenario=="Mirrored Fractions",
                        SpearmanRho(x = fraction.true, y = fraction.predicted, conf.level = 0.95),
                        by = plot_name]
frac_m_corrs$conf.bound <- rep(c("rho", "conf.lower", "conf.upper"), 23)
frac_m_corrs <- dcast(frac_m_corrs, formula = ... ~ conf.bound, value.var = "V1")

ggplot(frac_m_corrs, aes(plot_name, rho))+geom_col()+
  geom_errorbar(aes(ymin = conf.lower, ymax=conf.upper))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab("Cell Type")+
  ylab("Spearman Correlation")

p1 <- ggplot(frac_dt[scenario == "Mirrored Fractions"], aes(fraction.true, fraction.predicted))+
  geom_point()+
  stat_cor(method = "spearman")+
  geom_abline(slope=1)+
  # facet_wrap(vars(scenario), scales = "free")+
  xlab("True Fraction")+
  ylab("Predicted Fraction")


p2 <- ggplot(frac_m_corrs, aes(plot_name, rho))+geom_col()+
  geom_errorbar(aes(ymin = conf.lower, ymax=conf.upper))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab("Cell Type")+
  ylab("Spearman Correlation")


ggarrange(p1, p2, labels = "AUTO")


# QQ Plots
qq1 <- ggplot(frac_dt[scenario=="Even Fractions"], aes(sample=fraction.true))+
  stat_qq()+
  stat_qq_line()+
  xlab("Theoretical Quantiles")+
  ylab("Cell Type Fraction Quantiles")+
  labs(title = "True Fractions")+
  theme(plot.title = element_text(hjust = 0.5))

qq2 <- ggplot(frac_dt[scenario=="Even Fractions"], aes(sample=fraction.predicted))+
  stat_qq()+
  stat_qq_line()+
  xlab("Theoretical Quantiles")+
  ylab("Cell Type Fraction Quantiles")+
  labs(title = "Predicted Fractions")+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(qq1, qq2)




# pairwise correlation

frac_pred_wide <- frac_dt[, c("sample", "sample_new", "fraction.predicted", "plot_name", "scenario")] %>%
  dcast(., ... ~plot_name, value.var = "fraction.predicted")

pmat <- corr.test(frac_pred_wide[scenario == "Even Fractions", !c("sample", "sample_new", "scenario")])$p

corrplot(corr = cor(frac_pred_wide[scenario == "Even Fractions", !c("sample", "sample_new", "scenario")]),
         p.mat = pmat,
         insig = "blank",
         sig.level = 0.05,
         type = "lower")



# special cases

## Missing CD8 T_RM

deconv <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_final_no_TCELL_CD8_RMeven_tpm_ct0_rep0/deconvolution.rds")
metrics <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_final_no_TCELL_CD8_RMeven_tpm_ct0_rep0/results_metric.rds")
# add missing cell type
tmp_met <- metrics$facs_groud_truth %>% as.data.frame()
tmp_met["TCELL_CD8_RM", ] <- rep(0, 10)
metrics$facs_groud_truth <- tmp_met %>% as.matrix()


deconv_m <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_final_no_TCELL_CD8_RMmirror_tpm_ct0_rep0/deconvolution.rds")
metrics_m <- readRDS("/nfs/data/tcell_deconvolution/data/bassez/final_results/dwls_finalDWLS_counts_final_no_TCELL_CD8_RMmirror_tpm_ct0_rep0/results_metric.rds")
# add missing cell type
tmp_met <- metrics_m$facs_groud_truth %>% as.data.frame()
tmp_met["TCELL_CD8_RM", ] <- rep(0, 10)
metrics_m$facs_groud_truth <- tmp_met %>% as.matrix()

frac_dt <- create_frac_tbl(deconv, metrics)
frac_m_dt <- create_frac_tbl(deconv_m, metrics_m)


frac_dt$scenario <- "Even Fractions"
frac_m_dt$scenario <- "Mirrored Fractions"

frac_dt <- rbind(frac_dt, frac_m_dt)



frac_dt <- merge(frac_dt, sample_lookup, by= "sample")


# heatmap
ggplot(frac_dt[scenario == "Even Fractions"], aes(sample_new, plot_name, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "D")+
  xlab("Sample")+
  ylab("Cell Type")

# heatmap mirror
ggplot(frac_dt[scenario == "Mirrored Fractions"], aes(sample_new, plot_name, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "D")+
  xlab("Sample")+
  ylab("Cell Type")

# heatmap norm
scale_breaks <- c(0, -sqrt(1), -sqrt(10), -sqrt(25))

ggplot(frac_dt[scenario == "Mirrored Fractions"], aes(sample_new, plot_name, fill=sqdscore(diff.d)))+geom_tile()+
  geom_text(aes(label= round(sqdscore(diff.d), 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "ND")+
  xlab("Sample")+
  ylab("Cell Type")


