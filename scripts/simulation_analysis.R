

aggregate_under_results <- function(results_folder, method, scenario) {
  folders <- list.dirs(path = results_folder,
                      recursive = FALSE, full.names = FALSE)
  folders <- folders[grepl(paste0("^", method), folders)]
  folders <- folders[grepl(scenario, folders)]
  
  lapply(folders, function(folder){
    namebuild <- strsplit(folder, "_", fixed = TRUE)[[1]]
    
    celltype <- namebuild[(which(namebuild=="no")+1):(length(namebuild)-3)]
    celltype <- paste(celltype, collapse = "_")
    celltype <- substring(celltype, 1, (nchar(celltype)-nchar(scenario)))
    
    if (scenario == "mirror") {
      scenario_n <- "mirror_db"
    } else {
      scenario_n <- scenario
    }
  

    deconv <- readRDS(file.path(results_folder, folder, "deconvolution.rds"))
    metrics <- readRDS(file.path(results_folder, folder, "results_metric.rds"))
    
    frac_dt <- create_frac_tbl(deconv, metrics, missing = celltype, scenario = scenario_n)
    frac_dt$missing <- celltype
    
    return(frac_dt)
  })
}

results_folder <- "/nfs/data/tcell_deconvolution/data/bassez/final_results"


frac_ls <- aggregate_under_results(results_folder = results_folder, method = "bayesprism", scenario = "even")

fraggregate <- rbindlist(frac_ls)

fraggregate <- merge(fraggregate, pretty_names[, c(1,2)], by.x="missing", by.y = "data_name")

relevant_types <- c("TCELL_CD4_EM", "TCELL_CD4_EX", "TCELL_CD4_EX_Proliferating",
                    "TCELL_CD4_REG", "TCELL_CD4_REG_Proliferating", "TCELL_CD8_EM",
                    "TCELL_CD8_EX", "TCELL_CD8_EX_Proliferating", "TCELL_CD8_RM")

non_relevant_types <- c("B_cell", "Cancer_cell", "Endothelial_cell", "Fibroblast",
                        "Mast_cell", "Myeloid_cell", "pDC", "TCELL_CD4_N",
                        "TCELL_CD8_EMRA", "TCELL_CD8_N", "TCELL_gdT", "TCELL_NK_CYTO", 
                        "TCELL_NK_REST", "TCELL_Vg9Vd2_gdT")

ggplot(fraggregate[missing %in% relevant_types], aes(plot_name.x, diff))+
  geom_boxplot(aes(color = (plot_name.y==plot_name.x)))+
  scale_color_manual("Missing", values = c("black", "red"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  geom_abline(slope=0)+
  facet_wrap(vars(plot_name.y))+
  xlab("Cell Type")+
  ylab("D")



mean_fragg <- fraggregate[, lapply(.SD, mean),
            by=c("cell_type", "missing", "plot_name.x", "plot_name.y"),
            .SDcols = c("fraction.true", "fraction.predicted", "diff")]  


# heatmap
ggplot(mean_fragg[missing %in% relevant_types], aes(plot_name.y, plot_name.x, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "D")+
  xlab("Missing Cell Type")+
  ylab("Cell Type")+
  geom_tile(data=mean_fragg[missing %in% relevant_types],
            aes(color = (plot_name.y==plot_name.x), size = (plot_name.y==plot_name.x)), alpha=0)+
  scale_colour_manual("Missing", values = c("white", "red", "black"))+ 
  scale_size_manual("Missing", values = c(0, 0.5))

# non relevant
ggplot(mean_fragg[!(missing %in% relevant_types)], aes(plot_name.y, plot_name.x, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "D")+
  xlab("Missing Cell Type")+
  ylab("Cell Type")+
  geom_tile(data=mean_fragg[!(missing %in% relevant_types)],
            aes(color = (plot_name.y==plot_name.x), size = (plot_name.y==plot_name.x)), alpha=0)+
  scale_colour_manual("Missing", values = c("white", "red", "black"))+ 
  scale_size_manual("Missing", values = c(0, 0.5))


aggregate_over_results <- function(results_folder, method) {
  folders <- list.dirs(path = results_folder,
                       recursive = FALSE, full.names = FALSE)
  folders <- folders[grepl(paste0("^", method), folders)]
  folders <- folders[grepl("_only_", folders)]
  
  lapply(folders, function(folder){
    namebuild <- strsplit(folder, "_", fixed = TRUE)[[1]]
    
    celltype <- namebuild[(which(namebuild=="only")+1):(length(namebuild)-3)]
    celltype <- paste(celltype, collapse = "_")

    deconv <- readRDS(file.path(results_folder, folder, "deconvolution.rds"))
    metrics <- readRDS(file.path(results_folder, folder, "results_metric.rds"))
    
    frac_dt <- create_frac_tbl(deconv, metrics)
    frac_dt$overrep <- celltype
    
    return(frac_dt)
    
  })
}

results_folder <- "/nfs/data/tcell_deconvolution/data/bassez/final_results"


frac_ls <- aggregate_over_results(results_folder = results_folder, method = "dwls")

fraggregate <- rbindlist(frac_ls)

fraggregate <- merge(fraggregate, pretty_names[, c(1,2)], by.x="overrep", by.y = "data_name")

relevant_types <- c("TCELL_CD4_EM", "TCELL_CD4_EX", "TCELL_CD4_EX_Proliferating",
                    "TCELL_CD4_REG", "TCELL_CD4_REG_Proliferating", "TCELL_CD8_EM",
                    "TCELL_CD8_EX", "TCELL_CD8_EX_Proliferating", "TCELL_CD8_RM")

ggplot(fraggregate[overrep %in% relevant_types], aes(plot_name.x, diff))+
  geom_boxplot(aes(color = (plot_name.y==plot_name.x)))+
  scale_color_manual("Overrepresented", values = c("black", "red"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  geom_abline(slope=0)+
  facet_wrap(vars(plot_name.y))+
  xlab("Cell Type")+
  ylab("D")

mean_fragg <- fraggregate[, lapply(.SD, mean),
                          by=c("cell_type", "overrep", "plot_name.x", "plot_name.y"),
                          .SDcols = c("fraction.true", "fraction.predicted", "diff")]  


# heatmap
ggplot(mean_fragg[overrep %in% relevant_types], aes(plot_name.y, plot_name.x, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "D")+
  xlab("Overrepresented Cell Type")+
  ylab("Cell Type")+
  geom_tile(data=mean_fragg[overrep %in% relevant_types],
            aes(color = (plot_name.y==plot_name.x), size = (plot_name.y==plot_name.x)), alpha=0)+
  scale_colour_manual("Overrepresented", values = c("white", "red", "black"))+ 
  scale_size_manual("Overrepresented", values = c(0, 0.5))

# non relevant
ggplot(mean_fragg[!(overrep %in% relevant_types)], aes(plot_name.y, plot_name.x, fill=diff))+geom_tile()+
  geom_text(aes(label= round(diff, 3)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_gradient2(name = "D")+
  xlab("Overrepresented Cell Type")+
  ylab("Cell Type")+
  geom_tile(data=mean_fragg[!(overrep %in% relevant_types)],
            aes(color = (plot_name.y==plot_name.x), size = (plot_name.y==plot_name.x)), alpha=0)+
  scale_colour_manual("Overrepresented", values = c("white", "red", "black"))+ 
  scale_size_manual("Overrepresented", values = c(0, 0.5))


