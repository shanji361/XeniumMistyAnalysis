data_path <- "C:/Users/juesh/UROP/data/misty_test"

results_folder <- "C:/Users/juesh/UROP/results/misty_test"


python_path <- NULL 


instructions <- createGiottoInstructions(save_dir = results_folder, 
                                         save_plot = TRUE, 
                                         show_plot = TRUE, 
                                         return_plot = FALSE, 
                                         python_path = python_path)

visium_lungcancer_test <- createGiottoVisiumObject(visium_dir = data_path,
                                                   expr_data = "raw",
                                                   png_name = "tissue_lowres_image.png",
                                                   gene_column_index = 2,
                                                   instructions = instructions)



# To flip the y-axis in your Visium object
visium_lungcancer_test <- flip(visium_lungcancer_test, direction = "vertical")

spatPlot(gobject = visium_lungcancer_test, 
         cell_color = "in_tissue", 
         show_image = TRUE, 
         point_alpha = 0.7)





pDataDT(visium_lungcancer_test)

# check available image names
showGiottoImageNames(visium_lungcancer_test) # "image" is the default name

# show aligned image
spatPlot(gobject = visium_lungcancer_test, 
         cell_color = "in_tissue", 
         show_image = TRUE, 
         point_alpha = 0.7)



visium_lungcancer_test<- filterGiotto(gobject = visium_lungcancer_test,
                                      expression_threshold = 1,
                                      feat_det_in_min_cells = 20,
                                      min_det_feats_per_cell = 500,
                                      expression_values = "raw",
                                      verbose = TRUE)

visium_lungcancer_test<- normalizeGiotto(gobject = visium_lungcancer_test, 
                                         scalefactor = 6000, 
                                         verbose = TRUE)

visium_lungcancer_test<- addStatistics(gobject = visium_lungcancer_test)



visium_lungcancer_test<- calculateHVF(gobject = visium_lungcancer_test)

visium_lungcancer_test<- runPCA(gobject = visium_lungcancer_test)

screePlot(visium_lungcancer_test, 
          ncp = 30)
visium_lungcancer_test<- runUMAP(visium_lungcancer_test, 
                                 dimensions_to_use = 1:10)

plotUMAP(gobject = visium_lungcancer_test)


visium_lungcancer_test<- runtSNE(visium_lungcancer_test, 
                                 dimensions_to_use = 1:10)

plotTSNE(gobject = visium_lungcancer_test)


# Create shared nearest network (SNN) and perform leiden clustering
visium_lungcancer_test<- createNearestNetwork(gobject = visium_lungcancer_test, 
                                              dimensions_to_use = 1:10, 
                                              k = 30)

visium_lungcancer_test<- doLeidenCluster(gobject = visium_lungcancer_test, 
                                         spat_unit = "cell", 
                                         feat_type = "rna", 
                                         resolution = 0.4, 
                                         n_iterations = 1000)

# # visualize UMAP cluster results
# plotUMAP(gobject = visium_lungcancer_test, 
#          cell_color = "leiden_clus", 
#          show_NN_network = TRUE, 
#          point_size = 2)

# # visualize tSNE cluster results
# plotTSNE(gobject = visium_lungcancer_test, 
#          cell_color = "leiden_clus", 
#          show_NN_network = TRUE, 
#          point_size = 2)
# 
# # visualize expression and spatial results
# spatDimPlot(gobject = visium_lungcancer_test, 
#             cell_color = "leiden_clus",
#             dim_point_size = 2, 
#             spat_point_size = 2)
# 
# spatDimPlot(gobject = visium_lungcancer_test, 
#             cell_color = "nr_feats", 
#             color_as_factor = FALSE,
#             dim_point_size = 2, 
#             dim_show_legend = TRUE, 
#             spat_show_legend = TRUE,
#             spat_point_size = 2)


# Cell type marker detection
# Gini markers
gini_markers_subclusters <- findMarkers_one_vs_all(gobject = visium_lungcancer_test,
                                                   method = "gini",
                                                   expression_values = "normalized",
                                                   cluster_column = "leiden_clus",
                                                   min_featss = 20,
                                                   min_expr_gini_score = 0.5,
                                                   min_det_gini_score = 0.5)

# get top 2 genes per cluster and visualize with violin plot
topgenes_gini <- gini_markers_subclusters[, head(.SD, 2), by = "cluster"]$feats

# violinPlot(visium_lungcancer_test, 
#            feats = unique(topgenes_gini), 
#            cluster_column = "leiden_clus",
#            strip_text = 8, 
#            strip_position = "right")
# 
# 
# # cluster heatmap
# plotMetaDataHeatmap(visium_lungcancer_test,
#                     selected_feats = topgenes_gini,
#                     metadata_cols = "leiden_clus",
#                     x_text_size = 10, 
#                     y_text_size = 10)
# 
# # umap plots
# dimFeatPlot2D(visium_lungcancer_test,
#               expression_values = "scaled",
#               feats = gini_markers_subclusters[, head(.SD, 1), by = "cluster"]$feats,
#               cow_n_col = 3,
#               point_size = 1)

# Cell type marker detection
# Scran markers
scran_markers_subclusters <- findMarkers_one_vs_all(gobject = visium_lungcancer_test,
                                                    method = "scran",
                                                    expression_values = "normalized",
                                                    cluster_column = "leiden_clus")

# get top 2 genes per cluster and visualize with violin plot
topgenes_scran <- scran_markers_subclusters[, head(.SD, 2), by = "cluster"]$feats

# violinPlot(visium_lungcancer_test, 
#            feats = unique(topgenes_scran),
#            cluster_column = "leiden_clus",
#            strip_text = 10, 
#            strip_position = "right")

# 
# # cluster heatmap
# plotMetaDataHeatmap(visium_lungcancer_test,
#                     selected_feats = topgenes_scran,
#                     metadata_cols = "leiden_clus",
#                     x_text_size = 10, 
#                     y_text_size = 10)
# 


# Create PAGE matrix
# PAGE matrix: each row a marker gene and each column a cell type

Tcells_markers <- c("CD2", "CD3D", "CD3E", "CD3G")
macrophage_markers <- c("MARCO", "CSF1R", "CD68", "GLDN", 
                        "APOE", "CCL3L1", "TREM2", "C1QB", 
                        "NUPR1", "FOLR2", "RNASE1", "C1QA")
dendritic_markers <- c("CD1E", "CD1C", "FCER1A", "PKIB", "CYP2S1", "NDRG2")
mast_markers <- c("CMA1", "TPSAB1", "TPSB2")
Bcell_markers <- c("IGLL5", "MZB1", "JCHAIN", "DERL3", "SDC1", 
                   "MS$A1", "BANK1", "PAX5", "CD79A")
Bcell_PB_markers <- c("PRDM1", "XSP1", "IRF4")
Bcell_mem_markers <- c("MS4A1", "IRF8")

neutrophils_markers <- c("FCGR3B", "ALPL", "CXCR1", "CXCR2", 
                         "ADGRG3", "CMTM2", "PROK2", "MME", "MMP25", "TNFRSF10C")
pdcs_markers <- c("SLC32A1", "SHD", "LRRC26", "PACSIN1", 
                  "LILRA4", "CLEC4C", "DNASE1L3", "SCT", "LAMP5")

# newly added cell type markers:
# Smooth muscle cells
smooth_muscle_markers <- c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "DES", "CALD1", "TPM2")

# Fibroblasts
fibroblast_markers <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRA", "THY1", "FBLN1", "VCAN")

# Alveolar Type I cells
atI_markers <- c("AGER", "PDPN", "CAV1", "RTKN2", "AQP5", "CLIC5", "EMP2", "HOPX")

# Alveolar Type II cells
atII_markers <- c("SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "SFTPD", "NAPSA", "LAMP3", "ABCA3", "PGC", "CLDN18")

# Ciliated cells
ciliated_markers <- c("FOXJ1", "TPPP3", "SPEF2", "PIFO", "DNAH5", "RSPH1", "CAPS", "DNAI1", "DNAI2", "SPAG6")

# Endothelial cells
endothelial_markers <- c("PECAM1", "CDH5", "VWF", "CLDN5", "CD34", "MCAM", "TEK", "KDR", "FLT1", "EMCN")

# Squamous cell carcinoma (some markers that might be elevated in SCC)
squamous_markers <- c("KRT5", "KRT14", "KRT6A", "TP63", "SOX2", "DSG3", "PKP1", "EGFR", "FGFR1", "S100A8")



# Update signature matrix with all cell types
signature_matrix <- makeSignMatrixPAGE(
  sign_names = c("T_Cells", "Macrophage", "Dendritic", "Mast", 
                 "B_cell", "Bcell_PB", "Bcells_memory",
                 "Neutrophils", "pDCs", 
                 "Smooth_Muscle", "Fibroblast", "Alveolar_TypeI", 
                 "Alveolar_TypeII", "Ciliated", "Endothelial", "Squamous_Carcinoma"),
  sign_list = list(Tcells_markers,
                   macrophage_markers,
                   dendritic_markers,
                   mast_markers,
                   Bcell_markers,
                   Bcell_PB_markers,
                   Bcell_mem_markers,
                   neutrophils_markers,
                   pdcs_markers,
                   smooth_muscle_markers,
                   fibroblast_markers,
                   atI_markers,
                   atII_markers,
                   ciliated_markers,
                   endothelial_markers,
                   squamous_markers))



visium_lungcancer_test<- runPAGEEnrich(gobject = visium_lungcancer_test, 
                                       sign_matrix = signature_matrix, 
                                       min_overlap_genes = 1)

cell_types <- colnames(signature_matrix)

# partition cells into groups of 4 for PAGE plots:
cell_types_subset1 <- cell_types[1:4] 
cell_types_subset2 <- cell_types[5:8]
cell_types_subset3 <- cell_types[9:12] 
cell_types_subset4 <- cell_types[13:16]  

library(tibble)
composition <- as_tibble(visium_lungcancer_test@spatial_enrichment$cell$rna$PAGE@enrichDT)


# plotMetaDataCellsHeatmap(gobject = visium_lungcancer_test,
#                          metadata_cols = "leiden_clus",
#                          value_cols = cell_types,
#                          spat_enr_names = "PAGE",
#                          x_text_size = 8,
#                          y_text_size = 8,
#                          show_plot = TRUE)
# 
# # 

library(viridis) 
viridis_colors <- viridis(7, option= "turbo")


# # Visualize first group of 4 PAGE cell plots 
spatCellPlot(gobject = visium_lungcancer_test,
             spat_enr_names = "PAGE",
             cell_annotation_values = cell_types_subset1,
             cow_n_col = 2,
             coord_fix_ratio = NULL,
             cell_color_gradient = viridis_colors,
             point_size = 1.5)

# # Visualize second group of 4 PAGE cell plots 
spatCellPlot(gobject = visium_lungcancer_test,
             spat_enr_names = "PAGE",
             cell_annotation_values = cell_types_subset2,
             cow_n_col = 2,
             coord_fix_ratio = NULL,
             cell_color_gradient = viridis_colors,
             point_size = 1.5)

# # Visualize third group of 4 PAGE cell plots 
spatCellPlot(gobject = visium_lungcancer_test,
             spat_enr_names = "PAGE",
             cell_annotation_values = cell_types_subset3,
             cow_n_col = 2,
             coord_fix_ratio = NULL,
             cell_color_gradient = viridis_colors,
             point_size = 1.5)

# # Visualize fourth group of 4 PAGE cell plots 
spatCellPlot(gobject = visium_lungcancer_test,
             spat_enr_names = "PAGE",
             cell_annotation_values = cell_types_subset4,
             cow_n_col = 2,
             coord_fix_ratio = NULL,
             cell_color_gradient = viridis_colors,
             point_size = 1.5)


page_results <- visium_lungcancer_test@spatial_enrichment$cell$rna$PAGE$enrichDT






# enrichment test with PAGE
markers_scran <- findMarkers_one_vs_all(gobject = visium_lungcancer_test, 
                                        method = "scran",
                                        expression_values = "normalized", 
                                        cluster_column = "leiden_clus", 
                                        min_feats = 3)

topgenes_scran <- markers_scran[, head(.SD, 10), by = "cluster"]

celltypes <- levels(factor(markers_scran$cluster))

sign_list <- list()

for (i in 1:length(celltypes)){
  sign_list[[i]] <- topgenes_scran[which(topgenes_scran$cluster == celltypes[i]),]$feats
}

PAGE_matrix <- makeSignMatrixPAGE(sign_names = celltypes,
                                  sign_list = sign_list)


visium_lungcancer_test<- createSpatialGrid(gobject = visium_lungcancer_test,
                                           sdimx_stepsize = 400,
                                           sdimy_stepsize = 400,
                                           minimum_padding = 0)

# spatPlot(visium_lungcancer_test, 
#          cell_color = "leiden_clus", 
#          point_size = 2.5, 
#          show_grid = TRUE,
#          grid_color = "red", 
#          spatial_grid_name = "spatial_grid")


## Delaunay network: stats + creation
plotStatDelaunayNetwork(gobject = visium_lungcancer_test, 
                        maximum_distance = 400)

visium_lungcancer_test<- createSpatialNetwork(gobject = visium_lungcancer_test, 
                                              minimum_k = 0)

showGiottoSpatNetworks(visium_lungcancer_test)
# 
# spatPlot(gobject = visium_lungcancer_test, 
#          show_network = TRUE,
#          network_color = "blue", 
#          spatial_network_name = "Delaunay_network")

# kmeans binarization
# km_spatialfeats <- binSpect(visium_lungcancer_test)
# 
# spatFeatPlot2D(visium_lungcancer_test,
#                expression_values = "scaled",
#                feats = km_spatialfeats$feats[1:6],
#                cow_n_col = 2,
#                point_size = 1.5)






# MISTy ANALYSIS BEGINS:

library(mistyR)
#for parallel processing 
library(future)
# data manipulation
#set row col names:
library(data.table)
#data cleaning:
library(janitor)
library(Matrix)
library(dplyr)
library(reticulate)
library(tidyverse)
#create progeny model with proper gene names
library(msigdbr)
# normalization
library(sctransform)
# resource
library(decoupleR)
# setup parallel execution
plan(multisession)
library(SpatialExperiment)
library(SingleCellExperiment)
library(SummarizedExperiment)
# library(GSVA)
library(msigdbr)
library(progeny)
library(OmnipathR)
library(progeny)
library(decoupleR)
library(mistyR)




# visium_lungcancer_test<- loadGiotto("results/misty_test/visium_lungcancer_object_test")


#---------------------------------------------
#LOAD GIOTTO OBJECT

saveGiotto(visium_lungcancer_test, "results/misty_test/visium_lungcancer_object_test", overwrite = TRUE)

visium_lungcancer_test <- loadGiotto(path_to_folder = "C:/Users/juesh/UROP/results/misty_test/visium_lungcancer_object_test")





#---------------------------------------------
#EXTRACT AND PREPARE EXPRESSION DATA:

raw_matrix <- visium_lungcancer_test@expression$cell$rna$raw@exprMat
# Get the row names (gene symbols)
symbols <- rownames(raw_matrix)

# Identify duplicates
d.index <- which(duplicated(symbols))

# Add suffixes to make duplicates unique
if(length(d.index) > 0) {
  symbols[d.index] <- paste0(symbols[d.index], "_", seq_along(d.index))
}
# Filter for genes expressed in at least 5% of spots
min_spots_percentage <- 0.05
min_spots <- ncol(raw_matrix) * min_spots_percentage
expression <- raw_matrix[rowSums(raw_matrix > 0) >= min_spots, ]




#---------------------------------------------
#EXTRACT SPATIAL COORDS:

# Get spatial coordinates
geometry <- getSpatialLocations(visium_lungcancer_test, output = "data.table")
geometry <- geometry[, c("sdimx", "sdimy")]
setnames(geometry, c("sdimx", "sdimy"), c("row", "col"))



#---------------------------------------------
#PATHWAY ANALYSIS WITH PROGENy:


# if (!requireNamespace("OmnipathR", quietly = TRUE)) {
#   install.packages("BiocManager")
#   BiocManager::install("OmnipathR")
# }
# 


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

#BiocManager::install("progeny")


# progeny_data <- static_table(query = "annotations", resource = "PROGENy", organism = 9606)



#loaded in data from above cmd:
progeny_data <- read.csv(gzfile("C:\\Users\\juesh\\UROP\\results\\misty_test\\progeny_static_table.csv.gz"))

#note: output from this function was sent via slack: progeny_static_table.gz

model <- progeny_data %>%
  group_by(pathway) %>%
  arrange(desc(abs(weight))) %>%
  slice_head(n = 1000) %>%
  ungroup() %>%
  select(source = pathway, target = genesymbol, weight) %>%
  distinct(source, target, .keep_all = TRUE)  # Key fix for duplicates

#est_path_act <- run_mlm(expression, model, .mor = NULL)

est_path_act <- read.csv(gzfile("C:\\Users\\juesh\\UROP\\results\\misty_test\\est_path_act.csv.gz"))

#note: output from this function was sent via slack: est_path_act.csv.gz

est_path_act_wide <- est_path_act %>% 
  pivot_wider(id_cols = condition, names_from = source, values_from = score) %>%
  column_to_rownames("condition") 

colnames(est_path_act_wide)  <- est_path_act_wide %>% 
  clean_names(parsing_option = 0) %>% 
  colnames(.)



visium_lungcancer_test@expression[["cell"]][["rna"]][["progeny"]] <- t(est_path_act_wide)
#note: this line above is a workaround (bc we have a giotto object, not seurat) 

pathway_activity <-t(as.matrix(Giotto::getExpression(visium_lungcancer_test, values = "progeny", spat_unit = "cell", feat_type = "rna", output = "matrix")))

#note: this next section is a huge workaround lol, I referred to this tutorial to figure out how to access via url alongside perplexity: https://workflows.omnipathdb.org/intercell-networks-r.html
interactions_url <- "https://omnipathdb.org/interactions?genesymbols=yes&datasets=omnipath,pathwayextra,ligrecextra&organisms=9606&fields=sources,references,curation_effort&license=academic"
interactions <- read.delim(interactions_url)

ligands_url <- "https://omnipathdb.org/intercell?scope=generic&categories=ligand&causality=trans&license=academic"
ligands <- read.delim(ligands_url)
ligands <- unique(ligands[, c("genesymbol", "category", "scope", "aspect", "source")])
colnames(ligands)[1] <- "source_genesymbol"

receptors_url <- "https://omnipathdb.org/intercell?scope=generic&categories=receptor&causality=rec&license=academic"
receptors <- read.delim(receptors_url)
receptors <- unique(receptors[, c("genesymbol", "category", "scope", "aspect", "source")])
colnames(receptors)[1] <- "target_genesymbol"

# Join: interactions where source is a ligand and target is a receptor
ligrec <- merge(interactions, ligands, by = "source_genesymbol")
ligrec <- merge(ligrec, receptors, by = "target_genesymbol")

#Note: if you are unable to run the code above, I have send the file via slack called this: ligrec.csv.gz 

ligands <- unique(ligrec$source_genesymbol)

#back track to beginning of code 
expression <- as.matrix(Giotto::getExpression(visium_lungcancer_test, values = "raw", spat_unit = "cell", feat_type = "rna", output = "matrix"))
gene_names <- rownames(expression[(rowSums(expression > 0) / ncol(expression)) >= 0.05,]) 

slide_markers <- ligands[ligands %in% gene_names] 
ligand_expr <- t(as.matrix(expression[slide_markers,])) %>% clean_names()

#I did not think we need to create an equivalent to this line: #rownames(seurat_vs@assays$SCT@data) <- seurat_vs@assays$SCT@data %>% clean_names(parsing_option = 0) %>% rownames(.)



spatFeatPlot2D(
  gobject = visium_lungcancer_test,
  spat_unit = 'cell',
  feat_type = "rna", 
  show_image = TRUE,
  feats = c("jak.stat", "hypoxia"),  
  expression_values = "progeny",     
  cell_color_gradient = viridis_colors,
  gradient_style = "s",
  point_size = 1.5,
  legend_text = 7,
  cow_n_col = 2                     
)


saveGiotto(visium_lungcancer_test, "results/misty_test/visium_lungcancer_object_test", overwrite = TRUE)


colnames(composition)  <- composition %>% clean_names(parsing_option = 0) %>% colnames(.)

comp_views <- create_initial_view(composition) 

path_act_views <- create_initial_view(est_path_act_wide) %>%
  add_juxtaview(geometry,  neighbor.thr = 130) %>% 
  add_paraview(geometry, l= 200, family = "gaussian")

com_path_act_views <- comp_views %>%
  add_views(create_view("juxtaview.path.130", path_act_views[["juxtaview.130"]]$data, "juxta.path.130"))%>% 
  add_views(create_view("paraview.path.200", path_act_views[["paraview.200"]]$data, "para.path.200")) 

# Only process actual views, skip metadata fields
view_names <- names(com_path_act_views)
view_names <- view_names[view_names != "misty.uniqueid"]

for (view_name in view_names) {
  # this will only process actual view objects
  com_path_act_views[[view_name]]$data <- apply(com_path_act_views[[view_name]]$data, 2, 
                                                function(x) as.numeric(as.character(x)))
}

# missing_values_report <- list()
# non_numeric_report <- list()
# 
# for (view_name in view_names) {
#   data_matrix <- com_path_act_views[[view_name]]$data
#   missing_values_report[[view_name]] <- any(is.na(data_matrix))
#   # Check if all columns are numeric
#   non_numeric_report[[view_name]] <- any(!sapply(data_matrix, is.numeric))
# }
# 
# print("Missing values report:")
# print(missing_values_report)
# print("Non-numeric data report:")
# print(non_numeric_report)
# 


#converts all data to numeric and inputes any NA values with 0 
for (view_name in view_names) {
  data_matrix <- com_path_act_views[[view_name]]$data
  # Convert to numeric
  data_matrix <- apply(data_matrix, 2, function(x) as.numeric(as.character(x)))
  # Impute NAs with 0
  data_matrix[is.na(data_matrix)] <- 0
  com_path_act_views[[view_name]]$data <- data_matrix
}



# Calculate variance for each column
var_per_col <- apply(com_path_act_views$intraview$data, 2, var)
# Find columns with zero variance
zero_var_cols <- names(var_per_col[var_per_col == 0])
print(zero_var_cols)

cols_to_keep <- names(var_per_col[var_per_col > 0])
com_path_act_views$intraview$data <- com_path_act_views$intraview$data[, cols_to_keep, drop=FALSE]

# Print confirmation
print(paste('Removed zero variance column:', 'cell.id'))

for (view_name in view_names) {
  # Check if data is a matrix and convert to data frame 
  if (is.matrix(com_path_act_views[[view_name]]$data)) {
    com_path_act_views[[view_name]]$data <- as.data.frame(com_path_act_views[[view_name]]$data)
    print(paste("Converted", view_name, "from matrix to data frame"))
  }
}

# Now try running MISTy again
result <- run_misty(com_path_act_views, "results/misty_test/comp_path_act")


misty_results_com_path_act <- collect_results("results/misty_test/comp_path_act/")

misty_results_com_path_act %>%
  plot_improvement_stats("intra.R2")%>%
  plot_improvement_stats("gain.R2") 


misty_results_com_path_act %>% 
  plot_view_contributions()



#EMPTY JUXTA PATH 
misty_results_com_path_act %>%
  plot_interaction_heatmap("juxta.path.130", clean = TRUE)


#RUN AFTER JUXTA PATH no longer empty:

# # Visualize gene expression accord to juxta path
# spatFeatPlot2D(
#   gobject = visium_lungcancer_test,
#   spat_unit = 'cell',
#   feat_type = "rna", 
#   show_image = TRUE,
#   feats = c("change_for_gene_or_cell", "change_for_gene_or_cell"),  
#   expression_values = "progeny",   
#   cell_color_gradient = viridis_colors,
#   gradient_style = "s",
#   point_size = 1.5,
#   legend_text = 7,
#   cow_n_col = 2                     
# )


# # Visualize cell expression accord. to prev gene
# spatFeatPlot2D(
#   gobject = visium_lungcancer_test,
#   spat_unit = 'cell',
#   feat_type = "rna", 
#   show_image = TRUE,
#   feats = c("change_for_gene_or_cell", "change_for_gene_or_cell"),  
#   expression_values = "progeny",   
#   cell_color_gradient = viridis_colors,
#   gradient_style = "s",
#   point_size = 1.5,
#   legend_text = 7,
#   cow_n_col = 2                     
# )


run_misty(com_path_act_views, "results/misty_test/comp_path_act_linear", model.function = linear_model, bypass.intra = TRUE)

misty_results_com_path_act_linear <- collect_results("results/misty_test/comp_path_act_linear")



misty_results_com_path_act_linear %>%
  plot_improvement_stats("gain.R2") %>%
  plot_view_contributions()


#EMPTY JUXTA PATH:
misty_results_com_path_act_linear %>%
  plot_interaction_heatmap("juxta.path.130", clean = TRUE) 


#RUN AFTER JUXTA PATH no longer empty:

# # Visualize gene expression accord to juxta path
# spatFeatPlot2D(
#   gobject = visium_lungcancer_test,
#   spat_unit = 'cell',
#   feat_type = "rna", 
#   show_image = TRUE,
#   feats = c("change_for_gene_or_cell", "change_for_gene_or_cell"),  
#   expression_values = "progeny",   
#   cell_color_gradient = viridis_colors,
#   gradient_style = "s",
#   point_size = 1.5,
#   legend_text = 7,
#   cow_n_col = 2                     
# )


# # Visualize cell expression accord. to prev gene
# spatFeatPlot2D(
#   gobject = visium_lungcancer_test,
#   spat_unit = 'cell',
#   feat_type = "rna", 
#   show_image = TRUE,
#   feats = c("change_for_gene_or_cell", "change_for_gene_or_cell"),  
#   expression_values = "progeny",   
#   cell_color_gradient = viridis_colors,
#   gradient_style = "s",
#   point_size = 1.5,
#   legend_text = 7,
#   cow_n_col = 2                     
# )



pathway_act_view <- create_initial_view(as_tibble(pathway_activity) ) %>%
  add_paraview(geometry, l = 10, family = "constant")


ligand_view <- create_initial_view(as_tibble(ligand_expr)  %>% clean_names()) %>%
  add_paraview(geometry, l = 10, family = "constant")


combined_views <- pathway_act_view %>% add_views(create_view("paraview.ligand.10", ligand_view[["paraview.10"]]$data, "para.ligand.10"))


run_misty(combined_views, "results/misty_test/functional_ligand")

misty_results <- collect_results("results/misty_test/functional_ligand/")


misty_results %>%
  plot_improvement_stats("multi.R2") %>%
  plot_improvement_stats("gain.R2")

misty_results %>% plot_view_contributions()

misty_results %>%
  plot_interaction_heatmap("intra", clean = TRUE, cutoff = 1.5)




# # Visualize gene expression accord to interactions plot
# spatFeatPlot2D(
#   gobject = visium_lungcancer_test,
#   spat_unit = 'cell',
#   feat_type = "rna", 
#   show_image = TRUE,
#   feats = c("change_for_gene_or_cell", "change_for_gene_or_cell"),  
#   expression_values = "progeny",   
#   cell_color_gradient = viridis_colors,
#   gradient_style = "s",
#   point_size = 1.5,
#   legend_text = 7,
#   cow_n_col = 2                     
# )


# # Visualize cell expression accord. to prev gene
# spatFeatPlot2D(
#   gobject = visium_lungcancer_test,
#   spat_unit = 'cell',
#   feat_type = "rna", 
#   show_image = TRUE,
#   feats = c("change_for_gene_or_cell", "change_for_gene_or_cell"),  
#   expression_values = "progeny",   
#   cell_color_gradient = viridis_colors,
#   gradient_style = "s",
#   point_size = 1.5,
#   legend_text = 7,
#   cow_n_col = 2                     
# )


#We can observe that egfr is a significant predictor 
#for the activity of the mapk pathway when in the same 
#spot. Let’s take a look at the spatial distribution of 
#these pathway activities in the tissue slide:
spatFeatPlot2D(
  gobject = visium_lungcancer_test,
  spat_unit = 'cell',
  feat_type = "rna", 
  show_image = TRUE,
  feats = c("egfr", "mapk"),  
  expression_values = "progeny",   
  cell_color_gradient = viridis_colors,
  gradient_style = "s",
  point_size = 1.5,
  legend_text = 7,
  cow_n_col = 2                     
)



misty_results %>%
  plot_interaction_heatmap(view = "para.10", 
                           clean = TRUE, 
                           trim = 0.5,
                           trim.measure = "gain.R2",
                           cutoff = 1.25)


spatFeatPlot2D(
  gobject = visium_lungcancer_test,
  spat_unit = 'cell',
  feat_type = "rna", 
  show_image = TRUE,
  feats = c("androgen", "vegf"),  
  expression_values = "progeny",   
  cell_color_gradient = viridis_colors,
  gradient_style = "s",
  point_size = 1.5,
  legend_text = 7,
  cow_n_col = 2                     
)

misty_results %>%
  plot_interaction_heatmap(view = "para.ligand.10", clean = TRUE, trim = 0.5,
                           trim.measure = "gain.R2", cutoff=3)



#MYL9 is a predictor of pi3k and tgfb 

spatFeatPlot2D(
  gobject = visium_lungcancer_test,
  spat_unit = 'cell',
  feat_type = "rna",
  show_image = TRUE,
  feats = c("pi3k", "tgfb"),
  expression_values = "progeny",
  cell_color_gradient = viridis_colors,
  gradient_style = "s",
  point_size = 1.5,
  legend_text = 7,
  cow_n_col = 2
)



spatFeatPlot2D(
    gobject = visium_lungcancer_test,
    show_image = TRUE,
    point_alpha = 2,
    cell_color_gradient = viridis_colors,
    gradient_style = "s",
    point_size = 3,
    feats = "MYL9"
)


saveGiotto(visium_lungcancer_test, "results/misty_test/visium_lungcancer_object_test", overwrite = TRUE)




