setwd("C:/Users/juesh/UROP")
library(reticulate)
use_python("C:/Users/juesh/anaconda3/envs/giotto_env/python.exe", required = TRUE)
library(Giotto)


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
library(tidyr)
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



# set up paths
data_path <- "/Users/juesh/UROP/data/xenium_lung"
save_dir <- '/Users/juesh/UROP/results/xenium_lung_new_param'
dir.create(save_dir, recursive = TRUE)

# # download the mini dataset and untar
# options("timeout" = Inf)
# download.file(
#   url = "https://zenodo.org/records/13207308/files/workshop_xenium.zip?download=1",
#   destfile = file.path(save_dir, "workshop_xenium.zip"),
#   mode = "wb",           # Force binary mode
#   method = "libcurl",    # Use libcurl method
#   timeout = 300          # Set longer timeout
# )



unzip(file.path(save_dir, "workshop_xenium.zip"), 
      exdir = data_path)


#-------------------------------------------------------------------------- Giotto Workshop 2024 Tutorial Starts Here ----------------------------------------------------------------
g <- createGiottoXeniumObject(xenium_dir = data_path)

# set instructions for save directory and to save the plots to disk
instructions(g, "save_dir") <- save_dir
instructions(g, "save_plot") <- TRUE



g <- addSpatialCentroidLocations(g, poly_info = "cell")
g <- addSpatialCentroidLocations(g, poly_info = "nucleus")

# spatInSituPlotPoints(g,
#                      polygon_feat_type = "cell",
#                      feats = list(rna = head(featIDs(g))), # must be named list
#                      use_overlap = FALSE, 
#                      polygon_color = "cyan", 
#                      polygon_line_size = 0.1
# )

#####see if this removal still works######

x <- importXenium(data_path)

force(x)

x$qv <- 20 # default
tx <- x$load_transcripts()

plot(tx[[1]])


rm(tx) # remove to save space

x$filetype$expression <- "mtx" # change to mtx instead of .h5 which is not in the mini dataset

ex <- x$load_expression()
featType(ex)

force(g)


featType(ex[[2]]) <- c("NegControlProbe")
featType(ex[[3]]) <- c("NegControlCodeword")
featType(ex[[4]]) <- c("UnassignedCodeword")



g2 <- g
# append the expression info
g2 <- setGiotto(g2, ex)

# load cell metadata
cx <- x$load_cellmeta()
g2 <- setGiotto(g2, cx)

force(g2)



rm(g2) # save space


# Check if morphology images exist in your Xenium directory
xenium_dir <- "data/xenium_lung"  # Directory, not file
image_files <- list.files(xenium_dir, pattern = "morphology", full.names = TRUE)



img_paths <- c(
  sprintf("data/xenium_lung/morphology_focus/morphology_focus_%04d.tif", 0:3),
  "data/xenium_lung/he_mini.tif"
)

img_list <- createGiottoLargeImageList(
  img_paths, 
  # naming is based on the channel metadata above
  names = c("DAPI", "18S", "ATP1A1/CD45/E-Cadherin", "alphaSMA/Vimentin", "HE"),
  use_rast_ext = TRUE,
  verbose = FALSE
)

# make some images brighter
img_list[[1]]@max_window <- 5000
img_list[[2]]@max_window <- 5000
img_list[[3]]@max_window <- 5000

# append images to gobject
g <- setGiotto(g, img_list)



# # example plots
# spatInSituPlotPoints(g,
#                      show_image = TRUE,
#                      image_name = "HE",
#                      polygon_feat_type = "cell",
#                      polygon_color = "cyan",
#                      polygon_line_size = 0.1,
#                      polygon_alpha = 0
# )
# 
# 
# 
# 
# spatInSituPlotPoints(g,
#                      show_image = TRUE,
#                      image_name = "DAPI",
#                      polygon_feat_type = "nucleus",
#                      polygon_color = "cyan",
#                      polygon_line_size = 0.1,
#                      polygon_alpha = 0
# )
# 
# 
# spatInSituPlotPoints(g,
#                      show_image = TRUE,
#                      image_name = "18S",
#                      polygon_feat_type = "cell",
#                      polygon_color = "cyan",
#                      polygon_line_size = 0.1,
#                      polygon_alpha = 0
# )
# 
# 
# 
# spatInSituPlotPoints(g,
#                      show_image = TRUE,
#                      image_name = "ATP1A1/CD45/E-Cadherin",
#                      polygon_feat_type = "nucleus",
#                      polygon_color = "cyan",
#                      polygon_line_size = 0.1,
#                      polygon_alpha = 0
# )


g <- calculateOverlap(g,
                      spatial_info = "cell",
                      feat_info = "rna"
)

g <- overlapToMatrix(g)
g <- addStatistics(g, expression_values = "raw")

cell_stats <- pDataDT(g)
ggplot2::ggplot(cell_stats, ggplot2::aes(total_expr)) +
  ggplot2::geom_histogram(binwidth = 5)


# very permissive filtering, mainly for removing 0 values
g <- filterGiotto(g,
                  expression_threshold = 1,
                  feat_det_in_min_cells = 1,
                  min_det_feats_per_cell = 5
)


g <- normalizeGiotto(g)
# overwrite original results with those for normalized values
g <- addStatistics(g)

# spatInSituPlotPoints(g,
#                      polygon_fill = "nr_feats",
#                      polygon_fill_gradient_style = "sequential",
#                      polygon_fill_as_factor = FALSE
# )
# 
# 
# spatInSituPlotPoints(g,
#                      polygon_fill = "total_expr",
#                      polygon_fill_gradient_style = "sequential",
#                      polygon_fill_as_factor = FALSE
# )

g <- runPCA(g, feats_to_use = NULL)


screePlot(g, ncp = 30)


g <- runUMAP(g, 
             dimensions_to_use = seq(15), 
             n_neighbors = 40 # default
)

Giotto::plotPCA(g)
plotUMAP(g)


g <- createNearestNetwork(g,
                          dimensions_to_use = seq(15), 
                          k = 40
)


g <- doLeidenCluster(g)

# 
# plotPCA_3D(g, 
#            cell_color = "leiden_clus", 
#            point_size = 1
# )
# 
# 
# plotUMAP(g, 
#          cell_color = "leiden_clus", 
#          point_size = 0.1, 
#          point_shape = "no_border"
# )
# 
# spatInSituPlotPoints(g,
#                      polygon_fill = "leiden_clus",
#                      polygon_fill_as_factor = TRUE,
#                      polygon_alpha = 1,
#                      show_image = TRUE,
#                      image_name = "HE"
# )

#----------------------------------- Giotto Workshop 2024 Tutorial Ends Here ----------------------------------------------------------------


# scran marker gene
res_scran <- findMarkers_one_vs_all(g, 
                                    cluster_column = "leiden_clus", 
                                    method = "scran",
                                    expression_values = "normalized"
)

# top 2 genes per cluster
topgenes_scran <- res_scran[, head(.SD, 2), by = 'cluster']


if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "20_ViolinPlot.png"), 
    width = 12, height = 30, units = "in", res = 300)



violinPlot(g,
           feats = unique(topgenes_scran$feats),
           cluster_column = "leiden_clus",
           save_param = list(base_height = 30,base_width = 12, save_dir = save_dir)
)

dev.off()



# writeChatGPTquery(
#   DEG_output = res_scran,
#   top_n_genes = 25,
#   tissue_type = 'lung_cancer',
#   folder_name = "C:/Users/juesh/UROP/results/xenium_lung/",
#   file_name = 'GPTQuery_Scran.txt'
# )

# 
# GiottoVisuals::dotPlot(
#   g,
#   spat_unit = "cell",
#   feats = c(   "FOXJ1", "DNAAF1", "SCGB2A1", #Bronchial Epithelial (Ciliated/ Club)
#                "VWF", "PECAM1",              #Endothelial
#                "PDGFRA", "COL5A2",           #Fibroblast  
#                "PDGFRB",                     #Pericyte
#                "MYH11", "ACTA2",             #Smooth Muscle
#                "SOX2",                       #Basal cells
#                "SFTA2", "ACE2",              #Alveolar Epithelial Type 2 
#                "AGER", "PDPN",               #Alveolar Epithelial Type 1
#                "ERBB2", "EGFR", "EPCAM", "KRT7", #LUAD Cancer
#                "MET",  "MYC",                #Oncogenes
#                "NKG7", "CD3E",               #NKcell/ Tcells
#                "MS4A1", "CD19", "CD79A", "MZB1", #Bcell/ Plasma 
#                "CD68", "MRC1", "CD14",       #Macrophage (Tissue-Resident: Alveolar, Interstitial) / Monocytes
#                "CD83", "CD86",               #Dendritic 
#                "KIT", "MS4A2"                #Granulocytes (Mast/ etc...)
#   ), 
#   cluster_column = "leiden_clus",
#   dot_size = function(x) mean(x != 0) * 100,
#   dot_size_threshold = 0,
#   dot_scale = 6,
#   dot_color = mean,
#   dot_color_gradient = c("royalblue3", 'orangered', "yellow"),
#   gradient_style = "s",
#   expression_values = "normalized",
#   show_legend = TRUE,
#   legend_text = 10,
#   legend_symbol_size = 2,
#   background_color = "white",
#   axis_text = 10, 
#   default_save_name = "dotPlot",
#   save_param = list(base_height = 8,base_width = 5)
# )



if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "21_Dotplot.png"), 
    width = 12, height = 8, units = "in", res = 300)


GiottoVisuals::dotPlot(
  g,
  spat_unit = "cell",
  feats = c(   "FOXJ1", "DNAAF1", "SCGB2A1", #Bronchial Epithelial (Ciliated/ Club)
               "VWF", "PECAM1",              #Endothelial
               "PDGFRA", "COL5A2",           #Fibroblast  
               "PDGFRB",                     #Pericyte
               "MYH11", "ACTA2",             #Smooth Muscle
               "SOX2",                       #Basal cells
               "SFTA2", "ACE2",              #Alveolar Epithelial Type 2 
               "AGER", "PDPN",               #Alveolar Epithelial Type 1
               "ERBB2", "EGFR", "EPCAM", "KRT7", #LUAD Cancer
               "MET",  "MYC",                #Oncogenes
               "NKG7", "CD3E",               #NKcell/ Tcells
               "MS4A1", "CD19", "CD79A", "MZB1", #Bcell/ Plasma 
               "CD68", "MRC1", "CD14",       #Macrophage (Tissue-Resident: Alveolar, Interstitial) / Monocytes
               "CD83", "CD86",               #Dendritic 
               "KIT", "MS4A2"                #Granulocytes (Mast/ etc...)
  ), 
  cluster_column = "leiden_clus",
  dot_size = function(x) mean(x != 0) * 100,
  dot_size_threshold = 0,
  dot_scale = 6,
  dot_color = mean,
  dot_color_gradient = c("royalblue3", 'orangered', "yellow"),
  gradient_style = "s",
  expression_values = "normalized",
  show_legend = TRUE,
  legend_text = 10,
  legend_symbol_size = 2,
  background_color = "white",
  axis_text = 10, 
  default_save_name = "dotPlot",
  save_param = list(base_height = 8,base_width = 12)
)
dev.off()



# Representative single marker genes for broad cell types
single_marker_genes <- c(
  "CD3E",     # NK / T cells
  "SFTA2",    # Alveolar epithelial (AT2-like / LUAD cancer)
  "PDGFRA",   # Fibroblasts / Stromal
  "CD68",     # Myeloid (Macrophages/Monocytes/DCs)
  "MS4A1",    # B cells
  "FOXJ1",    # Ciliated epithelial
  "PECAM1",   # Endothelial cells
  "SOX2",     # Basal cells
  "ACTA2",    # Smooth muscle
  "MZB1",     # Plasma cells
  "KIT"       # Mast cells
)

# 
# 
# 

if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "10A_DimFeatPlot2D.png"), 
    width = 18, height = 20, units = "in", res = 300)


dimFeatPlot2D(g,
              expression_values = "normalized",
              feats = single_marker_genes,
              dim_reduction_to_use = "umap",
              cow_n_col = 3,
              point_size = 1,
              cell_color_gradient = c("blue", "green"),
              save_param = list(base_height = 20, base_width = 18))

dev.off()

cell_types <- c(
  "NK / T cells",
  "Alveolar Epithelial cells (LUAD CANCER)",  #Type 1, Type 2
  "Stromal (Fibroblasts/ Pericytes)",
  "Myeloid (Macrophages / Monocytes) and Dendritic cells",
  "B cells",
  "Bronchial Epithelial (Ciliated) cells", #tumor?  tumor possibility exists but needs further context
  "Endothelial cells",
  "Basal cells",      #could be transformed in some tumors (?)
  "Alveolar Epithelial cells (LUAD CANCER)",
  "Smooth muscle cells",
  "Plasma cells",
  "Granulocytes (Mast)"
)


# Assign numeric names to the vector
names(cell_types) <- 1:length(cell_types)

# Annotate object
g <- annotateGiotto(gobject = g, 
                    spat_unit = "cell", 
                    annotation_vector = cell_types,
                    cluster_column = 'leiden_clus', 
                    name = 'subannot_clus')


xenium_lungcancer_test <- g 


#raw_matrix <- xenium_lungcancer_test@expression$cell$rna$raw[]


#symbols <- rownames(raw_matrix)

#d.index <- which(duplicated(symbols))
#print(d.index) 
#no duplicates found 

norm_matrix <- Giotto::getExpression(xenium_lungcancer_test, values = "normalized", output = "matrix")



geometry <- getSpatialLocations(
  gobject = xenium_lungcancer_test, 
  spat_unit = "cell",  #default
  output = "data.table"
)
geometry <- geometry[, c("sdimx", "sdimy")]
setnames(geometry, c("sdimx", "sdimy"), c("row", "col"))




#expression <- as.matrix(xenium_lungcancer_test@expression$cell$rna$raw[])



#!!!!!!!!!!!!!!REFERENCE MISTY TUTORIAL (CODE DIRECTLY USED)!!!!!!!!!!!!!!!
model <- get_progeny(organism = "human", top = 500) 

# Use multivariate linear model to estimate activity
est_path_act <- run_mlm(norm_matrix, model,.mor = NULL) 


est_path_act_wide <- est_path_act %>% 
  pivot_wider(id_cols = condition, names_from = source, values_from = score) %>%
  column_to_rownames("condition") 

colnames(est_path_act_wide) <- est_path_act_wide %>% 
  clean_names(parsing_option = 0) %>% 
  colnames(.)


# Add progeny results to object
path_act_exprobj = createExprObj(t(est_path_act_wide), name = "progeny")
xenium_lungcancer_test <- setExpression(xenium_lungcancer_test, path_act_exprobj, name = "progeny") #cell and rna are default



metadata <- getCellMetadata(g, spat_unit = "cell")

# Extract cell types from your Xenium metadata
cell_types <- metadata$subannot_clus

# Create one-hot encoded matrix 
cell_type_factor <- factor(cell_types)
cell_type_onehot <- model.matrix(~ cell_type_factor - 1)
colnames(cell_type_onehot) <- gsub("cell_type_factor", "", colnames(cell_type_onehot))

# Set proper row names using Xenium cell IDs
actual_cell_ids <- colnames(xenium_lungcancer_test@expression$cell$rna$raw)
rownames(cell_type_onehot) <- actual_cell_ids

# Clean names 
colnames(cell_type_onehot) <- cell_type_onehot %>% 
  as_tibble() %>%
  clean_names(parsing_option = 0) %>% 
  colnames(.)

# Convert to tibble for Misty
composition_xenium <- as_tibble(cell_type_onehot)





#=================trying param function above===========

# ==============================================================================
# MISTY SPATIAL TRANSCRIPTOMICS ANALYSIS PIPELINE
# ==============================================================================

# ==============================================================================
# STEP 1: CREATE SPATIAL VIEWS
# ==============================================================================

# Create pathway activity spatial views
path_act_views <- create_initial_view(est_path_act_wide) %>%
  add_juxtaview(geometry, neighbor.thr = 10) %>% 
  add_paraview(geometry, l = 15, family = "gaussian")


# Create cell composition spatial views
comp_views <- create_initial_view(composition_xenium) %>%
  add_juxtaview(geometry, neighbor.thr = 10) %>%
  add_paraview(geometry, l = 15, family = "gaussian")

final_misty_views <- path_act_views %>%
  add_views(create_view("juxtaview.composition.10", 
                        comp_views[["juxtaview.10"]]$data)) %>% 
  add_views(create_view("paraview.composition.15", 
                        comp_views[["paraview.15"]]$data)) 




# ==============================================================================

#RUN MISTY ANALYSES
# ==============================================================================


# Standard MISTy analysis (with intrinsic view)
run_misty(
  views = final_misty_views,  # Updated to use complete views
  cv.folds = 10,  
  model.function = random_forest_model,
  results.folder = file.path(save_dir, "misty_results_complete") # Updated folder name
)

# Spatial-only analysis (bypass intrinsic view)
# Tests purely spatial predictive power without cell's own composition
run_misty(view = final_misty_views, 
          cv.folds = 10, 
          model.function = linear_model, 
          results.folder = file.path(save_dir, "misty_results_lm_complete"), 
          bypass.intra = TRUE)


#collect results 

misty_results_complete <- collect_results(file.path(save_dir, "misty_results_complete"))
misty_results_complete_linear <- collect_results(file.path(save_dir, "misty_results_lm_complete"))




graphics.off()

if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "1_CompleteIntra.png"), 
    width = 12, height = 8, units = "in", res = 300)

misty_results_complete %>%
  plot_improvement_stats("intra.R2") 

dev.off()

graphics.off()

if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "2_CompleteGain.png"), 
    width = 12, height = 8, units = "in", res = 300)

misty_results_complete %>%
  plot_improvement_stats("gain.R2")

dev.off()


graphics.off()

if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "3_SpatialGain.png"), 
    width = 12, height = 8, units = "in", res = 300)
misty_results_complete_linear %>%
  plot_improvement_stats("gain.R2")

dev.off()



graphics.off()

if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "4_CompleteContributions.png"), 
    width = 12, height = 8, units = "in", res = 300)
misty_results_complete %>% 
  plot_view_contributions() 


dev.off()


graphics.off()

if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "5_SpatialContributions.png"), 
    width = 12, height = 8, units = "in", res = 300)
misty_results_complete_linear %>%
  plot_view_contributions()


dev.off()


# ==============================================================================

#Interaction Heatmaps 
# ==============================================================================





# clear all graphics devices
graphics.off()


# Only if the above displays correctly, then save:
if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "6_CompletePathwayJuxta.png"), 
    width = 12, height = 8, units = "in", res = 300)
misty_results_complete %>%
  plot_interaction_heatmap("juxta.10", clean = TRUE)
dev.off()



# clear all graphics devices
graphics.off()



# Only if the above displays correctly, then save:
if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "7_CompleteCompositionJuxta.png"), 
    width = 12, height = 8, units = "in", res = 300)
misty_results_complete %>%
  plot_interaction_heatmap("juxtaview.composition.10", clean = TRUE)
dev.off()


# clear all graphics devices
graphics.off()


# Only if the above displays correctly, then save:
if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "8_CompletePathwayPara.png"), 
    width = 12, height = 8, units = "in", res = 300)
misty_results_complete %>%
  plot_interaction_heatmap("para.15", clean = TRUE)
dev.off()


# clear all graphics devices
graphics.off()



# Only if the above displays correctly, then save:
if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "9_CompleteCompositionPara.png"),
    width = 12, height = 8, units = "in", res = 300)
misty_results_complete %>%
  plot_interaction_heatmap("paraview.composition.15", clean = TRUE)
dev.off()

# ==============================================================================

# SPATIAL-ONLY ANALYSIS HEATMAPS
# ==============================================================================



if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "10_SpatialPathwayJuxta.png"),
    width = 12, height = 8, units = "in", res = 300)

misty_results_complete_linear %>%
  plot_interaction_heatmap("juxta.10", clean = TRUE) 

dev.off()



if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "11_SpatialCompositionJuxta.png"),
    width = 12, height = 8, units = "in", res = 300)

misty_results_complete_linear %>%
  plot_interaction_heatmap("juxtaview.composition.10", clean = TRUE)

dev.off()



if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "12_SpatialPathwayPara.png"),
    width = 12, height = 8, units = "in", res = 300)

misty_results_complete_linear %>%
  plot_interaction_heatmap("para.15", clean = TRUE) 

dev.off()





if (dev.cur() != 1) dev.off()  # Close any open device
png(file.path(save_dir, "13_SpatialCompositionPara.png"),
    width = 12, height = 8, units = "in", res = 300)

misty_results_complete_linear %>%
  plot_interaction_heatmap("paraview.composition.15", clean = TRUE)

dev.off()



#spatial visualization examples 

# Example 1 : B cells and TRAIL pathway
# TRAIL is involved in immune-mediated apoptosis - expect high activity near immune cells

# Visualize TRAIL pathway activity

spatFeatPlot2D(xenium_lungcancer_test,
               spat_unit = "cell", 
               expression_values = "progeny", 
               show_image = TRUE, 
               feats = "trail", 
               gradient_style = "sequential", 
               cell_color_gradient = c("blue", "orangered", "yellow"), 
               background_color = "black", 
               point_size = 1, 
               save_plot = TRUE,   # Enable automatic saving
               save_param = list(
                 base_height = 8, 
                 base_width = 12, 
                 dpi = 600,
                 units = "in",
                 save_format = "png",
                 save_name = "14_TRAILPathway",
                 save_dir = save_dir  # Save to current working directory
               ))

# Visualize B cell locations
spatPlot2D(xenium_lungcancer_test,
           spat_unit = "cell", 
           cell_color = "subannot_clus", 
           show_image = TRUE, 
           select_cell_groups = "B cells", 
           point_size = 1, 
           other_point_size = 0.7, 
           other_cell_color = "#434343", 
           background_color = "black", 
           save_plot = TRUE,   # Enable automatic saving
           save_param = list(
             base_height = 8, 
             base_width = 12, 
             dpi = 600,
             units = "in",
             save_format = "png",
             save_name = "15_BCellLocations",
             save_dir = save_dir  # Save to current working directory
           ))


# Example 2: LUAD cancer cells and EGFR pathway
# EGFR is frequently dysregulated in lung adenocarcinoma - expect high activity in cancer regions

# Visualize EGFR pathway activity

spatFeatPlot2D(xenium_lungcancer_test,
               spat_unit = "cell", 
               expression_values = "progeny", 
               show_image = TRUE, 
               feats = "egfr",  
               gradient_style = "sequential", 
               cell_color_gradient = c("blue", "orangered", "green"), 
               background_color = "black", 
               point_size = 1, 
               save_plot = TRUE,   # Enable automatic saving
               save_param = list(
                 base_height = 8, 
                 base_width = 12, 
                 dpi = 600,
                 units = "in",
                 save_format = "png",
                 save_name = "16_EGFRPathway",
                 save_dir = save_dir  # Save to current working directory
               ))

# Visualize LUAD cancer cell locations

spatPlot2D(xenium_lungcancer_test,
           spat_unit = "cell", 
           cell_color = "subannot_clus", 
           show_image = TRUE, 
           select_cell_groups = "Alveolar Epithelial cells (LUAD CANCER)", 
           point_size = 1.1, 
           other_point_size = 0.7, 
           other_cell_color = "#434343", 
           background_color = "black", 
           save_plot = TRUE,   # Enable automatic saving
           save_param = list(
             base_height = 8, 
             base_width = 12, 
             dpi = 600,
             units = "in",
             save_format = "png",
             save_name = "17_LUADCancerLocations",
             save_dir = save_dir  # Save to current working directory
           ))




# Example 3 : NK/T cells and NFKB PATHWAY
# NFκB is key in immune activation - expect high activity in immune cell regions
# In lung cancer, NK/T cells often cluster near tumor boundaries

# Visualize NFκB pathway activity

spatFeatPlot2D(xenium_lungcancer_test,
               spat_unit = "cell", 
               show_image = TRUE, 
               expression_values = "progeny", 
               feats = "nfkb",  
               gradient_style = "sequential", 
               cell_color_gradient = c("blue","purple", "red", "yellow"), 
               background_color = "black", 
               point_size = 1, 
               save_plot = TRUE,   # Enable automatic saving
               save_param = list(
                 base_height = 8, 
                 base_width = 12, 
                 dpi = 600,
                 units = "in",
                 save_format = "png",
                 save_name = "18_NFKBPathway",
                 save_dir = save_dir  # Save to current working directory
               ))

# Visualize NK/T cell locations

spatPlot2D(xenium_lungcancer_test,
           spat_unit = "cell", 
           show_image = TRUE, 
           cell_color = "subannot_clus", 
           select_cell_groups = "NK / T cells", 
           point_size = 1, 
           other_point_size = 0.6, 
           other_cell_color = "#434343", 
           background_color = "black", 
           save_plot = TRUE,   # Enable automatic saving
           save_param = list(
             base_height = 8, 
             base_width = 12, 
             dpi = 600,
             units = "in",
             save_format = "png",
             save_name = "19_NKTCellLocations",
             save_dir = save_dir  # Save to current working directory
           ))





#analysis notes 

# Interpretation:
# - High intra.R2: Cell type identity strongly predicts its pathway activity
# - High gain.R2: Spatial neighborhood adds significant predictive power beyond cell identity
# - In Xenium data, we expect strong spatial patterns due to preserved tissue architecture


# Z-score normalization is important for pathway activities because:
# Different pathways have different baseline expression ranges
# Normalization ensures fair comparison across pathways in MISTy models
# Standard practice in pathway analysis to use z-scored activities

# Bypass intra = TRUE analysis:
# Tests how well spatial views (neighbors) predict pathway activity 
# WITHOUT the cell's own intrinsic composition as a predictor
# More challenging test - can we predict pathway activity purely from spatial context?




