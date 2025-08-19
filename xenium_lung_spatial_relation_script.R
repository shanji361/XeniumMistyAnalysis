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


# set up paths
data_path <- "/Users/juesh/UROP/data/xenium_lung"
save_dir <- '/Users/juesh/UROP/results/xenium_lung'
dir.create(save_dir, recursive = TRUE)

# download the mini dataset and untar
options("timeout" = Inf)
download.file(
  url = "https://zenodo.org/records/13207308/files/workshop_xenium.zip?download=1",
  destfile = file.path(save_dir, "workshop_xenium.zip"),
  mode = "wb",           # Force binary mode
  method = "libcurl",    # Use libcurl method
  timeout = 300          # Set longer timeout
)



unzip(file.path(save_dir, "workshop_xenium.zip"), 
      exdir = data_path)



g <- createGiottoXeniumObject(xenium_dir = data_path)

# set instructions for save directory and to save the plots to disk
instructions(g, "save_dir") <- save_dir
instructions(g, "save_plot") <- TRUE



g <- addSpatialCentroidLocations(g, poly_info = "cell")
g <- addSpatialCentroidLocations(g, poly_info = "nucleus")

spatInSituPlotPoints(g,
                     polygon_feat_type = "cell",
                     feats = list(rna = head(featIDs(g))), # must be named list
                     use_overlap = FALSE, 
                     polygon_color = "cyan", 
                     polygon_line_size = 0.1
)



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



# example plots
spatInSituPlotPoints(g,
                     show_image = TRUE,
                     image_name = "HE",
                     polygon_feat_type = "cell",
                     polygon_color = "cyan",
                     polygon_line_size = 0.1,
                     polygon_alpha = 0
)




spatInSituPlotPoints(g,
                     show_image = TRUE,
                     image_name = "DAPI",
                     polygon_feat_type = "nucleus",
                     polygon_color = "cyan",
                     polygon_line_size = 0.1,
                     polygon_alpha = 0
)


spatInSituPlotPoints(g,
                     show_image = TRUE,
                     image_name = "18S",
                     polygon_feat_type = "cell",
                     polygon_color = "cyan",
                     polygon_line_size = 0.1,
                     polygon_alpha = 0
)



spatInSituPlotPoints(g,
                     show_image = TRUE,
                     image_name = "ATP1A1/CD45/E-Cadherin",
                     polygon_feat_type = "nucleus",
                     polygon_color = "cyan",
                     polygon_line_size = 0.1,
                     polygon_alpha = 0
)


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

spatInSituPlotPoints(g,
                     polygon_fill = "nr_feats",
                     polygon_fill_gradient_style = "sequential",
                     polygon_fill_as_factor = FALSE
)


spatInSituPlotPoints(g,
                     polygon_fill = "total_expr",
                     polygon_fill_gradient_style = "sequential",
                     polygon_fill_as_factor = FALSE
)

g <- runPCA(g, feats_to_use = NULL)


screePlot(g, ncp = 30)


g <- runUMAP(g, 
             dimensions_to_use = seq(15), 
             n_neighbors = 40 # default
)

plotPCA(g)
plotUMAP(g)


g <- createNearestNetwork(g,
                          dimensions_to_use = seq(15), 
                          k = 40
)


g <- doLeidenCluster(g)


plotPCA_3D(g, 
           cell_color = "leiden_clus", 
           point_size = 1
)


plotUMAP(g, 
         cell_color = "leiden_clus", 
         point_size = 0.1, 
         point_shape = "no_border"
)

spatInSituPlotPoints(g,
                     polygon_fill = "leiden_clus",
                     polygon_fill_as_factor = TRUE,
                     polygon_alpha = 1,
                     show_image = TRUE,
                     image_name = "HE"
)

# scran marker gene
res_scran <- findMarkers_one_vs_all(g, 
                                    cluster_column = "leiden_clus", 
                                    method = "scran",
                                    expression_values = "normalized"
)

# top 2 genes per cluster
topgenes_scran <- res_scran[, head(.SD, 2), by = 'cluster']



violinPlot(g, 
           feats = unique(rankgenes_scran$feats), 
           cluster_column = "leiden_clus", 
           save_param = list(base_height = 20,base_width = 10)
)


# writeChatGPTquery(
#   DEG_output = res_scran,
#   top_n_genes = 25,
#   tissue_type = 'lung_cancer',
#   folder_name = "C:/Users/juesh/UROP/results/xenium_lung/",
#   file_name = 'GPTQuery_Scran.txt'
# )


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
  save_param = list(base_height = 8,base_width = 5)
)



cell_types <- c(
  "CD8+ T cell",      # IL7R, TRAC, CD3E, CD8A, GZMA, PRDM1
  "Epithelial cell",  # EPCAM, KRT7, GPRC5A, CYP2B6, CFTR
  "Fibroblast",       # PDGFRA, FBN1, LTBP2, THY1, VCAN
  "Macrophage",       # CD68, CD14, MS4A6A, CD163, AIF1
  "B cell",           # MS4A1, CD19, CD79A, BANK1, IRF8
  "Ciliated cell",    # FOXJ1, DNAAF1, CFAP53, CCDC39, SOX2
  "Endothelial cell", # VWF, PECAM1, CD34, EGFL7, ADGRL4
  "Club cell",        # SCGB1A1, GPRC5A, CYP2F1, ADAM28, AGR3, EHF
  "Epithelial cell 2", # EPCAM, KRT7, GPRC5A, CFTR, CYP2B6
  "Smooth muscle cell", # ACTA2, MYH11, CNN1, MYLK, DES
  "Plasma cell",      # MZB1, SLAMF7, PRDM1, CD27, TNFRSF17
  "Mast cell"         # KIT, CPA3, MS4A2, HPGDS, IL1RL1
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


raw_matrix <- xenium_lungcancer_test@expression$cell$rna$raw[]
symbols <- rownames(raw_matrix)

d.index <- which(duplicated(symbols))
print(d.index) 
#no duplicates found 



geometry <- getSpatialLocations(
  gobject = xenium_lungcancer_test, 
  spat_unit = "cell",  # or "aggregate"
  output = "data.table"
)
geometry <- geometry[, c("sdimx", "sdimy")]
setnames(geometry, c("sdimx", "sdimy"), c("row", "col"))




expression <- as.matrix(xenium_lungcancer_test@expression$cell$rna$raw[])

model <- get_progeny(organism = "human", top = 500)

# Use multivariate linear model to estimate activity
est_path_act <- run_mlm(expression, model,.mor = NULL) 


est_path_act_wide <- est_path_act %>% 
  pivot_wider(id_cols = condition, names_from = source, values_from = score) %>%
  column_to_rownames("condition") 

colnames(est_path_act_wide) <- est_path_act_wide %>% 
  clean_names(parsing_option = 0) %>% 
  colnames(.)

# Add progeny results to object
xenium_lungcancer_test@expression[["cell"]][["rna"]][["progeny"]] <- t(est_path_act_wide)



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





# Create intraview from cell-type identity 
comp_views <- create_initial_view(composition_xenium) 

# Juxta & para from pathway activity 
path_act_views <- create_initial_view(est_path_act_wide) %>%
  add_juxtaview(geometry, neighbor.thr = 20) %>% 
  add_paraview(geometry, l = 50, family = "gaussian")


com_path_act_views <- comp_views %>%
  add_views(create_view("juxtaview.path.20", 
                        path_act_views[["juxtaview.20"]]$data, 
                        "juxta.path.20")) %>% 
  add_views(create_view("paraview.path.50", 
                        path_act_views[["paraview.50"]]$data, 
                        "para.path.50"))




run_misty(com_path_act_views, "result/xenium_lung/comp_path_act")

misty_results_com_path_act <- collect_results("result/xenium_lung/comp_path_act/")


misty_results_com_path_act %>%
  plot_improvement_stats("intra.R2")%>%
  plot_improvement_stats("gain.R2") 


misty_results_com_path_act %>% 
  plot_view_contributions()

misty_results_com_path_act %>%
  plot_interaction_heatmap("juxta.path.20", clean = TRUE)




xenium_lungcancer_test@expression[["cell"]][["progeny"]][["raw"]] <- t(est_path_act_wide)

# Create normalized version (z-score normalization)
est_path_act_normalized <- scale(est_path_act_wide)
xenium_lungcancer_test@expression[["cell"]][["progeny"]][["normalized"]] <- t(est_path_act_normalized)




spatPlot2D(xenium_lungcancer_test,
           spat_unit = "cell", 
           cell_color = "subannot_clus",
           show_image = TRUE,
           point_size = 1.1,
           point_alpha = 0.8, 
           cell_color_gradient =  c("blue", "orange"))




#bypass intra = TRUE:
#see how well other views explained the intraview without intra's self prediction 
run_misty(com_path_act_views, "result/xenium_lung", model.function = linear_model, bypass.intra = TRUE)

misty_results_com_path_act_linear <- collect_results("result/xenium_lung")

misty_results_com_path_act_linear %>%
  plot_improvement_stats("gain.R2") %>%
  plot_view_contributions()


misty_results_com_path_act_linear %>%
  plot_interaction_heatmap("juxta.path.20", clean = TRUE) 


#default to normalized progeny 
spatFeatPlot2D(xenium_lungcancer_test,
               spat_unit = "cell", 
               feat_type = "progeny",
               feats = c("trail"),
               show_image = TRUE,
               point_size = 1.3,
               cow_n_col = 1, 
               cell_color_gradient =  c("blue", "red"))

spatPlot2D(xenium_lungcancer_test,
                       spat_unit = "cell",
                       cell_color = "b.cell",
                       show_image = TRUE,
                       point_size = 1.3)


spatFeatPlot2D(xenium_lungcancer_test,
               spat_unit = "cell", 
               feat_type = "progeny",
               feats = c("vegf", "wnt"),
               show_image = TRUE,
               point_size = 0.8,
               cow_n_col = 2,  
               cell_color_gradient =  c("blue", "green"))


spatFeatPlot2D(xenium_lungcancer_test,
               spat_unit = "cell", 
               feat_type = "progeny",
               feats = c("tnfa", "trail"),
               show_image = TRUE,
               point_size = 0.8,
               cow_n_col = 2,  
               cell_color_gradient =  c("blue", "green"))




spatFeatPlot2D(xenium_lungcancer_test,
               spat_unit = "cell", 
               feat_type = "progeny",
               feats = c( "tgfb", "androgen"),
               show_image = TRUE,
               point_size = 0.8,
               cow_n_col = 2,  
               cell_color_gradient =  c("blue", "green"))

spatFeatPlot2D(xenium_lungcancer_test,
               spat_unit = "cell", 
               feat_type = "progeny",
               feats = c( "nfkb", "p53"),
               show_image = TRUE,
               point_size = 0.8,
               cow_n_col = 2,  
               cell_color_gradient =  c("blue", "green"))


spatFeatPlot2D(xenium_lungcancer_test,
               spat_unit = "cell", 
               feat_type = "progeny",
               feats = c( "egfr", "estrogen"),
               show_image = TRUE,
               point_size = 0.8,
               cow_n_col = 2,  
               cell_color_gradient =  c("blue", "green"))




