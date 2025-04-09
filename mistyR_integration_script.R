data_path <- "C:/Users/juesh/UROP/data/misty"

results_folder <- "C:/Users/juesh/UROP/results/misty"


python_path <- NULL 


instructions <- createGiottoInstructions(save_dir = results_folder, 
                                         save_plot = TRUE, 
                                         show_plot = FALSE, 
                                         return_plot = FALSE, 
                                         python_path = python_path)

visium_lungcancer <- createGiottoVisiumObject(visium_dir = data_path,
                                              expr_data = "raw",
                                              png_name = "tissue_lowres_image.png",
                                              gene_column_index = 2,
                                              instructions = instructions)


pDataDT(visium_lungcancer)

# check available image names
showGiottoImageNames(visium_lungcancer) # "image" is the default name

# show aligned image
spatPlot(gobject = visium_lungcancer, 
         cell_color = "in_tissue", 
         show_image = TRUE, 
         point_alpha = 0.7)



visium_lungcancer <- filterGiotto(gobject = visium_lungcancer,
                                  expression_threshold = 1,
                                  feat_det_in_min_cells = 20,
                                  min_det_feats_per_cell = 500,
                                  expression_values = "raw",
                                  verbose = TRUE)

visium_lungcancer <- normalizeGiotto(gobject = visium_lungcancer, 
                                     scalefactor = 6000, 
                                     verbose = TRUE)

visium_lungcancer <- addStatistics(gobject = visium_lungcancer)



visium_lungcancer <- calculateHVF(gobject = visium_lungcancer)

visium_lungcancer <- runPCA(gobject = visium_lungcancer)

screePlot(visium_lungcancer, 
          ncp = 30)
visium_lungcancer <- runUMAP(visium_lungcancer, 
                             dimensions_to_use = 1:10)

plotUMAP(gobject = visium_lungcancer)


visium_lungcancer <- runtSNE(visium_lungcancer, 
                             dimensions_to_use = 1:10)

plotTSNE(gobject = visium_lungcancer)


# Create shared nearest network (SNN) and perform leiden clustering
visium_lungcancer <- createNearestNetwork(gobject = visium_lungcancer, 
                                          dimensions_to_use = 1:10, 
                                          k = 30)

visium_lungcancer <- doLeidenCluster(gobject = visium_lungcancer, 
                                     spat_unit = "cell", 
                                     feat_type = "rna", 
                                     resolution = 0.4, 
                                     n_iterations = 1000)

# visualize UMAP cluster results
plotUMAP(gobject = visium_lungcancer, 
         cell_color = "leiden_clus", 
         show_NN_network = TRUE, 
         point_size = 2)

# visualize tSNE cluster results
plotTSNE(gobject = visium_lungcancer, 
         cell_color = "leiden_clus", 
         show_NN_network = TRUE, 
         point_size = 2)

# visualize expression and spatial results
spatDimPlot(gobject = visium_lungcancer, 
            cell_color = "leiden_clus",
            dim_point_size = 2, 
            spat_point_size = 2)

spatDimPlot(gobject = visium_lungcancer, 
            cell_color = "nr_feats", 
            color_as_factor = FALSE,
            dim_point_size = 2, 
            dim_show_legend = TRUE, 
            spat_show_legend = TRUE,
            spat_point_size = 2)


# Cell type marker detection
# Gini markers
gini_markers_subclusters <- findMarkers_one_vs_all(gobject = visium_lungcancer,
                                                   method = "gini",
                                                   expression_values = "normalized",
                                                   cluster_column = "leiden_clus",
                                                   min_featss = 20,
                                                   min_expr_gini_score = 0.5,
                                                   min_det_gini_score = 0.5)

# get top 2 genes per cluster and visualize with violin plot
topgenes_gini <- gini_markers_subclusters[, head(.SD, 2), by = "cluster"]$feats

violinPlot(visium_lungcancer, 
           feats = unique(topgenes_gini), 
           cluster_column = "leiden_clus",
           strip_text = 8, 
           strip_position = "right")


# cluster heatmap
plotMetaDataHeatmap(visium_lungcancer,
                    selected_feats = topgenes_gini,
                    metadata_cols = "leiden_clus",
                    x_text_size = 10, 
                    y_text_size = 10)

# umap plots
dimFeatPlot2D(visium_lungcancer,
              expression_values = "scaled",
              feats = gini_markers_subclusters[, head(.SD, 1), by = "cluster"]$feats,
              cow_n_col = 3,
              point_size = 1)

# Cell type marker detection
# Scran markers
scran_markers_subclusters <- findMarkers_one_vs_all(gobject = visium_lungcancer,
                                                    method = "scran",
                                                    expression_values = "normalized",
                                                    cluster_column = "leiden_clus")

# get top 2 genes per cluster and visualize with violin plot
topgenes_scran <- scran_markers_subclusters[, head(.SD, 2), by = "cluster"]$feats

violinPlot(visium_lungcancer, 
           feats = unique(topgenes_scran),
           cluster_column = "leiden_clus",
           strip_text = 10, 
           strip_position = "right")


# cluster heatmap
plotMetaDataHeatmap(visium_lungcancer,
                    selected_feats = topgenes_scran,
                    metadata_cols = "leiden_clus",
                    x_text_size = 10, 
                    y_text_size = 10)



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

visium_lungcancer <- runPAGEEnrich(gobject = visium_lungcancer, 
                                   sign_matrix = signature_matrix, 
                                   min_overlap_genes = 1)

cell_types <- colnames(signature_matrix)

cell_types_subset1 <- cell_types[1:8]  # First set of cell types
cell_types_subset2 <- cell_types[9:16]  # Second set of cell types

plotMetaDataCellsHeatmap(gobject = visium_lungcancer,
                         metadata_cols = "leiden_clus",
                         value_cols = cell_types,
                         spat_enr_names = "PAGE",
                         x_text_size = 8,
                         y_text_size = 8,
                         show_plot = TRUE)


# Visualize first set of cell types
spatCellPlot(gobject = visium_lungcancer,
             spat_enr_names = "PAGE",
             cell_annotation_values = cell_types_subset1,
             cow_n_col = 4, 
             coord_fix_ratio = NULL, 
             point_size = 0.75)

# Visualize second set of cell types
spatCellPlot(gobject = visium_lungcancer,
             spat_enr_names = "PAGE",
             cell_annotation_values = cell_types_subset2,
             cow_n_col = 4, 
             coord_fix_ratio = NULL, 
             point_size = 0.75)



# For example, visualizing fibroblasts, epithelial cells, and cancer cells:
spatDimCellPlot(gobject = visium_lungcancer,
                spat_enr_names = "PAGE",
                cell_annotation_values = c("Fibroblast", "Alveolar_TypeII", "Squamous_Carcinoma"),
                cow_n_col = 1, 
                spat_point_size = 1.2,
                plot_alignment = "horizontal")


# 1.3 enrichment test with PAGE
markers_scran <- findMarkers_one_vs_all(gobject = giotto_SC, 
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


visium_lungcancer <- createSpatialGrid(gobject = visium_lungcancer,
                                       sdimx_stepsize = 400,
                                       sdimy_stepsize = 400,
                                       minimum_padding = 0)

spatPlot(visium_lungcancer, 
         cell_color = "leiden_clus", 
         point_size = 2.5, 
         show_grid = TRUE,
         grid_color = "red", 
         spatial_grid_name = "spatial_grid")


## Delaunay network: stats + creation
plotStatDelaunayNetwork(gobject = visium_lungcancer, 
                        maximum_distance = 400)
         
visium_lungcancer <- createSpatialNetwork(gobject = visium_lungcancer, 
                                          minimum_k = 0)

showGiottoSpatNetworks(visium_lungcancer)

spatPlot(gobject = visium_lungcancer, 
         show_network = TRUE,
         network_color = "blue", 
         spatial_network_name = "Delaunay_network")

# kmeans binarization
km_spatialfeats <- binSpect(visium_lungcancer)

spatFeatPlot2D(visium_lungcancer, 
               expression_values = "scaled",
               feats = km_spatialfeats$feats[1:6], 
               cow_n_col = 2, 
               point_size = 1.5)


saveGiotto(visium_lungcancer, "results/misty/visium_lungcancer_object", overwrite = TRUE)

#visium_lungcancer <- loadGiotto("results/misty/visium_lungcancer_object")







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
library(tibble)
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

# SpatialExperiment
library(SpatialExperiment)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(GSVA)
library(msigdbr)



#---------------------------------------------
#LOAD GIOTTO OBJECT

visium_lungcancer_object <- loadGiotto(path_to_folder = "C:/Users/juesh/UROP/results/misty/visium_lungcancer_object")





#---------------------------------------------
#EXTRACT AND PREPARE EXPRESSION DATA:

raw_matrix <- visium_lungcancer_object@expression$cell$rna$raw@exprMat
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
geometry <- getSpatialLocations(visium_lungcancer, output = "data.table")
geometry <- geometry[, c("sdimx", "sdimy")]
setnames(geometry, c("sdimx", "sdimy"), c("row", "col"))





#---------------------------------------------
#PATHWAY ANALYSIS WITH GSVA:

# Transpose the signature matrix (cell types as rows, genes as columns)
signature_matrix_transposed <- t(signature_matrix)

# Access the PAGE element 
page_results <- visium_lungcancer@spatial_enrichment$cell$rna$PAGE
page_data <- page_results@enrichDT

# Display the structure of the PAGE results
str(page_results, max.level = 2)

# Convert to data frame and then tibble
page_data_df <- as.data.frame(page_data)
composition <- tibble::as_tibble(page_data_df)

# Clean column names
colnames(composition) <- composition %>% 
  clean_names(parsing_option = 0) %>% 
  colnames(.)

# Get MSigDB Hallmark pathways for human
msigdb_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")

# Convert MSigDB data to list format for GSVA
gene_sets_list <- split(msigdb_gene_sets$gene_symbol, msigdb_gene_sets$gs_name)

# Check overlap between expression genes and MSigDB genes
all_msigdb_genes <- unique(msigdb_gene_sets$gene_symbol)
overlap_count <- length(intersect(rownames(expression), all_msigdb_genes))
print(paste("Number of overlapping genes:", overlap_count))

expression_matrix <- as.matrix(expression)

# Create the parameter object with the correct constructor
gsva_param <- GSVA::gsvaParam(
  exprData = expression_matrix, 
  geneSets = gene_sets_list,
  minSize = 5,
  maxSize = 500,
  kcdf = "Gaussian"  # For continuous data like normalized expression
)

# Run GSVA with the parameter object
gsva_result <- GSVA::gsva(param = gsva_param, verbose = TRUE)


# Check the dimensions of your GSVA result
dim(gsva_result)

# View the first few rows and columns
gsva_result[1:5, 1:5]

# Convert to a data frame or tibble for easier manipulation
gsva_df <- as.data.frame(gsva_result)

# Convert GSVA results to a format compatible with misty 
pathway_scores <- t(gsva_result) 

# Create a tibble for MISTy analysis
pathway_view_data <- as_tibble(pathway_scores)





#---------------------------------------------
#PREP DATA FOR MISTY ANALYSIS:

# order of cells in pathway_view_data matches expression data
if(length(rownames(pathway_view_data)) == 0) {
      rownames(pathway_view_data) <- colnames(expression)
}

#needed for misty: expression, spatial barcodes, page result barcodes (composition)

if(nrow(geometry) == length(expression_barcodes)) {
  geometry$full_barcode <- expression_barcodes
} else {
  warning("Number of barcodes doesn't match number of rows in geometry.")
}

expression_barcodes <- colnames(expression)
composition_barcodes <- page_data$cell_ID
spatial_barcodes <- geometry$full_barcode


pathway_scores_fixed <- t(gsva_result)
rownames(pathway_scores_fixed) <- colnames(expression)
pathway_view_data_new <- as_tibble(pathway_scores_fixed, rownames = "cell_ID")

# Convert pathway_view_data_new from tibble to dataframe
pathway_df <- as.data.frame(pathway_view_data_new)
rownames(pathway_df) <- pathway_df$cell_ID
pathway_df$cell_ID <- NULL  # Remove the cell_ID column

expression_barcodes <- unname(expression_barcodes)

# Recalculate common_cells (should be all cells since all barcodes match)
common_cells <- expression_barcodes

# Now filter pathway_df using the common_cells
pathway_df_filtered <- pathway_df[common_cells, ]
print(paste("Dimensions of pathway_df_filtered:", paste(dim(pathway_df_filtered), collapse="×")))

# Create geometry for MISTy
geometry_filtered <- data.frame(
  row = geometry$row,
  col = geometry$col,
  row.names = spatial_barcodes
)





#---------------------------------------------
#CREATE MISTY VIEWS

comp_views <- create_initial_view(page_data[, c("Alveolar_TypeI", "Alveolar_TypeII", "B_cell", 
                                                "Bcell_PB", "Bcells_memory", "Ciliated", 
                                                "Dendritic", "Endothelial", "Fibroblast", 
                                                "Macrophage", "Mast", "Neutrophils", 
                                                "Smooth_Muscle", "Squamous_Carcinoma", 
                                                "T_Cells", "pDCs")])
#add spatial views
comp_views <- comp_views %>%
  add_juxtaview(geometry_filtered, neighbor.thr = 55) %>%
  add_paraview(geometry_filtered, l = 200, family = "gaussian")

# Combine with pathway views
comp_views <- comp_views %>%
  add_views(create_view("juxtaview.path.55", path_act_views[["juxtaview.55"]]$data, "juxta.path.55")) %>% 
  add_views(create_view("paraview.path.200", path_act_views[["paraview.200"]]$data, "para.path.200"))





#---------------------------------------------
# RUN MISTY ANALYSIS:

run_misty(comp_views, "results/misty/comp_path_act")

# Collect and visualize results
misty_results <- collect_results(paste0(results_folder, "/comp_path_act"))





#---------------------------------------------
#VISUALIZE RESULTS:

pdf("C:/Users/juesh/UROP/results/misty/misty_plots/misty_plots.pdf", width = 12, height = 10)

# Plot improvement statistics
misty_results %>%
  plot_improvement_stats("intra.R2") %>%
  plot_improvement_stats("gain.R2")


# Plot view contributions
misty_results %>%
  plot_view_contributions()

#plot interaction heatmaps
misty_results %>%
  plot_interaction_heatmap("juxta.path.55", clean = TRUE)

misty_results %>%
  plot_interaction_heatmap("para.path.200", clean = TRUE)

#close PDF
dev.off()








#---------------------------------------------
#TF Activity Estimation



py_install(c("decoupler", "omnipath"), pip = TRUE)

# Load Python module
dc <- import("decoupler")

# Convert the sparse matrix to a dense matrix
dense_expression <- as.matrix(expression)

# Transpose the dense matrix
transposed_expression <- t(dense_expression)

# Convert the transposed matrix to a data.frame
expression_df <- as.data.frame(transposed_expression)

# Get CollecTRI network
net <- dc$get_collectri()

# Run ULM to estimate TF activities
acts_tfs <- dc$run_ulm(
  mat = expression_df,
  net = net,
  verbose = TRUE,
  use_raw = FALSE
)

# Convert back to R
est_TF <- py_to_r(acts_tfs)

# Convert expression to a sparse matrix
#temp_expr_matrix <- Matrix::Matrix(as.matrix(expression), sparse = TRUE)







#---------------------------------------------
#Identify highly variable genes and TFs 


# Calculate variance for each gene
gene_var <- apply(as.matrix(expression), 1, var)

# Sort genes by variance
sorted_var <- sort(gene_var, decreasing = TRUE)

# Take top 1000 variable genes
hvg_names <- names(sorted_var)[1:1000]

# Extract expression data for these HVGs
hvg_expr <- expression[hvg_names, ]

# Extract TFs that are also in the highly variable genes
hvg_TF <- est_TF[[1]][, colnames(est_TF[[1]]) %in% rownames(hvg_expr)]






#---------------------------------------------
# Create TF view and integrate with other views


# Create TF view with paraview
TF_view <- mistyR::create_initial_view(hvg_TF)
TF_view <- mistyR::add_paraview(TF_view, geometry_filtered, l = 200, family = "gaussian")

# Combine Views (cell type composition, TF, and pathway views)
comp_TF_path_views <- comp_views %>% 
  add_views(create_view("paraview.TF.200", TF_view[["paraview.200"]]$data, "para.TF.200"))

# Run MISTy
run_misty(comp_TF_path_views, paste0(results_folder, "/comp_TF_path"), 
          model.function = linear_model, bypass.intra = TRUE)

# Collect results
misty_results_comp_TF_pathway <- collect_results(paste0(results_folder, "/comp_TF_path"))

# Visualize results
pdf(paste0(results_folder, "/misty_plots/TF_pathway_interactions.pdf"), width = 12, height = 10)
misty_results_comp_TF_pathway %>%
  plot_improvement_stats("gain.R2") %>%
  plot_view_contributions()

misty_results_comp_TF_pathway %>%
  plot_interaction_heatmap("para.TF.200", 
                           clean = TRUE,
                           trim.measure = "gain.R2",
                           trim = 20)
dev.off()






#---------------------------------------------
# Ligand-Receptor Analysis


# Download OmniPath resource
download.file("https://raw.githubusercontent.com/saezlab/liana-py/main/liana/resource/omni_resource.csv", 
              destfile = paste0(results_folder, "/omni_resource.csv"), method = "curl")

# Load Ligand Receptor Resource
omni_resource <- read_csv(paste0(results_folder, "/omni_resource.csv")) %>% 
  filter(resource == "consensus")

# Get highly variable ligands
ligands <- omni_resource %>% 
  pull(source_genesymbol) %>% 
  unique()
hvg_lig <- hvg_expr[rownames(hvg_expr) %in% ligands,]

# Get highly variable receptors
receptors <- omni_resource %>% 
  pull(target_genesymbol) %>% 
  unique()
hvg_recep <- hvg_expr[rownames(hvg_expr) %in% receptors,]

# Clean names
rownames(hvg_lig) <- hvg_lig %>% 
  clean_names(parsing_option = 0) %>% 
  rownames(.)

rownames(hvg_recep) <- hvg_recep %>% 
  clean_names(parsing_option = 0) %>% 
  rownames(.)


# Create views for ligands and receptors
#hvg_recep is a sparse matrix of class dgtMatrix so we have to 
#convert it to a dense matrix before transposing and converting

receptor_view <- create_initial_view(as.data.frame(as.matrix(t(hvg_recep))))

ligand_view <- create_initial_view(as.data.frame(as.matrix(t(hvg_lig)))) %>% 
  add_paraview(geometry_filtered, l = 200, family = "gaussian")

# Combine views
lig_recep_view <- receptor_view %>% 
  add_views(create_view("paraview.ligand.200", ligand_view[["paraview.200"]]$data, "para.lig.200"))

# Run MISTy
run_misty(lig_recep_view, paste0(results_folder, "/lig_recep"), bypass.intra = TRUE)

# Collect and visualize results
misty_results_lig_recep <- collect_results(paste0(results_folder, "/lig_recep"))

pdf(paste0(results_folder, "/misty_plots/ligand_receptor_interactions.pdf"), width = 12, height = 10)
misty_results_lig_recep %>%
  plot_interaction_heatmap("para.lig.200", clean = TRUE, cutoff = 2, 
                           trim.measure = "gain.R2", trim = 10)
dev.off()


saveGiotto(visium_lungcancer, "results/misty/visium_lungcancer_object", overwrite = TRUE)





#---------------------------------------------

#VISUALIZE HALLMARK PATHWAYS SPATIALLY:


# Calculate variance for each hallmark pathway
hallmark_variance <- apply(gsva_result, 1, var)
top_hallmarks <- names(sort(hallmark_variance, decreasing = TRUE))[1:12]

# In Giotto 4.2.1, we need to use addGeneMetadata to add custom features
# Create a data.frame with cell IDs and hallmark scores
all_cell_ids <- pDataDT(visium_lungcancer)$cell_ID

# Start PDF device for hallmark spatial plots
pdf("results/misty/hallmarks_spatial_plots.pdf", width = 12, height = 10)

# First, create a spatial plot showing leiden clusters for reference
spatPlot(gobject = visium_lungcancer, 
         cell_color = "leiden_clus", 
         point_size = 2.2,
         title = "Leiden Clusters",
         show_image = TRUE)


# For each hallmark, add it to the Giotto object as cell metadata and create a plot
for (hallmark in top_hallmarks) {
  # Extract values for this hallmark
  cell_values <- gsva_result[hallmark, ]
  
  # Create data.frame for cell metadata
  cell_metadata <- data.frame(
    cell_ID = all_cell_ids,
    stringsAsFactors = FALSE
  )
  cell_metadata[[hallmark]] <- cell_values[match(all_cell_ids, colnames(gsva_result))]
  
  # Add cell metadata to Giotto object
  visium_lungcancer <- addCellMetadata(visium_lungcancer, 
                                       new_metadata = cell_metadata, 
                                       by_column = TRUE,
                                       column_cell_ID = "cell_ID")
  
  # Create spatial plot for this hallmark
  spatPlot(gobject = visium_lungcancer, 
           cell_color = hallmark, 
           color_as_factor = FALSE,
           point_size = 2.2,
           title = hallmark,
           show_image = TRUE)
  
}

# Close the PDF device
dev.off()


# Create a second PDF with multiple hallmarks in one plot
pdf("results/misty/multi_hallmarks_plots.pdf", width = 12, height = 10)

# Use spatPlot with cowplot for multiple plots
for (i in seq(1, length(top_hallmarks), by = 4)) {
  # Get a subset of hallmarks (up to 4 at a time)
  end_idx <- min(i+3, length(top_hallmarks))
  subset_hallmarks <- top_hallmarks[i:end_idx]
  
  # Create a list to store plots
  plot_list <- list()
  
  # Generate plots for each hallmark in the subset
  for (j in 1:length(subset_hallmarks)) {
    hallmark <- subset_hallmarks[j]
    p <- spatPlot(gobject = visium_lungcancer, 
                  cell_color = hallmark, 
                  color_as_factor = FALSE,
                  point_size = 2.2,
                  title = hallmark,
                  show_image = TRUE,
                  return_plot = TRUE)
    plot_list[[j]] <- p
  }
  
  # Combine plots with cowplot
  library(cowplot)
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)
  print(combined_plot)
}

# create plots for the cell types from PAGE analysis
cell_types <- c("Alveolar_TypeI", "Alveolar_TypeII", "B_cell", 
                "Bcell_PB", "Bcells_memory", "Ciliated", 
                "Dendritic", "Endothelial", "Fibroblast", 
                "Macrophage", "Mast", "Neutrophils", 
                "Smooth_Muscle", "Squamous_Carcinoma", 
                "T_Cells", "pDCs")

# Create spatial cell plots for each cell type group (4 at a time)
for (i in seq(1, length(cell_types), by = 4)) {
  end_idx <- min(i+3, length(cell_types))
  subset_cell_types <- cell_types[i:end_idx]
  
  spatCellPlot(gobject = visium_lungcancer,
               spat_enr_names = "PAGE",
               cell_annotation_values = subset_cell_types,
               cow_n_col = 2,
               point_size = 2.2)
}

dev.off()


#---------------------------------------------

#SAVE UPDATED GIOTTO OBJECT 

final_visium_lungcancer <- visium_lungcancer
saveGiotto(final_visium_lungcancer, "results/misty/final_visium_lungcancer_object", overwrite = TRUE)

