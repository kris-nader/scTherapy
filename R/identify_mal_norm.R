######################################################################################################
##         Get functions for checkpoint 1: identification of healthy/malignant clusters             ##
######################################################################################################
#
# GNU General Public License v3.0 (https://github.com/kris-nader/TBD/blob/main/LICENSE)
#
# Written by Kristen Michelle Nader <kristen.nader@helsinki.fi>
# and Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi> , April 2021
#
# Functions on this page:
# read_zip,sctype_source,copykat_source,run_scType,run_copyKat,run_SCEVAN,run_ensemble,visualize_ensemble_step
#

#' @title Read data from a zip file
#' @name read_zip
#' @description Read data from a zip file.
#' @param url The URL of the zip file.
#' @return The data file in dataframe format.
#' @export
#' @examples
#' data = read_zip("https://raw.githubusercontent.com/kris-nader/TBD/main/sample_x_exp.rawdata.txt.zip")
#' 

read_zip <- function(url) {
  temp_zip <- tempfile(); temp_dir <- tempdir(); download.file(url, temp_zip)
  file_inside_zip <- unzip(temp_zip, list = TRUE)$Name[1]
  unzip(temp_zip, files = file_inside_zip, exdir = temp_dir)
  data <- suppressWarnings(data.table::fread(file.path(temp_dir, file_inside_zip)) %>% as.data.frame())
  rownames(data) <- data[[1]]; data <- data[,-1]
  return(data)
}




#' @title sctype source files
#' @name sctype_source
#' @description loads sctype functions needed for an automated cell type annotation . 
#' @details none
#' @param none 
#' @return original ScType database
#' @export
#' @examples
#' db_=sctype_source()
#' 

sctype_source <- function(){
  # load tissue auto detect
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  # load ScType database
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  return(db_)
}

#' @title Load modified CopyKat functions
#' @name copykat_source
#' @description loads modified copykat functions needed for an quick classification of aneuploid/diploid cluster . 
#' @details none
#' @param none 
#' @return 
#' @export
#' @examples
#' copykat_source()
#' 

#copykat_source <- function(){
  # load Rcpp helper functions for faster run time
 # #Rcpp::sourceCpp("https://raw.githubusercontent.com/kris-nader/copykat_faster/main/helper_file.cpp")
  # load modified copykat 
 # source("https://raw.githubusercontent.com/kris-nader/copykat_faster/main/faster-copykat.R")
#}

#' @title Load modified scevan functions
#' @name scecan_source
#' @description loads modified scevan functions needed for an quick classification of aneuploid/diploid cluster . 
#' @details none
#' @param none 
#' @return 
#' @export
#' @examples
#' scevan_source()
#' 

scevan_source <- function(){
  # testing for scevan--no plotting 
  # load modified scevan
  source("https://raw.githubusercontent.com/kris-nader/scevan_edit/main/scevan_mod.R")
}



#' @title Run sctype analysis on Seurat object
#' @name run_scType
#' @description run an automated cell type annotation 
#' @details Useful to get an idea of different cells in the sample
#' @details implemented in the ensemble method(1. cell type annotation 2. mal/healthy using cell marker 2.0 --marker based technique) 
#' @details also used for standard cell type annotation (e.g: T cells, B cells...)
#' @param seurat_object A Seurat object
#' @param known_tissue_type The tissue type of the input data (optional)
#' @param custom_marker_file Path to the custom marker file (optional)
#' @param plot Whether to plot the results (default is FALSE)
#' @param name The name of the metadata column to store the scType results (default is "sctype_classification")
#' @return A modified copy of the input Seurat object with a new metadata column
#' 
#' @import sctype source code
#' @import Seurat DimPlot
#' 
#' @examples
#' seurat_object=run_scType(seurat_object,"Immune system)
#' 
#' @export
#' 





## tested -yes
run_sctype <- function(seurat_object, known_tissue_type = NULL, custom_marker_file = NULL, plot = FALSE, name = "sctype_classification") {
  db_=sctype_source()
  # Check for missing arguments
  if (is.null(seurat_object)) {
    stop("Argument 'seurat_object' is missing")
  }
  if (!inherits(seurat_object, "Seurat")) {
    stop("Argument 'seurat_object' must be a Seurat object")
  }
  # Set default custom marker file
  if (is.null(custom_marker_file)) {
    custom_marker_file = db_
  }
  # Auto-detect tissue type if not provided
  if (is.null(known_tissue_type)) {
    tissue_type = auto_detect_tissue_type(path_to_db_file = custom_marker_file, 
                                          seuratObject = seurat_object, 
                                          scaled = TRUE, assay = "RNA")
    rownames(tissue_type)=NULL
    tissue_type=tissue_type$tissue[1]
  } else {
    tissue_type = known_tissue_type
  }
  
  # Prepare gene sets
  gs_list = gene_sets_prepare(custom_marker_file, tissue_type)
  # Calculate scType scores
  es.max = sctype_score(scRNAseqData = seurat_object[["RNA"]]@scale.data,
                        scaled = TRUE,gs = gs_list$gs_positive, 
                        gs2 = gs_list$gs_negative)
  
  # Extract top cell types for each cluster
  cL_resutls = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  seurat_object_res=seurat_object
  seurat_object_res@meta.data[name] = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seurat_object_res@meta.data[seurat_object_res@meta.data$seurat_clusters == j,name] = as.character(cl_type$type[1])
  }
  if(plot){
    plot_ = DimPlot(seurat_object_res, reduction = "umap", group.by = name)   
    print(plot_)
  }
  text_=paste("New metadata added: ",name)
  print(text_)
  return(seurat_object_res)
}        


#' @title Get normal cells as reference
#' @name get_normal_cells
#' @description Extracts cell names that are known to be normal based on a previously done cell type annotation.
#' @param seurat_object Seurat object to analyze.
#' @param names_of_cell_types A character vector of cell type names that are known to be normal.
#' @param column_name The name of the column in the metadata that contains the cell type annotations. Default is "sctype_classification".
#' @return A character vector of cell types that are known to be normal.
#' @import Seurat
#' @examples
#' norm = c("CD8+ NKT-like cells","Naive CD4+ T cells")
#' normal_cells = get_normal_cells(seurat_object, names_of_cell_types = norm)
#' @export

get_normal_cells <- function(seurat_object, names_of_cell_types, column_name = "sctype_classification") {
  cells = rownames(seurat_object@meta.data[seurat_object@meta.data[[column_name]] %in% names_of_cell_types, ])
  return(cells)
}


#' @title Run CopyKat analysis on Seurat object
#' @name run_copyKat
#' @description runs CopyKAT analysis on the gene expression count matrix of a Seurat object and updates the Seurat metadata with the CopyKAT output.
#' @description identifies aneuploid and diploid clusters in a Seurat object using karyotype similarity-based clustering.
#' @details This function is useful for identifying tumor populations and works best when normal/reference cells are provided.
#' @param seurat_object Seurat object to analyze.
#' @param known_normal_cells A character vector of cell names that are known to be normal. Default is NULL (optional)
#' @param plot plot the CopyKAT results on a UMAP plot. Default is False
#' @return A Seurat object with CopyKat results added to metadata.
#' 
#' @import copykat copykat
#' @import Seurat DimPlot
#' 
#' @examples
#' norm = c("CD8+ NKT-like cells","Naive CD4+ T cells")
#' normal_cells = get_normal_cells(seurat_object, names_of_cell_types = norm)
#' seurat_object <- run_copyKat(seurat_object, known_normal_cells = normal_cells, plot = FALSE)
#' 
#' @export
#' 

# tested- yes
run_copyKat <- function(seurat_object, known_normal_cells=NULL, plot=FALSE,ncores = 4,genome=NULL){
  #copykat_source()
  # Extract count matrix
  count_mtx = as.matrix(seurat_object@assays$RNA@counts)
  # Run CopyKAT analysis
  copykat.test = copykatRcpp:::rcpp_copykat(rawmat = count_mtx, id.type = "S", ngene.chr = 5, 
                              win.size = 25, KS.cut = 0.1, sam.name = "test", 
                              distance = "euclidean",norm.cell.names = known_normal_cells,
                              output.seg = FALSE,plot.genes = FALSE, genome = genome,
                              n.cores = ncores)
  # Add annotations -- remove undefined cells
  pred.test = data.frame(copykat.test$prediction)
  pred.test = pred.test[pred.test$copykat.pred != "not.defined", ]  
  ## copykat pred.test rownames may look different from users rownames in the seurat object so we will fix the rownames of the df to match
  rownames(pred.test) = pred.test$cell.names
  # Update metadata with CopyKAT output
  seurat_object_res=seurat_object
  seurat_object_res@meta.data[rownames(seurat_object_res@meta.data),"copyKat_output"] = pred.test[rownames(seurat_object_res@meta.data),"copykat.pred"]
  seurat_object_res@meta.data[which(is.na(seurat_object_res@meta.data$copyKat_output)),"copyKat_output"]="aneuploid"
  seurat_object_res@meta.data[rownames(seurat_object_res@meta.data),"copyKat_output"] = ifelse(seurat_object_res@meta.data[rownames(seurat_object_res@meta.data),"copyKat_output"]=="diploid","healthy","malignant")

  if(plot){
    plot_ = DimPlot(seurat_object_res, reduction = "umap", group.by = 'copyKat_output')    
    print(plot_)
  }
  text_=paste("New metadata added: ","copyKat_output")
  print(text_)
  return(seurat_object_res)
}

## tested- yes
#' @title Run SCEVAN analysis on Seurat object
#' @name run_SCEVAN
#' @description Runs Variational Aneuploidy analysis to identify clusters of aneuploid and diploid cells :  a CNA estimation-based technique, the third approach in the ensemble method.
#' @details This function is used to identify tumor populations and works best when normal/reference cells are provided.
#' @param seurat_object Seurat object to analyze.
#' @param known_normal_cells Names of cells known to be normal--refer to sctype for accurate cell type annotation. Default is NULL.
#' @param plot Whether to plot the results using Seurat's `DimPlot`. Default is FALSE.
#' @return A Seurat object with SCEVAN results added to metadata.
#' 
#' @import SCEVAN pipelineCNA
#' @import Seurat DimPlot
#' 
#' @examples
#' norm = c("CD8+ NKT-like cells","Naive CD4+ T cells")
#' normal_cells = get_normal_cells(seurat_object, names_of_cell_types = norm)
#' seurat_object = run_SCEVAN(seurat_object, known_normal_cells = normal_cells, plot = FALSE)
#' 
#' @export
#' 

run_SCEVAN <- function(seurat_object, known_normal_cells = NULL, plot = FALSE,ncores = 4) {
  scevan_source()
  # Extract count matrix
  count_mtx = seurat_object@assays$RNA@counts
  # Run SCEVAN analysis
  results = pipelineCNA1(count_mtx, sample = "temp", par_cores = ncores, 
                         SUBCLONES = FALSE, plotTree = FALSE, 
                         norm_cell = known_normal_cells)
  # Identify confident normal cells
  confident_normal = rownames(results[which(results$confidentNormal == "yes" ), ])
  # Update metadata with SCEVAN output
  seurat_object_res=seurat_object
  seurat_object_res@meta.data$SCEVAN_output = "malignant"
  seurat_object_res@meta.data[confident_normal, "SCEVAN_output"] = "healthy"
  # Plot results if requested
  if (plot) {
    plot_ = DimPlot(seurat_object_res, reduction = "umap", group.by = "SCEVAN_output")
    print(plot_)
  }
  text_=paste("New metadata added: ","SCEVAN_output")
  print(text_)
  return(seurat_object_res)
}


#' @title Run Ensemble model on seurat object
#' @name run_ensemble
#' @description  This function is the main function for step 1: Identification of malignant/normal clusters
#' run copykat, scevan, and sctype to confidently identify clusters of aneuploid and diploid clusters 
#' a majority vote will be used to finalize classification.
#' @details Classify tumor populations
#' @param seurat_object: Seurat object with raw counts slot
#' @param plot: default=FALSE
#' @param known_normal_cells: default=NULL set of normal cells identified from previous study(see function run_scType)
#' @param disease: default=NULL disease of interest, will be used to identify malignant/healthy cells from new scTypeDOM: sctype disease specific oncogenic markers
#' @return A Seurat object with ensemble results added to metadata.
#' 
#' @import ggplot2
#' @import Seurat DimPlot
#' @import CopyKAT 
#' @import SCEVAN 
#' 
#' @examples
#' norm = c("CD8+ NKT-like cells","Naive CD4+ T cells")
#' normal_cells = get_normal_cells(seurat_object, names_of_cell_types = norm)
#' seurat_object=run_ensemble(seurat_object,known_normal_cells=normal_cells,disease="AML",plot=FALSE)
#' 
#' @export
#' 



run_ensemble <- function(seurat_object, disease=NULL,known_normal_cells=NULL,genome=genome1,plot=FALSE){
  if (is.null(seurat_object)) {
    stop("Argument 'seurat_object' is missing")
  }
  if (!inherits(seurat_object, "Seurat")) {
    stop("Argument 'seurat_object' must be a Seurat object")
  } 
  custom_marker="https://raw.githubusercontent.com/kris-nader/TBD/main/sctype_aml_cellmarker20_cosmic.xlsx";
  
  # run modified sctype-- marker based approach
  seurat_object = run_sctype(seurat_object,known_tissue_type = disease,
                             plot=FALSE,
                             custom_marker_file =custom_marker,
                             name="sctype_malignant_healthy")
  # run copykat analysis-- CNA estimation approach
  seurat_object=run_copyKat(seurat_object,known_normal_cells=known_normal_cells,plot=FALSE,genome="hg20")
  # run SCEVAN analysis-- CNA estimation approach
  seurat_object = run_SCEVAN(seurat_object, 
                             known_normal_cells=known_normal_cells,
                             plot=FALSE)
  seurat_object1=seurat_object
  seurat_object1@meta.data[which(seurat_object1@meta.data$sctype_malignant_healthy=="Unknown"),"sctype_malignant_healthy"]="malignant"
  temp = data.frame(seurat_object1@meta.data$copyKat_output,seurat_object1@meta.data$SCEVAN_output,seurat_object1@meta.data$sctype_malignant_healthy)
  rownames(temp)=rownames(seurat_object1@meta.data)
  # create an ensemble method to take a majority vote
  majority <- apply(temp, 1, function(x) {
    levels_x = unique(x)
    freq_x = tabulate(match(x, levels_x))
    levels_x[which.max(freq_x)]
  })
  seurat_object1@meta.data$ensemble_output = majority
  text_=paste("New metadata added: ","ensemble_output")
  print(text_)
  return(seurat_object1)
}


#' @title Visualize Ensemble step of scRNA-seq analysis
#' @name visualize_ensemble_step
#' @description Generates a plot grid of 5 Seurat DimPlot visualizations with various grouping variables.
#' @param seurat_object a seurat object.
#' @return A plot grid of 5 Seurat DimPlot visualizations.
#'
#' @import Seurat DimPlot
#' @import cowplot plot_grid
#' 
#' @export
#' 
#' @examples
#' visualize_ensemble_step(seurat_object)
#' 
#' 
#' 

visualize_ensemble_step <- function(seurat_object) {
  # Define plot colors
  plot_cols = c("malignant" = "#F8756D", "healthy" = "#03BFC4")
  # Generate Seurat DimPlot visualizations
  p1 = DimPlot(seurat_object, group.by = "ensemble_output", cols = plot_cols,label.size = 3)
  if(identical(integer(0),which(colnames(seurat_object@meta.data)=="sctype_classification"))){
    print("sctype_classification does not exist. For an automated cell type annotation: 'run_scType'")
    p2=NULL
  }
  else{
    p2 = DimPlot(seurat_object, group.by = "sctype_classification", label = TRUE,label.size = 2) + NoLegend()
  }
  p3 = DimPlot(seurat_object, group.by = "sctype_malignant_healthy", cols = plot_cols) + NoLegend()
  p4 = DimPlot(seurat_object, group.by = "copyKat_output", cols = plot_cols) + NoLegend()
  p5 = DimPlot(seurat_object, group.by = "SCEVAN_output", cols = plot_cols) + NoLegend()

  # Combine plots into a grid
  plot_grid(
    p1, plot_grid(p2, p3, nrow = 1), plot_grid(p4, p5, nrow = 1),
    nrow = 3, rel_heights = c(1, 0.8, 0.8),
    labels = c("A", "B", "C"), label_size = 15
  )
}
