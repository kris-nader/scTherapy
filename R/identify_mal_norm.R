######################################################################################################
##         Get functions for checkpoint 1: identification of healthy/malignant clusters             ##
######################################################################################################
#
# GNU General Public License v3.0 (https://github.com/kris-nader/scTherapy/blob/main/LICENSE)
#
# Written by Kristen Michelle Nader <kristen.nader@helsinki.fi>
# and Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi> , April 2021
#
# Functions on this page:
# read_zip,sctype_source,copykat_source,run_scType,run_copyKat,run_SCEVAN,run_ensemble,visualize_ensemble_step
# og_copykat

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
run_copyKat <- function(seurat_object, known_normal_cells="", plot=FALSE,ncores = 4,genome=NULL){
  #source("https://raw.githubusercontent.com/kris-nader/cpp_copykat/main/copykat_noplot_original.R");
  # Extract count matrix
  count_mtx = as.matrix(seurat_object@assays$RNA@counts)
  # Run CopyKAT analysis
  copykat.test = og_copykat(rawmat = count_mtx, id.type = "S", ngene.chr = 5, 
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
  check= openxlsx::read.xlsx(custom_marker)
  sctype_e=FALSE
  if ( disease %in% unique(check$tissueType) {
    sctype_e=TRUE
    seurat_object = run_sctype(seurat_object,known_tissue_type = disease,
                             plot=FALSE,
                             custom_marker_file =custom_marker,
                             name="sctype_malignant_healthy")
    }
      
  # run copykat analysis-- CNA estimation approach
  # Determine the number of cores to use based on OS
  ncoree <- ifelse(Sys.info()["sysname"] != "Windows", max(1, parallel::detectCores() - 2), 1)
  seurat_object=run_copyKat(seurat_object,known_normal_cells=known_normal_cells,plot=FALSE,genome="hg20",ncores=ncoree)
  # run SCEVAN analysis-- CNA estimation approach
  seurat_object = run_SCEVAN(seurat_object, 
                             known_normal_cells=known_normal_cells,
                             plot=FALSE)
  seurat_object1=seurat_object
  if (sctype_e ==TRUE){
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
    }
      else{
        temp = data.frame(seurat_object1@meta.data$copyKat_output,seurat_object1@meta.data$SCEVAN_output)
        rownames(temp)=rownames(seurat_object1@meta.data)
        temp$ensemble_output <- ifelse(temp$SCEVAN_output == temp$copyKat_output, temp$copyKat_output, "malignant")
        seurat_object1@meta.data$ensemble_output[rownames(temp)] = temp[rownames(temp),"ensemble_output"]
        text_=paste("New metadata added: ","ensemble_output")
        print(text_)
        }
      
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


# copykat without plotting

## modified copykat::copykat 
og_copykat <- function(rawmat=rawdata, id.type="S", cell.line="no", ngene.chr=5,LOW.DR=0.05, UP.DR=0.1, win.size=25, norm.cell.names="", KS.cut=0.1, sam.name="", distance="euclidean", output.seg="FALSE", plot.genes="TRUE", genome="hg20", n.cores=1){

  print("running copykat modified -- Using Rcppcluster-- no plotting")
  start_time <- Sys.time()
  set.seed(1234)
  sample.name <- paste(sam.name,"_copykat_", sep="")

  print("running copykat v1.1.0")

  print("step1: read and filter data ...")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", sep=""))

  genes.raw <- apply(rawmat, 2, function(x)(sum(x>0)))
  #genes.raw <-apply_cpp_col(rawmat)

  if(sum(genes.raw> 200)==0) stop("none cells have more than 200 genes")
  if(sum(genes.raw<100)>1){
    rawmat <- rawmat[, -which(genes.raw< 200)]
    #print(paste("filtered out ", sum(genes.raw<=200), " cells with less than 200 genes; remaining ", ncol(rawmat), " cells", sep=""))
  }
  ##
  der<- apply(rawmat,1,function(x)(sum(x>0)))/ncol(rawmat)
  #der=apply_cpp_row_der(rawmat)

  if(sum(der>LOW.DR)>=1){
    rawmat <- rawmat[which(der > LOW.DR), ]; print(paste(nrow(rawmat)," genes past LOW.DR filtering", sep=""))
  }

  WNS1 <- "data quality is ok"
  if(nrow(rawmat) < 7000){
    WNS1 <- "low data quality"
    UP.DR<- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }

  print("step 2: annotations gene coordinates ...")
  if(genome=="hg20"){
  anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
  } else if(genome=="mm10"){
  anno.mat <- annotateGenes.mm10(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
  dim(rawmat)
  }
  anno.mat <- anno.mat[order(as.numeric(anno.mat$abspos), decreasing = FALSE),]

# print(paste(nrow(anno.mat)," genes annotated", sep=""))

  ### module 3 removing genes that are involved in cell cycling

  if(genome=="hg20"){
  HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
  toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), HLAs))
  if(length(toRev)>0){
    anno.mat <- anno.mat[-toRev, ]
  }
  }
#  print(paste(nrow(anno.mat)," genes after rm cell cycle genes", sep=""))
  ### secondary filtering
  ToRemov2 <- NULL
  for(i in 8:ncol(anno.mat)){
    cell <- cbind(anno.mat$chromosome_name, anno.mat[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    } else if(length(rle(cell[,1])$length)<length(unique((anno.mat$chromosome_name)))|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    i<- i+1
  }

  if(length(ToRemov2)==(ncol(anno.mat)-7)) stop("all cells are filtered")
  if(length(ToRemov2)>0){
    anno.mat <-anno.mat[, -which(colnames(anno.mat) %in% ToRemov2)]
  }

  # print(paste("filtered out ", length(ToRemov2), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat<- log(sqrt(rawmat3)+sqrt(rawmat3+1))
  norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x)))
  #norm.mat=apply_cpp_col_norm(norm.mat)
  colnames(norm.mat) <-  colnames(rawmat3)

  #print(paste("A total of ", ncol(norm.mat), " cells, ", nrow(norm.mat), " genes after preprocessing", sep=""))

  ##smooth data
  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c){
    model <- dlm::dlmModPoly(order=1, dV=0.16, dW=0.001)
    x <- dlm::dlmSmooth(norm.mat[, c], model)$s
    x<- x[2:length(x)]
    x <- x-mean(x)
  }

  test.mc <-parallel::mclapply(1:ncol(norm.mat), dlm.sm, mc.cores = n.cores)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE)

  colnames(norm.mat.smooth) <- colnames(norm.mat)

  print("step 4: measuring baselines ...")
  if (cell.line=="yes"){
    print("running pure cell line mode")
        relt <- baseline.synthetic(norm.mat=norm.mat.smooth, min.cells=10, n.cores=n.cores)
        norm.mat.relat <- relt$expr.relat
        CL <- relt$cl
        WNS <- "run with cell line mode"
        preN <- NULL

      } else if(length(norm.cell.names)>1){

        #print(paste(length(norm.cell.names), "normal cells provided", sep=""))
         NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in% norm.cell.names)])
         print(paste(NNN, " known normal cells found in dataset", sep=""))

         if (NNN==0) stop("known normal cells provided; however none existing in testing dataset")
         print("run with known normal...")

         basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% norm.cell.names)],1,median); print("baseline is from known input")
         #basel =apply_cpp_row_basel(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% norm.cell.names)])
          
          #d <- parallelDist::parDist(t(norm.mat.smooth),threads =n.cores, method="euclidean") ##use smooth and segmented data to detect intra-normal cells

          km <- 6
          fit=Rclusterpp.hclust(t(norm.mat.smooth),method="ward",distance="euclidean")
          #fit <- hclust(d, method="ward.D2")
           CL <- cutree(fit, km)

           while(!all(table(CL)>5)){
          km <- km -1
          CL <- cutree(fit, k=km)
         if(km==2){
         break
         }
         }

        WNS <- "run with known normal"
        preN <- norm.cell.names
         ##relative expression using pred.normal cells
        norm.mat.relat <- norm.mat.smooth-basel

        }else {
         basa <- baseline.norm.cl(norm.mat.smooth=norm.mat.smooth, min.cells=5, n.cores=n.cores)
          basel <- basa$basel
          WNS <- basa$WNS
          preN <- basa$preN
          CL <- basa$cl
          if (WNS =="unclassified.prediction"){

                    basa <- baseline.GMM(CNA.mat=norm.mat.smooth, max.normal=5, mu.cut=0.05, Nfraq.cut=0.99,RE.before=basa,n.cores=n.cores)
                    basel <-basa$basel
                    WNS <- basa$WNS

                    preN <- basa$preN

              }
          ##relative expression using pred.normal cells
             norm.mat.relat <- norm.mat.smooth-basel

             }

  ###use a smaller set of genes to perform segmentation
  DR2 <- apply(rawmat3,1,function(x)(sum(x>0)))/ncol(rawmat3)
  #DR2=apply_cpp_row_DR2(rawmat3)

  ##relative expression using pred.normal cells
  norm.mat.relat <- norm.mat.relat[which(DR2>=UP.DR),]

  ###filter cells
  anno.mat2 <- anno.mat[which(DR2>=UP.DR), ]

  ToRemov3 <- NULL
  for(i in 8:ncol(anno.mat2)){
    cell <- cbind(anno.mat2$chromosome_name, anno.mat2[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    } else if(length(rle(cell[,1])$length)<length(unique((anno.mat$chromosome_name)))|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    i<- i+1
  }

  if(length(ToRemov3)==ncol(norm.mat.relat)) stop ("all cells are filtered")

  if(length(ToRemov3)>0){
    norm.mat.relat <-norm.mat.relat[, -which(colnames(norm.mat.relat) %in% ToRemov3)]
   #print(paste("filtered out ", length(ToRemov3), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  }

  #print(paste("final segmentation: ", nrow(norm.mat.relat), " genes; ", ncol(norm.mat.relat), " cells", sep=""))

  CL <- CL[which(names(CL) %in% colnames(norm.mat.relat))]
  CL <- CL[order(match(names(CL), colnames(norm.mat.relat)))]

  print("step 5: segmentation...")
  #results <- rcpp_CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = KS.cut, n.cores=n.cores)
  results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = KS.cut, n.cores=n.cores)

  if(length(results$breaks)<25){
    print("too few breakpoints detected; decreased KS.cut to 50%")
    #results <- rcpp_CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*KS.cut, n.cores=n.cores)
    results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*KS.cut, n.cores=n.cores)
  }

  if(length(results$breaks)<25){
    print("too few breakpoints detected; decreased KS.cut to 75%")
    #results <- rcpp_CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*0.5*KS.cut, n.cores=n.cores)
    results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*0.5*KS.cut, n.cores=n.cores)
  }

  if(length(results$breaks)<25) stop ("too few segments; try to decrease KS.cut; or improve data")

  colnames(results$logCNA) <- colnames(norm.mat.relat)
  results.com <- apply(results$logCNA,2, function(x)(x <- x-mean(x)))

  #results.com=apply_cpp_col_norm(results$logCNA)
  #colnames(results.com)=colnames(results$logCNA)
    

  RNA.copycat <- cbind(anno.mat2[, 1:7], results.com)

  #write.table(RNA.copycat, paste(sample.name, "CNA_raw_results_gene_by_cell.txt", sep=""), sep="\t", row.names = FALSE, quote = F)

  if(genome=="hg20"){
  print("step 6: convert to genomic bins...") ###need multi-core
  #Aj <- rcpp_convert.all.bins.hg20(DNA.mat = DNA.hg20, RNA.mat=RNA.copycat, n.cores = n.cores)
  Aj <- convert.all.bins.hg20(DNA.mat = DNA.hg20, RNA.mat=RNA.copycat, n.cores = n.cores)

  uber.mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])

  print("step 7: adjust baseline ...")

    if(cell.line=="yes"){

               mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
               write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)

                if(distance=="euclidean"){
                 #hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
                 hcc=Rclusterpp.hclust(t(mat.adj),method="ward",distance="euclidean")
 
                  }else {
                 hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
                   }


                  #saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))

                   #plot heatmap
                   #print("step 8: ploting heatmap ...")
                  # my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

                  #  chr <- as.numeric(Aj$DNA.adj$chrom) %% 2+1
                  #  rbPal1 <- colorRampPalette(c('black','grey'))
                  #  CHR <- rbPal1(2)[as.numeric(chr)]
                  #  chr1 <- cbind(CHR,CHR)


                  #  if (ncol(mat.adj)< 3000){
                  #  h <- 10
                  #  } else {
                  #  h <- 15
                  #    }

                  # col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
                  # #library(parallelDist)

                   #if(distance=="euclidean"){
                          #jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
                          #heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
                          #ColSideColors=chr1,Colv=NA, Rowv=TRUE,
                          #notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                          #keysize=1, density.info="none", trace="none",
                          #cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                          #symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
                          #dev.off()
                          ### add a step to plot out gene by cell matrix
                          #if(plot.genes=="TRUE"){

                          #rownames(results.com) <- anno.mat2$hgnc_symbol
                          #chrg <- as.numeric(anno.mat2$chrom) %% 2+1
                          #rbPal1g <- colorRampPalette(c('black','grey'))
                          #CHRg <- rbPal1(2)[as.numeric(chrg)]
                          #chr1g <- cbind(CHRg,CHRg)

                          #pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
                          #heatmap.3(t(results.com),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
                          #ColSideColors=chr1g,Colv=NA, Rowv=TRUE,
                          #notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                          #keysize=1, density.info="none", trace="none",
                          #cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                          #symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
                          #dev.off()
                         #  }
                         #end of ploting gene by cell matrix

               # } else {
                         # jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
                         # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                         # ColSideColors=chr1,Colv=NA, Rowv=TRUE,
                         # notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                         # keysize=1, density.info="none", trace="none",
                         # cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                         # symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
                         # dev.off()
                           ### add a step to plot out gene by cell matrix
            # if(plot.genes=="TRUE"){

                     #     rownames(results.com) <- anno.mat2$hgnc_symbol
                     #     chrg <- as.numeric(anno.mat2$chrom) %% 2+1
                     #     rbPal1g <- colorRampPalette(c('black','grey'))
                     #     CHRg <- rbPal1(2)[as.numeric(chrg)]
                     #     chr1g <- cbind(CHRg,CHRg)

                     #     pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
                     #     heatmap.3(t(results.com),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                     #     ColSideColors=chr1g,Colv=NA, Rowv=TRUE,
                     #     notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                     #     keysize=1, density.info="none", trace="none",
                     #     cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                     #     symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
                     #     dev.off()
                     #     }
                         #end of ploting gene by cell matrix
                      #    }

                          end_time<- Sys.time()
                          print(end_time -start_time)

                         reslts <- list(cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
                         names(reslts) <- c("CNAmat","hclustering")
                         return(reslts)
    } else {
      ########## cell line mode ends here ####################

      #removed baseline adjustment
        if(distance=="euclidean"){
        #hcc <- hclust(parallelDist::parDist(t(uber.mat.adj),threads =n.cores, method = distance), method = "ward.D")
        hcc=Rclusterpp.hclust(t(uber.mat.adj),method="ward",distance="euclidean")

        }else {
        hcc <- hclust(as.dist(1-cor(uber.mat.adj, method = distance)), method = "ward.D")
        }
        hc.umap <- cutree(hcc,2)
        names(hc.umap) <- colnames(results.com)

        cl.ID <- NULL
        for(i in 1:max(hc.umap)){
        cli <- names(hc.umap)[which(hc.umap==i)]
        pid <- length(intersect(cli, preN))/length(cli)
        cl.ID <- c(cl.ID, pid)
        i<- i+1
         }

        com.pred <- names(hc.umap)
        com.pred[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
        com.pred[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
        names(com.pred) <- names(hc.umap)

  ################removed baseline adjustment
        results.com.rat <- uber.mat.adj-apply(uber.mat.adj[,which(com.pred=="diploid")], 1, mean)
        results.com.rat <- apply(results.com.rat,2,function(x)(x <- x-mean(x)))
        results.com.rat.norm <- results.com.rat[,which(com.pred=="diploid")]; dim(results.com.rat.norm)

        cf.h <- apply(results.com.rat.norm, 1, sd)
        base <- apply(results.com.rat.norm, 1, mean)

        adjN <- function(j){
        a <- results.com.rat[, j]
        a[abs(a-base) <= 0.25*cf.h] <- mean(a)
        a
        }


        mc.adjN <-  parallel::mclapply(1:ncol(results.com.rat),adjN, mc.cores = n.cores)
        adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
        colnames(adj.results) <- colnames(results.com.rat)

        #rang <- 0.5*(max(adj.results)-min(adj.results))
        #mat.adj <- adj.results/rang
        mat.adj <- t(t(adj.results)-apply(adj.results,2,mean))

        print("step 8: final prediction ...")

        if(distance=="euclidean"){
         #hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
         hcc=Rclusterpp.hclust(t(mat.adj),method="ward",distance="euclidean")
         }else {
         hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
         }

         hc.umap <- cutree(hcc,2)
         names(hc.umap) <- colnames(results.com)

        #saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))

        cl.ID <- NULL
        for(i in 1:max(hc.umap)){
        cli <- names(hc.umap)[which(hc.umap==i)]
        pid <- length(intersect(cli, preN))/length(cli)
        cl.ID <- c(cl.ID, pid)
        i<- i+1
         }

        com.preN <- names(hc.umap)
        com.preN[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
        com.preN[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
        names(com.preN) <- names(hc.umap)

        if(WNS=="unclassified.prediction"){
        com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
        com.preN[which(com.preN == "aneuploid")] <- "c2:aneuploid:low.conf"
        }

      print("step 9: saving results...")

  ##add back filtered cells as not defined in prediction results
  '%!in%' <- function(x,y)!('%in%'(x,y))

  ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% names(com.preN))]
  if(length(ndef)>0){
    res <- data.frame(cbind(c(names(com.preN),ndef), c(com.preN, rep("not.defined",length(ndef)))))
    colnames(res) <- c("cell.names", "copykat.pred")
  } else {
    res <- data.frame(cbind(names(com.preN), com.preN))
    colnames(res) <- c("cell.names", "copykat.pred")
  }
  ##end
  #write.table(res, paste(sample.name, "prediction.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)

  ####save copycat CNA
  #write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)

  ####%%%%%%%%%%%%%%%%%next heatmaps, subpopulations and tSNE overlay
  # print("step 10: ploting heatmap ...")
  # my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

  # chr <- as.numeric(Aj$DNA.adj$chrom) %% 2+1
  # rbPal1 <- colorRampPalette(c('black','grey'))
  # CHR <- rbPal1(2)[as.numeric(chr)]
  # chr1 <- cbind(CHR,CHR)

  # rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
  # compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]

  # cells <- rbind(compreN_pred,compreN_pred)

  # if (ncol(mat.adj)< 3000){
  #   h <- 10
  # } else {
  #   h <- 15
  # }

  # col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

  # if(distance=="euclidean"){
  # jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
  #  heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
  #           ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
  #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
  #           keysize=1, density.info="none", trace="none",
  #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
  #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

  # legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)
  # dev.off()

  ### add a step to plot out gene by cell matrix
  # if(plot.genes=="TRUE"){
  #   dim(results.com)
  #   rownames(results.com) <- anno.mat2$hgnc_symbol
  #   chrg <- as.numeric(anno.mat2$chrom) %% 2+1
  #   rbPal1g <- colorRampPalette(c('black','grey'))
  #   CHRg <- rbPal1(2)[as.numeric(chrg)]
  #   chr1g <- cbind(CHRg,CHRg)

  #   pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
  #   heatmap.3(t(results.com),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
  #             ColSideColors=chr1g,RowSideColors=cells,Colv=NA, Rowv=TRUE,
  #             notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
  #             keysize=1, density.info="none", trace="none",
  #             cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
  #             symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
  #   dev.off()
  # }
  #end of ploting gene by cell matrix



  # } else {
  #   jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
  #   heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
  #                ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
  #             notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
  #             keysize=1, density.info="none", trace="none",
  #             cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
  #             symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

  #   legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)

  #   dev.off()
  #   ### add a step to plot out gene by cell matrix
  #   if(plot.genes=="TRUE"){
  #     dim(results.com)
  #     rownames(results.com) <- anno.mat2$hgnc_symbol
  #     chrg <- as.numeric(anno.mat2$chrom) %% 2+1
  #     rbPal1g <- colorRampPalette(c('black','grey'))
  #     CHRg <- rbPal1(2)[as.numeric(chrg)]
  #     chr1g <- cbind(CHRg,CHRg)

  #     pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
  #     heatmap.3(t(results.com),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
  #               ColSideColors=chr1g,RowSideColors=cells, Colv=NA, Rowv=TRUE,
  #               notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
  #               keysize=1, density.info="none", trace="none",
  #               cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
  #               symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
  #     dev.off()
  #   }
  #   #end of ploting gene by cell matrix
  # }

#  if(output.seg=="TRUE"){
#   print("generating seg files for IGV viewer")

#   thisRatio <- cbind(Aj$RNA.adj[, 1:3], mat.adj)
#   Short <- NULL
#   chr <- rle(thisRatio$chrom)[[2]]

#   for(c in 4:ncol(thisRatio))
#   {
#     for (x in 1:length(chr)){
#       thisRatio.sub <- thisRatio[which(thisRatio$chrom==chr[x]), ]
#       seg.mean.sub <- rle(thisRatio.sub[,c])[[2]]

#       rle.length.sub <- rle(thisRatio.sub[,c])[[1]]

#       num.mark.sub <- seq(1,length(rle.length.sub),1)
#       loc.start.sub <-seq(1,length(rle.length.sub),1)
#       loc.end.sub <- seq(1,length(rle.length.sub),1)

#       len <-0
#       j <-1

#       for (j in 1: length(rle.length.sub)){
#         num.mark.sub[j] <- rle.length.sub[j]
#         loc.start.sub[j] <- thisRatio.sub$chrompos[len+1]
#         len <- num.mark.sub[j]+len
#         loc.end.sub[j] <- thisRatio.sub$chrompos[len]
#         j <- j+1
#       }

#       ID <- rep(colnames(thisRatio[c]), times=length(rle.length.sub))
#       chrom <- rep(chr[x], times=length(rle.length.sub))
#       Short.sub <- cbind(ID,chrom,loc.start.sub,loc.end.sub,num.mark.sub,seg.mean.sub)
#       Short <- rbind(Short, Short.sub)
#       x <- x+1
#     }
#     c<- c+1
#   }

#   colnames(Short) <- c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
#   head(Short)
#   write.table(Short, paste(sample.name, "CNA_results.seg", sep=""), row.names = FALSE, quote=FALSE, sep="\t")

# }
  end_time<- Sys.time()
  print(end_time -start_time)

  reslts <- list(res, cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
  names(reslts) <- c("prediction", "CNAmat","hclustering")
  return(reslts)
}

  }

  if(genome=="mm10") {
    uber.mat.adj <- data.matrix(results.com)
    dim(uber.mat.adj)
    if(distance=="euclidean"){
      hcc=Rclusterpp.hclust(t(uber.mat.adj),method="ward",distance="euclidean")
      #hcc <- hclust(parallelDist::parDist(t(uber.mat.adj),threads =n.cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(uber.mat.adj, method = distance)), method = "ward.D")
    }
    hc.umap <- cutree(hcc,2)
    names(hc.umap) <- colnames(results.com)

    cl.ID <- NULL
    for(i in 1:max(hc.umap)){
      cli <- names(hc.umap)[which(hc.umap==i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i<- i+1
    }

    com.pred <- names(hc.umap)
    com.pred[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
    com.pred[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
    names(com.pred) <- names(hc.umap)

    ################removed baseline adjustment
    results.com.rat <- uber.mat.adj-apply(uber.mat.adj[,which(com.pred=="diploid")], 1, mean)

    results.com.rat <- apply(results.com.rat,2,function(x)(x <- x-mean(x)))
    results.com.rat.norm <- results.com.rat[,which(com.pred=="diploid")]; dim(results.com.rat.norm)

    cf.h <- apply(results.com.rat.norm, 1, sd)
    base <- apply(results.com.rat.norm, 1, mean)

    adjN <- function(j){
      a <- results.com.rat[, j]
      a[abs(a-base) <= 0.25*cf.h] <- mean(a)
      a
    }


    mc.adjN <-  parallel::mclapply(1:ncol(results.com.rat),adjN, mc.cores = n.cores)
    adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
    colnames(adj.results) <- colnames(results.com.rat)

    #rang <- 0.5*(max(adj.results)-min(adj.results))
    #mat.adj <- adj.results/rang
    mat.adj <- t(t(adj.results)-apply(adj.results,2,mean))

    print("step 8: final prediction ...")

    if(distance=="euclidean"){
      hcc=Rclusterpp.hclust(t(mat.adj),method="ward",distance="euclidean")
      #hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
    }

    hc.umap <- cutree(hcc,2)
    names(hc.umap) <- colnames(results.com)

    saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))

    cl.ID <- NULL
    for(i in 1:max(hc.umap)){
      cli <- names(hc.umap)[which(hc.umap==i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i<- i+1
    }

    com.preN <- names(hc.umap)
    com.preN[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
    com.preN[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
    names(com.preN) <- names(hc.umap)

    if(WNS=="unclassified.prediction"){
      com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
      com.preN[which(com.preN == "aneuploid")] <- "c2:aneuploid:low.conf"
    }

    print("step 9: saving results...")

    ##add back filtered cells as not defined in prediction results
    '%!in%' <- function(x,y)!('%in%'(x,y))
    ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% names(com.preN))]
    if(length(ndef)>0){
      res <- data.frame(cbind(c(names(com.preN),ndef), c(com.preN, rep("not.defined",length(ndef)))))
      colnames(res) <- c("cell.names", "copykat.pred")
    } else {
      res <- data.frame(cbind(names(com.preN), com.preN))
      colnames(res) <- c("cell.names", "copykat.pred")
    }
    ##end
    write.table(res, paste(sample.name, "prediction.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)

    ####save copycat CNA
    write.table(cbind(anno.mat2[, 1:7], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)

    ####%%%%%%%%%%%%%%%%%next heatmaps, subpopulations and tSNE overlay
    # print("step 10: ploting heatmap ...")
    # my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

    # rownames(mat.adj) <- anno.mat2$mgi_symbol
    # chrg <- as.numeric(anno.mat2$chromosome_name) %% 2+1
    # rle(as.numeric(anno.mat2$chromosome_name))
    # rbPal1g <- colorRampPalette(c('black','grey'))
    # CHRg <- rbPal1g(2)[as.numeric(chrg)]
    # chr1g <- cbind(CHRg,CHRg)


    # rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
    # compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]

    # cells <- rbind(compreN_pred,compreN_pred)

    # if (ncol(mat.adj)< 3000){
    #   h <- 10
    # } else {
    #   h <- 15
    # }

    # col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

    # if(distance=="euclidean"){

    #     pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
    #     heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
    #               ColSideColors=chr1g,RowSideColors=cells,Colv=NA, Rowv=TRUE,
    #               notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
    #               keysize=1, density.info="none", trace="none",
    #               cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
    #               symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
    #     dev.off()


    # } else {

    #     pdf(paste(sample.name,"with_genes_heatmap1.pdf",sep=""), height=h*2.5, width=40)
    #     heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
    #               ColSideColors=chr1g,RowSideColors=cells,Colv=NA, Rowv=TRUE,
    #               notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
    #               keysize=1, density.info="none", trace="none",
    #               cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
    #               symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
    #     dev.off()

    #   #end of ploting gene by cell matrix
    # }

    # if(output.seg=="TRUE"){
    #   print("generating seg files for IGV viewer")

    #   thisRatio <- cbind(anno.mat2[, c(2,3,1)], mat.adj)
    #   Short <- NULL
    #   chr <- rle(thisRatio$chromosome_name)[[2]]

    #   for(c in 4:ncol(thisRatio))
    #   {
    #     for (x in 1:length(chr)){
    #       thisRatio.sub <- thisRatio[which(thisRatio$chromosome_name==chr[x]), ]
    #       seg.mean.sub <- rle(thisRatio.sub[,c])[[2]]

    #       rle.length.sub <- rle(thisRatio.sub[,c])[[1]]

    #       num.mark.sub <- seq(1,length(rle.length.sub),1)
    #       loc.start.sub <-seq(1,length(rle.length.sub),1)
    #       loc.end.sub <- seq(1,length(rle.length.sub),1)

    #       len <-0
    #       j <-1

    #       for (j in 1: length(rle.length.sub)){
    #         num.mark.sub[j] <- rle.length.sub[j]
    #         loc.start.sub[j] <- thisRatio.sub$start_position[len+1]
    #         len <- num.mark.sub[j]+len
    #         loc.end.sub[j] <- thisRatio.sub$start_position[len]
    #         j <- j+1
    #       }

    #       ID <- rep(colnames(thisRatio[c]), times=length(rle.length.sub))
    #       chrom <- rep(chr[x], times=length(rle.length.sub))
    #       Short.sub <- cbind(ID,chrom,loc.start.sub,loc.end.sub,num.mark.sub,seg.mean.sub)
    #       Short <- rbind(Short, Short.sub)
    #       x <- x+1
    #     }
    #     c<- c+1
    #   }

    #   colnames(Short) <- c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")

    #   write.table(Short, paste(sample.name, "CNA_results.seg", sep=""), row.names = FALSE, quote=FALSE, sep="\t")

    # }
    end_time<- Sys.time()
    print(end_time -start_time)

    reslts <- list(res, cbind(anno.mat2[, 1:7], mat.adj), hcc)
    names(reslts) <- c("prediction", "CNAmat","hclustering")
    return(reslts)

  }
}
