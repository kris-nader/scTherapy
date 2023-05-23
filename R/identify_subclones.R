###########################################################################
## Get functions for checkpoint 2: identification of malignant subclones ##
############################################################################
#
# GNU General Public License v3.0 (https://github.com/kris-nader/TBD/blob/main/LICENSE)

# Written by Kristen Michelle Nader <kristen.nader@helsinki.fi>
# and Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi> , April 2021
#
# Functions on this page:
# get_gene_annotation,get_normal_cells,run_infercnv
#


lapply(c("dplyr","Seurat","HGNChelper","biomaRt","infercnv"), library, character.only = T)


#' @title Get gene annotation information from biomaRt
#' @name get_gene_annotation
#' @description Query biomaRt to retrieve chromosome name, HGNC symbol, start and end position of genes in a Seurat object. 
#' @details Used with infercnv run
#' @param seurat_object A Seurat object with gene names as row names.
#' @return  A data frame with gene name, chromosome, start and end positions.
#' 
#' @import biomaRt
#' @import dplyr
#' 
#' @examples
#' genes_annotation=get_gene_annotation(seurat_object)
#' 
#' @export
#' 
get_gene_annotation <- function(seurat_object){
    requireNamespace("biomaRt")
    ensembl = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    genes = biomaRt::getBM(
        attributes = c("hgnc_symbol","chromosome_name","start_position","end_position"),
        filters = "hgnc_symbol", 
        values = rownames(seurat_object), 
        mart = ensembl)
    
    colnames(genes) = c("hgnc_symbol","chr","start","stop")
    #remove duplicated genes (keep those with only one position)
    genes = genes[!(duplicated(genes$hgnc_symbol) | duplicated(genes$hgnc_symbol, fromLast = !0)), ]
    rownames(genes) = genes$hgnc_symbol
    # filter & merge
    genes = genes[genes$chr %in% c(1:22,"X", "Y"),]
    inters_ = intersect(rownames(genes), rownames(seurat_object))
    seurat_object = seurat_object[inters_, ]
    genes = genes[inters_, ]
    genes$chr=paste0("chr",genes$chr)
    genes$hgnc_symbol = NULL
    colnames(genes) <- NULL
    return(genes)
}

#' @title Get normal cells
#' @name get_normal_cells
#' @description Identify normal cells in a Seurat object and return a dataframe with cell names and their classification
#' @details Used with infercnv run
#' @param seurat_object A Seurat object
#' @param known_normal_cells A vector of cell names known to be normal (optional)
#' @return A dataframe with two columns: "cell" (cell names) and "class" (cell classification)
#' 
#' @import Seurat
#' 
#' @examples
#' annotation_file_for_infer <- get_normal_cells(seurat_object, known_normal_cells)
#' @export

get_normal_cells <- function(seurat_object, known_normal_cells = NULL) {
    # Get cell classifications from Seurat object
    annotation_file_for_infer = seurat_object@meta.data[,"ensemble_output", drop=!1]
    # If no known normal cells were provided, use the information from step 1: ensemble_output
    if (is.null(known_normal_cells)) {
        known_normal_cells=rownames(seurat_object@meta.data[which(seurat_object@meta.data$ensemble_output=="healthy"),])
    }
    annotation_file_for_infer[known_normal_cells,"ensemble_output"] = "healthy"
    colnames(annotation_file_for_infer)=NULL
    # Return dataframe with cell names and classifications
    return(annotation_file_for_infer)
}



#' @title  Run infercnv analysis
#' @name run_infercnv
#' @description Runs infercnv to identify tumor subclones
#' @details 
#' @param seurat_object a seurat object containing raw count matrix and metadata 
#' @param known_normal_cells optional vector of known normal cell names default: will use those identified in step 1 
#' @return a Seurat object with infercnv subclones defined in metadata
#' 
#' @import Seurat
#' @import infercnv
#' @import biomaRt
#' 
#' @export
#' @examples
#' seurat_object <- run_infercnv(seurat_object, known_normal_cells = NULL, ncores = 4)
#' 
#' @export


run_infercnv <- function(seurat_object, known_normal_cells = NULL, ncores = 4) {
    # Get gene annotation
    genes_annot <- get_gene_annotation(seurat_object)
    # Get normal cells for annotation
    annotation_file_for_infer <- get_normal_cells(seurat_object, known_normal_cells)
    
    # Create InferCNV object
    infercnv_obj <- infercnv::CreateInfercnvObject(
        raw_counts_matrix = as.matrix(seurat_object@assays$RNA@counts), 
        gene_order_file = genes_annot, 
        annotations_file = annotation_file_for_infer, 
        ref_group_names = "healthy"
    )
    
    # Run InferCNV
    infercnv_obj = infercnv::run(infercnv_obj, cutoff=0.1,
                                 out_dir="./", cluster_by_groups=!1, 
                                 plot_steps=FALSE, scale_data=!0,
                                 denoise=!0, noise_filter=0.12, 
                                 analysis_mode='subclusters',
                                 tumor_subcluster_partition_method = "random_trees",
                                 HMM=!0, HMM_type='i6',
                                 num_threads=ncores, 
                                 up_to_step = 17, save_rds = !1, 
                                 save_final_rds = !0, no_plot = TRUE, 
                                 no_prelim_plot = TRUE,output_format=NA,inspect_subclusters=FALSE)
    ## read file
    
    grouping_structure=read.table("./17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings",header=TRUE,row.names = "cell")
    #rownames(grouping_structure)=grouping_structure$cell
    
    
    # Find subclones
    subclones=rownames(grouping_structure %>% filter(startsWith(cell_group_name,"all_observations.all_observations")))
    healthy=rownames(grouping_structure %>% filter(startsWith(cell_group_name,"healthy.healthy")))
    # Update cell group names for subclones
    if (length(subclones) > 0) {
        grouping_structure[subclones,"cell_group_name"]=substr(grouping_structure[subclones,"cell_group_name"], 35, 41)
    }
    if (length(healthy) > 0) {
        grouping_structure[healthy,"cell_group_name"]="healthy"
    }
    seurat_object1=seurat_object
    seurat_object1@meta.data[rownames(grouping_structure),"infercnvgroupings"]=grouping_structure[rownames(grouping_structure),"cell_group_name"]
    
    # Add subclones to Seurat object metadata
    seurat_object1@meta.data[which(startsWith(seurat_object1@meta.data$infercnvgroupings,"1.1.")),"infercnv_broad_groupings"]="A"
    seurat_object1@meta.data[which(startsWith(seurat_object1@meta.data$infercnvgroupings,"1.2.")),"infercnv_broad_groupings"]="B"
    seurat_object1@meta.data[which(startsWith(seurat_object1@meta.data$infercnvgroupings,"healthy")),"infercnv_broad_groupings"]="healthy"

    return(seurat_object1)
}


#' @title  Run DEG analysis
#' @name subclone_DEG
#' @description Identifies differentially expressed clone specific genes
#' @details 
#' @param seurat_object a seurat object after identification of subclones( see run_infercnv for more information)
#' @param subclone: subclone of interest
#' @param known_normal_cells: known and annotated normal cell names default: will use those identified in step 1 
#' @return clone specific DEGs
#' 
#' @import Seurat
#' 
#' @export
#' @examples
#' seurat_object <- subclone_DEG(seurat_object, subclone="A",known_normal_cells ="healthy")
#' 
#' @export

subclone_DEG <- function(seurat_object, subclone = NULL, known_normal_cells="healthy",save=FALSE){
    temp=seurat_object
    Idents(temp)=temp@meta.data$infercnv_broad_groupings
    degRes_subcome = FindMarkers(object = temp, ident.1 = subclone, ident.2 = known_normal_cells, min.pct = -Inf, logfc.threshold = -Inf,
                                 min.cells.feature = 1, min.cells.group = 1,test.use="wilcox")
    if(save==TRUE){
        file_name=paste0("./DEG_",subclone,"_healthy.txt")
        write.table(degRes_subcome,file=file_name,quote=FALSE)
    }
    return(degRes_subcome)                           
}



