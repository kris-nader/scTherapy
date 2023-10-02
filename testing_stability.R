## testing stability of detecting normal and malignant cells at different %

# Load required libraries and source functions
invisible(lapply(c("dplyr","Seurat","HGNChelper","openxlsx","copykat","copykatRcpp","ggplot2", "SCEVAN","yaGST","cowplot","Rcpp","Rclusterpp",
                   "parallel","biomaRt","logger","httr", "jsonlite", "readr","future"), library, character.only = !0))

invisible(lapply(c("https://raw.githubusercontent.com/kris-nader/TBD/main/R/identify_mal_norm.R",
                   "https://raw.githubusercontent.com/kris-nader/TBD/main/R/identify_subclones.R",
                   "https://raw.githubusercontent.com/kris-nader/TBD/main/R/predict_compounds.R"),source))


## randomly remove subset of cells
percentages_to_remove = c(25, 50, 75)

data = read_zip("https://raw.githubusercontent.com/kris-nader/scTherapy/main/sample_x_exp.rawdata.txt.zip")
# Loop through each percentage and remove data
for (percentage in percentages_to_remove) {
    # Calculate the number of cells to keep
    cols_to_keep = sample(1:ncol(data), size = floor((100 - percentage) / 100 * ncol(data)))
    # Create a new data frame with the selected cells
    new_data = data[,cols_to_keep]
    # Print information about the removed percentage
    cat(sprintf("Removed %d%% of data. Remaining cells: %d\n", percentage, nrow(new_data)))
    # Update the data frame for the next iteration
    patient_sample = CreateSeuratObject(counts = new_data)	
    # simple filtering
    patient_sample[["percent.mt"]] = PercentageFeatureSet(patient_sample, pattern = "^MT-")
    patient_sample = subset(patient_sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    # normalize data
    patient_sample = NormalizeData(patient_sample, normalization.method = "LogNormalize", scale.factor = 10000)
    patient_sample = FindVariableFeatures(patient_sample, selection.method = "vst", nfeatures = 2000)
    # scale and PCA
    patient_sample = ScaleData(patient_sample, features = rownames(patient_sample))
    patient_sample = RunPCA(patient_sample, features = VariableFeatures(object = patient_sample))
    # check number of PC (optional)
    ElbowPlot(patient_sample)
    # Cluster and visualize
    patient_sample = FindNeighbors(patient_sample, dims = 1:10)
    patient_sample = FindClusters(patient_sample, resolution = 0.8)
    patient_sample = RunUMAP(patient_sample, dims = 1:10)
 
    patient_sample=run_sctype(patient_sample,known_tissue_type="Immune system",plot=TRUE)
    # please note that this step is time consuming (~10 minutes for example data), consider running on faster multi-core Linux or MacOS-based PC to speed up this process
    norm_cells=get_normal_cells(patient_sample,c("Memory CD4+ T cells"))
    patient_sample=run_ensemble(patient_sample,disease="AML",known_normal_cells=norm_cells,plot=FALSE)
    
    #write to a file to compare later
    write.table(patient_sample@meta.data,file=paste0("robust_analys_",percentage,"_metadata.txt"))
    visualize_ensemble_step(patient_sample)

}
