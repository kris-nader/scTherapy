# Identification of clone-specific and cancer-selective multi-targeting therapeutic options based on single-cell transcriptomic profiles of cancer patients
<br>

**Article**: [to be filled]
  
<p style="text-align:justify;"> <b>XXX</b> is a computational framework to predict personalized monotherapies and multi-targeting drug combinations for cancer patients using only single-cell transcriptomic data.</p>

##
<br>
<p align="center">
  <img src="https://github.com/kris-nader/TBD/blob/main/ensemble_pred_workflow.png" width="80%">
</p>
<br>

The tool consists of two main steps:
1. Our tool utilizes a semi-automated approach to categorize cells into healthy and malignant, followed by the separation of malignant cells into subpopulations/clones <i>(Fig 1, panel a)</i>.
2. Next, it identifies differentially expressed genes between healthy and malignant cells, which are subsequently used as input to our pre-trained machine learning model. The model predicts potential treatment options that effectively target malignant cells while minimizing toxicity to healthy cells <i>(Fig 1, panel b)</i>.

For more information, please refer to original publication [to be filled].
<br><br>


##

### Step 1.0: Load libraries and process the data

<p>First, let's load all the necessary libraries and source functions required for the analysis.</p>


<details>
  <summary>Install libraries by clicking <b>HERE</b></summary>
	
	
  ```R
# run this code to install required libraries
  packages <- c("dplyr","Seurat","HGNChelper","openxlsx","copykat","copykatRcpp","ggplot2","SCEVAN", "cowplot",
			  "Rcpp","Rclusterpp","parallel","biomaRt","logger","httr", "jsonlite", "readr")
	      
install_load_packages <- function(packages){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    if (!requireNamespace("devtools", quietly = TRUE))
        install.packages("devtools")
    
    sapply(packages, function(pkg){
        if (!require(pkg, character.only = TRUE)){
            if (pkg %in% c("copykat", "yaGST", "SCEVAN", "Rclusterpp")) {
                tryCatch({
                    if (pkg == "copykat") {
                        devtools::install_github("navinlabcode/copykat")
                    } else if (pkg == "yaGST") {
                        devtools::install_github("miccec/yaGST")
                    } else if (pkg == "SCEVAN") {
                        devtools::install_github("AntonioDeFalco/SCEVAN")
                    } else if (pkg == "Rclusterpp") {
                        devtools::install_github("nolanlab/Rclusterpp")
                    }else if (pkg == "copykatRcpp") {
                        devtools::install_github('IanevskiAleksandr/copykatRcpp')
                    }
                    library(pkg, character.only = TRUE)
                }, error = function(e){
                    install_from_CRAN_or_Bioconductor(pkg)
                })
            } else {
                install_from_CRAN_or_Bioconductor(pkg)
            }
        }
    })
}

install_from_CRAN_or_Bioconductor <- function(pkg) {
  tryCatch({
    install.packages(pkg); library(pkg, character.only = TRUE)
  }, error = function(e){
    BiocManager::install(pkg); library(pkg, character.only = TRUE)
  })
}

install_load_packages(packages)

  ```
				     
</details>



```R
# Load required libraries and source functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx","copykat","copykatRcpp","ggplot2", "yaGST", "SCEVAN", "cowplot","Rcpp","Rclusterpp",
          "parallel","biomaRt","logger","httr", "jsonlite", "readr"), library, character.only = !0)
         
source("https://raw.githubusercontent.com/kris-nader/TBD/main/R/identify_mal_norm.R"); 
source("https://raw.githubusercontent.com/kris-nader/TBD/main/R/identify_subclones.R"); 
source("https://raw.githubusercontent.com/kris-nader/TBD/main/R/predict_compounds.R"); 
```


<br>	
<p>Next, let's load an example PBMC scRNA-seq dataset, consisting of ~3000 cells obtained from a human AML patient. The count matrix is normalized and clustered using Seurat (see <a href="https://satijalab.org/seurat/articles/pbmc3k_tutorial.html" target="_blank">Seurat tutorial for more details</a>). An example raw count matrix data can be found <a href='https://raw.githubusercontent.com/kris-nader/TBD/main/sample_x_exp.rawdata.txt.zip' target="_blank">here</a>.</p>
	

 ```R
# Load example dataset 
# You can also upload your own expression matrix: data = read.table("your.exp.rawdata.txt", header = TRUE, row.names = 1, sep = "\t").
data <- read_zip("https://raw.githubusercontent.com/kris-nader/TBD/main/sample_x_exp.rawdata.txt.zip")
patient_sample = CreateSeuratObject(counts = data)	
	
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
DimPlot(patient_sample, reduction = "umap")
```
<br>

### Step 1.1: Automated Cell type annotation with ScType
<p>In this step, we use our method for fully-automated cell type annotation called ScType. It only requires an scRNA-seq object and the name of the tissue type as input, please see ScType GitHub for more details: <a href="https://github.com/IanevskiAleksandr/sc-type">ScType</a>. </p>
<p>For our AML patient sample, we specify <code>known_tissue_type</code> as <code>Immune system</code>, but other possible tissue types include: Immune system, Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus.</p>

```R
patient_sample=run_sctype(patient_sample,known_tissue_type="Immune system",plot=TRUE)
```

<br>

### Step 1.2: Identification of malignant/healthy clusters
<p>In this step, we use an ensemble of three tools <i>(CopyKat, scType+new markers, and SCEVAN)</i> to confidently classify cells into malignant and healthy groups. To enhance accuracy, we recommend providing prior knowledge of a confident healthy cell cluster as input to <code>runEnsemble</code> function.</p>
<p>Here, we provide <code>T cells</code> as confident healthy cell cluster, given that <code>T cells</code> are considered as normal cells in AML - <a href="https://doi.org/10.1016/j.cell.2019.01.031">PMID:30827681</a></p>

```R
# please note that this step is time consuming (~10 minutes for example data), consider running on faster multi-core Linux or MacOS-based PC to speed up this process
norm_cells=get_normal_cells(patient_sample,c("Memory CD4+ T cells"))
patient_sample=run_ensemble(patient_sample,disease="AML",known_normal_cells=norm_cells,plot=FALSE)
visualize_ensemble_step(patient_sample)
```

<p align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/example_ensemble.png">
</p>
	
To predict single agent therapies. we refer users to the [Tutorial on predicting monotherpies using malignant specific DEG](#tutorial-on-predicting-monotherpies-using-malignant-specific-DEG)

To predict combination therapies, we refer users to the [Tutorial on predicting combination therapies using subclone specific DEG](#tutorial-on-predicting-combination-therapies-using-subclone-specific-deg). This includes infering clonal architecture within the malignant cluster identified in step 1.2.

	
	
## Tutorial on predicting monotherpies using malignant specific DEG
### Step 1.3: Extract malignant cluster specific DEG
At this point, users can use TBD to predict monotherapies to target the malignant cluster identified in step 1.2
```R
malignant_cells_DEG=subclone_DEG(patient_sample,"malignant","healthy",monotherapies=TRUE)
```
### Step 1.4: Use malignant specific DEG as input to the pre-trained LightGBM model.
First, we begin by filtering the differentially expressed genes based on 2 critera with <code>process_DEG</code>:
1. (p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) 
2. (avg_log2FC > -0.1 & avg_log2FC < 0.1)

Then, we can run <code>predict_drugs</code> using the malignant DEGs as input to predict monotherpies.
					 
```R
# filter DEGS
DEG_malignant=process_DEG(malignant_cells_DEG)
#predict monotherpies for malignant cluster
monotherapy_drugs=predict_drugs(DEG_malignant)
```
	
## Tutorial on predicting combination therapies using subclone specific DEG
If investigating tumor subclones is of interest, users can continue from step 1.2 with the following steps
	
### Step 1.3: Identification of genetically distinct subclones
This step uses healthy/reference cells identified by step 1.2(ensemble model) to identify genetically distinct sublcones. Note that this step may be computationally intensive. 
```R
patient_sample=run_infercnv(patient_sample)
```
<p align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/example_infercnv.png">
</p>

### Step 1.4: Extract subclone specific DEG
We will focus on broad levels subclones in this tutorial, but more specific subclones can be used in this step for more specific analysis. For subclones A and B:
```R
subcloneA=subclone_DEG(patient_sample,"A","healthy")
subcloneB=subclone_DEG(patient_sample,"B","healthy")
```

### Step 1.5: Use subclone specific DEG as input to the pre-trained LightGBM model.
First, we begin by filtering the differentially expressed genes based on 2 critera with :
1. (p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) 
2. (avg_log2FC > -0.1 & avg_log2FC < 0.1)

```R
DEG_A=process_DEG(subcloneA)
DEG_B=process_DEG(subcloneB)
```
For each run of `predict_drugs`, the model predicts drug:dose based on a predefined set of drug:dose:response integrated from LINCS L1000 and PharmacoDB. 

```R
subcloneA_drugs=predict_drugs(DEG_A)
subcloneB_drugs=predict_drugs(DEG_B)
```


 

## Contact information
For any questions please contact **Aleksandr Ianevski** [aleksandr.ianevski@helsinki.fi] and  **Kristen Nader** [kristen.nader@helsinki.fi]

## Copyright and license

Code copyright 2023 TBD, [https://github.com/kris/sc-type/blob/master/LICENSE](https://github.com/kris-nader/TBD/blob/main/LICENSE)

## Reference Papers
1. Ianevski, A., Giri, A. K. &amp; Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Nature Communications 13, (2022). 
2. Gao, R. et al. Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes. Nature Biotechnology 39, 599–608 (2021). 
3. De Falco, A., Caruso, F., Su, X.-D., Iavarone, A. &amp; Ceccarelli, M. A variational algorithm to detect the clonal copy number substructure of tumors from scrna-seq data. Nature Communications 14, (2023). 
4. Hu, C. et al. CellMarker 2.0: An updated database of manually curated cell markers in human/mouse and web tools based on scRNA-Seq Data. Nucleic Acids Research 51, (2022). 
5. Subramanian, A. et al. A Next Generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. Cell 171, (2017). 
6. Smirnov, P. et al. PharmacoDB: An integrative database for mining in vitro anticancer drug screening studies. Nucleic Acids Research 46, (2017). 

