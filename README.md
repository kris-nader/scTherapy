# Identification of clone-specific and cancer-selective multi-targeting therapeutic options based on single-cell transcriptomic profiles of cancer patients (scTherapy)
<br>

**Article**: [to be filled]
  
<p style="text-align:justify;"> <b>scTherapy</b> is a computational framework to predict personalized monotherapies and multi-targeting drug combinations for cancer patients using only single-cell transcriptomic data.</p>

##
<br>
<p align="center">
  <img src="https://github.com/kris-nader/TBD/blob/main/ensemble_pred_workflow.png" width="80%">
</p>
<br>

The tool consists of two main steps:
1. A semi-automated approach to categorize cells into healthy and malignant, followed by the separation of malignant cells into subpopulations/clones <i>(Fig 1, panel a)</i>.
2. Next, it identifies differentially expressed genes between healthy and malignant cells, which are subsequently used as input to our pre-trained machine learning model. The model predicts potential treatment options that effectively target malignant cells while minimizing toxicity to healthy cells <i>(Fig 1, panel b)</i>.

For more information, please refer to original publication [to be filled].
<br><br>


## Quick Start 
<p>Prepare a <b>'malignant_cells_DEG'</b> data.frame showing differential gene expression results between cancer and normal cells. A data.frame row names should correspond to <b>gene symbols</b> and data.frame itself should feature 2 columns:</p>
<ul><li>log2 fold-change ('<b>avg_log2FC</b>') indicating up-regulation <i>(positive values)</i> or down-regulation <i>(neg. values)</i> in cancer cells</li> <li>adjusted p-values ('<b>p_val_adj</b>')</li></ul>

  ```R
# load example data
malignant_cells_DEG = readRDS(url("https://raw.githubusercontent.com/kris-nader/scTherapy/main/malignant_cells_DEG_example.RDS"))

# > head(malignant_cells_DEG)
#           avg_log2FC    p_val_adj
# LYZ       5.254787      1.512117e-107
# S100A9    4.530434      9.073954e-39
# 	         .............
# IL32 	   -3.599663      3.197208e-164
# GNLY     -4.417322      8.259763e-81
  ```

<p>Predict multi-targeting monotherapies with scTherapy:</p>

  ```R
 # load libraries
invisible(lapply(c("dplyr","logger","httr", "jsonlite", "readr"), library, character.only = !0));
invisible(source("https://raw.githubusercontent.com/kris-nader/scTherapy/main/R/predict_compounds.R"))
 
# load data for making drug:dose predictions
gene_list <- "https://raw.githubusercontent.com/kris-nader/scTherapy/main/geneinfo_beta_input.txt"
gene_info <- data.table::fread(gene_list) %>% as.data.frame()

# filter DEG
DEG_malignant <- malignant_cells_DEG %>%
    mutate(gene_symbol = rownames(.)) %>% inner_join(gene_info, by = "gene_symbol") %>%
    filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
DEG_malignant_list <- setNames(as.list(DEG_malignant$avg_log2FC), DEG_malignant$gene_symbol)	

#predict monotherapies for malignant cluster
monotherapy_drugs <- predict_drugs(DEG_malignant_list)

  ```


## Predicting therapies using scRNAseq

### Optional: Set up docker
This step is <b>optional</b> and used to replicate the exact package versioning (environment) utilized in the publication. The Docker image is pre-configured with the following libraries: Seurat4, inferCNV, copyKat, SCEVAN, and biomaRt. Utilizing Docker eliminates the need for manual package installations, allowing you to proceed directly to step 1 for loading packages and processing data, by running R from docker. 

```posh
# Pull the scTherapy Docker image from Docker Hub
docker pull kmnader/sctherapy
# Run docker (Ensure that the directory /absolute/path/sctherapy_dir exists on your machine)
docker run -it --name sctherapy_docker -v /absolute/path/sctherapy_dir:/home kmnader/sctherapy:latest bash  
```


### Step 1.0: Load libraries and process the data

<p>First, let's load all the necessary libraries and source functions required for the analysis.</p>


<details>
  <summary>Install libraries by clicking <b>HERE</b></summary>
	
	
  ```R
# run this code to install required libraries
  packages = c("dplyr","Seurat","HGNChelper","openxlsx","copykat","copykatRcpp","ggplot2","SCEVAN","yaGST","cowplot",
              "Rcpp","Rclusterpp","parallel","logger","httr", "jsonlite", "readr")

install_load_packages <- function(packages){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        library("BiocManager")
    
    if (!requireNamespace("devtools", quietly = TRUE))
        install.packages("devtools")
        library("devtools")
    
    sapply(packages, function(pkg){
        if (!require(pkg, character.only = TRUE)){
            if (pkg %in% c("copykat", "yaGST", "SCEVAN", "Rclusterpp","copykatRcpp")) {
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
invisible(lapply(c("dplyr","Seurat","HGNChelper","openxlsx","copykat","copykatRcpp","ggplot2", "SCEVAN","yaGST","cowplot","Rcpp","Rclusterpp",
          "parallel","biomaRt","logger","httr", "jsonlite", "readr","future"), library, character.only = !0))
         
invisible(lapply(c("https://raw.githubusercontent.com/kris-nader/scTherapy/main/R/identify_mal_norm.R",
		   "https://raw.githubusercontent.com/kris-nader/scTherapy/main/R/identify_subclones.R",
		   "https://raw.githubusercontent.com/kris-nader/scTherapy/main/R/predict_compounds.R"),source))
```


<br>	
<p>Next, let's load an example PBMC scRNA-seq dataset, consisting of ~3000 cells obtained from a human AML patient. The count matrix is normalized and clustered using Seurat (see <a href="https://satijalab.org/seurat/articles/pbmc3k_tutorial.html" target="_blank">Seurat tutorial for more details</a>). An example raw count matrix data can be found <a href='https://raw.githubusercontent.com/kris-nader/scTherapy/main/sample_x_exp.rawdata.txt.zip' target="_blank">here</a>.</p>
	

 ```R
# Load example dataset 
# You can also upload your own expression matrix: data = read.table("your.exp.rawdata.txt", header = TRUE, row.names = 1, sep = "\t").
data <- read_zip("https://raw.githubusercontent.com/kris-nader/scTherapy/main/sample_x_exp.rawdata.txt.zip")
patient_sample <- CreateSeuratObject(counts = data)	
	
# simple filtering
patient_sample[["percent.mt"]] <- PercentageFeatureSet(patient_sample, pattern = "^MT-")
patient_sample <- subset(patient_sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalize data
patient_sample <- NormalizeData(patient_sample, normalization.method = "LogNormalize", scale.factor = 10000)
patient_sample <- FindVariableFeatures(patient_sample, selection.method = "vst", nfeatures = 2000)

# scale and PCA
patient_sample <- ScaleData(patient_sample, features = rownames(patient_sample))
patient_sample <- RunPCA(patient_sample, features = VariableFeatures(object = patient_sample))

# check number of PC (optional)
ElbowPlot(patient_sample)

# Cluster and visualize
patient_sample <- FindNeighbors(patient_sample, dims = 1:10)
patient_sample <- FindClusters(patient_sample, resolution = 0.8)
patient_sample <- RunUMAP(patient_sample, dims = 1:10)
DimPlot(patient_sample, reduction = "umap")
```
<br>

### Step 1.1: Automated Cell type annotation with ScType
<p>In this step, we use our method for fully-automated cell type annotation called ScType. It only requires an scRNA-seq object and the name of the tissue type as input, please see ScType GitHub for more details: <a href="https://github.com/IanevskiAleksandr/sc-type">ScType</a>. </p>
<p>For our AML patient sample, we specify <code>known_tissue_type</code> as <code>Immune system</code>, but other possible tissue types include: Immune system, Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus.</p>

```R	
patient_sample <- run_sctype(patient_sample,known_tissue_type="Immune system",plot=TRUE)
```
<br>

### Step 1.2: Identification of malignant/healthy clusters
<p>In this step, we use an ensemble of three tools <i>(CopyKat, scType+new markers, and SCEVAN)</i> to confidently classify cells into malignant and healthy groups. To enhance accuracy, we recommend providing prior knowledge of a confident healthy cell cluster as input to <code>runEnsemble</code> function.</p>
<p>Here, we provide <code>T cells</code> as confident healthy cell cluster, given that <code>T cells</code> are considered as normal cells in AML - <a href="https://doi.org/10.1016/j.cell.2019.01.031">PMID:30827681</a></p>

```R
# please note that this step is time consuming (~10 minutes for example data), consider running on faster multi-core Linux or MacOS-based PC to speed up this process
norm_cells <- get_normal_cells(patient_sample,c("Memory CD4+ T cells"))
patient_sample <- run_ensemble(patient_sample,disease="AML",known_normal_cells=norm_cells,plot=FALSE)
visualize_ensemble_step(patient_sample)
```

<p align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/ensembl_.png">
</p>


Users have the option to choose between two options: making predictions to target the entire malignant tumor mass or exploring clonal lineages to predict combination therapies aimed at targeting specific subclones.
1. A tutorial on prediciting single agent therapies: [Predict monotherapies using malignant specific DEG](#predict-monotherapies-using-malignant-specific-deg)
2. A tutorial on prediciting combination therapies: [Predict combination therapies using subclone specific DEG](#predict-combination-therapies-using-subclone-specific-deg). This includes infering clonal architecture within the malignant cluster identified in step 1.2.

	
	
## Predict monotherapies using malignant specific DEG
### Step 2.1: Extract malignant cluster specific DEG
At this point, users can use TBD to predict monotherapies to target the malignant cluster identified in step 1.2
```R
plan("multisession", workers = 4)
malignant_cells_DEG <- clone_DEG(patient_sample,malignant_identifier="malignant",known_normal_cells="healthy",save=FALSE)
```
<br>
	
### Step 2.2: Use malignant specific DEG as input to the pre-trained LightGBM model.
First, we begin by filtering the differentially expressed genes based on 2 critera:
1. (p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) 
2. (avg_log2FC > -0.1 & avg_log2FC < 0.1)

Then, we can run <code>predict_drugs</code> using the malignant DEGs as input to predict monotherpies.
					 
```R
# filter DEGS
gene_list <- "https://raw.githubusercontent.com/kris-nader/scTherapy/main/geneinfo_beta_input.txt"
gene_info <- data.table::fread(gene_list) %>% as.data.frame()

DEG_malignant <- malignant_cells_DEG %>%
    mutate(gene_symbol = rownames(.)) %>% inner_join(gene_info, by = "gene_symbol") %>%
    filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
DEG_malignant_list <- setNames(as.list(DEG_malignant$avg_log2FC), DEG_malignant$gene_symbol)
	
#predict monotherapies for malignant cluster
monotherapy_drugs <- predict_drugs(DEG_malignant_list)
```
	
## Predict combination therapies using subclone specific DEG
If users are interested in exploring tumor subclones and targeting them with specific compounds, they can follow the subsequent steps outlined below.
	
### Step 3.1: Identification of genetically distinct subclones
This step uses healthy/reference cells identified by step 1.2(ensemble model) to identify genetically distinct sublcones. Note that this step may be computationally intensive. 
	
<details>
  <summary>Help installing infercnv by clicking <b>HERE</b></summary>
	
  ```R
  packages=c("biomaRt","rjags","infercnv")
  install_load_packages(packages)
  ```
				     
</details>
	
```R
# Load extra required libraries 
lapply(c("rjags","biomaRt","infercnv"), library, character.only = !0)
# runs infercnv and all helper functions for the analysis 
patient_sample <- run_infercnv(patient_sample)
```
<p align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/example_infercnv.png">
</p>

<br>

### Step 3.2: Extract subclone specific DEG
We will focus on broad levels subclones in this tutorial, but more specific subclones can be used in this step for more specific analysis. For subclones A and B:
```R
plan("multisession", workers = 4)
subcloneA <- clone_DEG(patient_sample,malignant_identifier="A",known_normal_cells="healthy",save=FALSE)
subcloneB <- clone_DEG(patient_sample,malignant_identifier="B",known_normal_cells="healthy",save=FALSE)
```
<br>

### Step 3.3: Use subclone specific DEG as input to the pre-trained LightGBM model.
First, we begin by filtering the differentially expressed genes based on 2 critera with :
1. (p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) 
2. (avg_log2FC > -0.1 & avg_log2FC < 0.1)

```R
# filter subclone "A" specific DEG					 
DEG_subclone_A <- subcloneA %>%
    mutate(gene_symbol = rownames(.)) %>% inner_join(gene_info, by = "gene_symbol") %>%
    filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
DEG_subclone_A_list <- setNames(as.list(DEG_subclone_A$avg_log2FC), DEG_subclone_A$gene_symbol)      
	
# filter subclone "B" specific DEG	
DEG_subclone_B <- subcloneB %>%
    mutate(gene_symbol = rownames(.)) %>% inner_join(gene_info, by = "gene_symbol") %>%
    filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
DEG_subclone_B_list <- setNames(as.list(DEG_subclone_B$avg_log2FC), DEG_subclone_B$gene_symbol)   
```
For each run of `predict_drugs`, the model predicts drug:dose based on a predefined set of drug:dose:response integrated from LINCS L1000 and PharmacoDB. 

```R
# predict subclone A specific drug:dose
subcloneA_drugs <- predict_drugs(degs_list=DEG_subclone_A_list)

# predict subclone B specific drug:dose
subcloneB_drugs <- predict_drugs(degs_list=DEG_subclone_B_list)
```

## Notes for reproducibility 
```R
sessionInfo();
R version 4.3.0 (2023-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] future_1.32.0      readr_2.1.4        jsonlite_1.8.4     httr_1.4.5        
 [5] logger_0.2.2       biomaRt_2.56.1     Rclusterpp_0.2.6   Rcpp_1.0.10       
 [9] cowplot_1.1.1      yaGST_2017.08.25   doParallel_1.0.17  iterators_1.0.14  
[13] foreach_1.5.2      SCEVAN_1.0.1       ggplot2_3.4.2      copykatRcpp_0.1.0 
[17] copykat_1.1.0      openxlsx_4.2.5.2   HGNChelper_0.8.1   SeuratObject_4.1.3
[21] Seurat_4.3.0       dplyr_1.1.2       

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.20        splines_4.3.0           later_1.3.0            
  [4] bitops_1.0-7            filelock_1.0.3          tibble_3.2.1           
  [7] polyclip_1.10-4         XML_3.99-0.16           lifecycle_1.0.3        
 [10] globals_0.16.2          lattice_0.21-8          MASS_7.3-59            
 [13] magrittr_2.0.3          plotly_4.10.1           httpuv_1.6.9           
 [16] sctransform_0.3.5       zip_2.3.0               sp_1.6-0               
 [19] spatstat.sparse_3.0-1   reticulate_1.28         pbapply_1.7-0          
 [22] DBI_1.1.3               RColorBrewer_1.1-3      abind_1.4-5            
 [25] zlibbioc_1.46.0         Rtsne_0.16              purrr_1.0.1            
 [28] BiocGenerics_0.46.0     RCurl_1.98-1.12         rappdirs_0.3.3         
 [31] GenomeInfoDbData_1.2.10 IRanges_2.34.0          S4Vectors_0.38.1       
 [34] ggrepel_0.9.3           irlba_2.3.5.1           listenv_0.9.0          
 [37] spatstat.utils_3.0-3    goftest_1.2-3           spatstat.random_3.1-5  
 [40] fitdistrplus_1.1-11     parallelly_1.35.0       leiden_0.4.3           
 [43] codetools_0.2-19        xml2_1.3.4              tidyselect_1.2.0       
 [46] matrixStats_1.2.0       stats4_4.3.0            BiocFileCache_2.8.0    
 [49] spatstat.explore_3.2-1  ellipsis_0.3.2          progressr_0.13.0       
 [52] ggridges_0.5.4          survival_3.5-5          tools_4.3.0            
 [55] progress_1.2.3          ica_1.0-3               glue_1.6.2             
 [58] gridExtra_2.3           GenomeInfoDb_1.36.0     withr_2.5.0            
 [61] fastmap_1.1.1           fansi_1.0.4             digest_0.6.31          
 [64] R6_2.5.1                mime_0.12               colorspace_2.1-0       
 [67] scattermore_1.1         tensor_1.5              spatstat.data_3.0-1    
 [70] RSQLite_2.3.4           utf8_1.2.3              tidyr_1.3.0            
 [73] generics_0.1.3          data.table_1.14.8       prettyunits_1.1.1      
 [76] htmlwidgets_1.6.2       uwot_0.1.14             pkgconfig_2.0.3        
 [79] gtable_0.3.3            blob_1.2.4              lmtest_0.9-40          
 [82] XVector_0.40.0          htmltools_0.5.5         scales_1.2.1           
 [85] Biobase_2.60.0          png_0.1-8               tzdb_0.4.0             
 [88] reshape2_1.4.4          nlme_3.1-162            curl_5.0.0             
 [91] zoo_1.8-12              cachem_1.0.7            stringr_1.5.0          
 [94] KernSmooth_2.23-20      miniUI_0.1.1.1          AnnotationDbi_1.62.2   
 [97] pillar_1.9.0            grid_4.3.0              vctrs_0.6.5            
[100] RANN_2.6.1              promises_1.2.0.1        dbplyr_2.3.4           
[103] xtable_1.8-4            cluster_2.1.4           cli_3.6.1              
[106] compiler_4.3.0          rlang_1.1.1             crayon_1.5.2           
[109] future.apply_1.11.0     plyr_1.8.8              stringi_1.7.12         
[112] viridisLite_0.4.2       deldir_1.0-9            munsell_0.5.0          
[115] Biostrings_2.68.1       lazyeval_0.2.2          spatstat.geom_3.2-1    
[118] Matrix_1.6-1            hms_1.1.3               patchwork_1.1.2        
[121] bit64_4.0.5             KEGGREST_1.40.1         shiny_1.7.4            
[124] ROCR_1.0-11             igraph_1.4.2            memoise_2.0.1          
[127] bit_4.0.5 
```
## Contact information
For any questions please contact **Aleksandr Ianevski** [aleksandr.ianevski@helsinki.fi] and  **Kristen Nader** [kristen.nader@helsinki.fi]

## Copyright and license

Code copyright 2023 scTherapy, [https://github.com/kris-nader/scTherapy/blob/main/LICENSE](https://github.com/kris-nader/scTherapy/blob/main/LICENSE)

## Reference Papers
1. Ianevski, A., Giri, A. K. &amp; Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Nature Communications 13, (2022). 
2. Gao, R. et al. Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes. Nature Biotechnology 39, 599–608 (2021). 
3. De Falco, A., Caruso, F., Su, X.-D., Iavarone, A. &amp; Ceccarelli, M. A variational algorithm to detect the clonal copy number substructure of tumors from scrna-seq data. Nature Communications 14, (2023). 
4. Hu, C. et al. CellMarker 2.0: An updated database of manually curated cell markers in human/mouse and web tools based on scRNA-Seq Data. Nucleic Acids Research 51, (2022). 
5. Subramanian, A. et al. A Next Generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. Cell 171, (2017). 
6. Smirnov, P. et al. PharmacoDB: An integrative database for mining in vitro anticancer drug screening studies. Nucleic Acids Research 46, (2017). 

