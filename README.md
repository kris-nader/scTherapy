# Identification of clone-specific and cancer-selective multi-targeting therapeutic options based on single-cell transcriptomic profiles of cancer patients
<br>
  
We introduce NAME, a computational framework that uses only single-cell transcriptomic data to predict personalized monotherapies and multi-targeting drug combinations for cancer patients.

The tool consists of two main steps:
1. Our tool utilizes a semi-automated approach to categorize cells into normal and malignant, followed by the separation of malignant cells into subpopulations/clones.
2. Next, it identifies differentially expressed genes between normal and malignant cells, which are subsequently used as input to our pre-trained machine learning model. The model predicts potential treatment options that effectively target malignant cells while minimizing toxicity to normal cells.

For more information, please refer to original publication [to be filled].
<br><br>
<b><h2>TBD workflow</h2></b>
<span align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/full_workflow_ensemble.png">
</span>

Prediction of subclone-specific and cancer-selective compounds is performed in two major steps using only the expression count matrix.

1. An automated cell type annotation tool, [ScType](https://github.com/IanevskiAleksandr/sc-type) is used  to accurately identify the cell types and  determine which clusters can be used as reference/healthy cells for the next step. Then an ensemble prediction is done using three different approaches to ensure confident calling of healthy cells([CopyKAT](https://github.com/navinlabcode/copykat),[ScType](https://github.com/IanevskiAleksandr/sc-type) + [CellMarker2.0](http://117.50.127.228/CellMarker/CellMarker_download.html) and [SCEVAN](https://github.com/AntonioDeFalco/SCEVAN)) using reference cells identified through cell type annotation. Finally, [inferCNV](https://github.com/broadinstitute/infercnv) is applied to identify genetically distinct subclones from malignant cells. 

2. Subclone-specific differentially-expressed genes are identified through differential expression analysis. These identified genes, along with drug information such as molecular fingerprints and drug doses , are used as inputs for the trained LightGBM model. This model then predicts the most active compounds and their effective doses for each subclone, based on the provided inputs.

We created a comprehensive dataset that combines transcriptional changes observed in small molecule perturbation experiments ([LINCS L1000 dataset](https://clue.io/about)), drug chemical structures represented as fingerprints, and drug-dose response data ([PharmacoDB resource](http://pharmacodb.ca/)). By matching doses from the LINCS L1000 dataset with dose-response curves from PharmacoDB, we obtained interpolated cell viability data as the outcome variable for our prediction model and trained a LightGBM model.


<br>


## Prediciting subclone specific drug combinations 
<br>

### Step -1: Load required functions
```R

lapply(c("dplyr","Seurat","HGNChelper","openxlsx","copykat","ggplot2","SCEVAN", "cowplot","Rcpp","Rclusterpp",
          "parallel","biomaRt","infercnv","logger","httr", "jsonlite", "readr"), library, character.only = !0)
         
source("https://raw.githubusercontent.com/kris-nader/R/identify_mal_norm.R"); 
source("https://raw.githubusercontent.com/kris-nader/R/identify_subclones.R"); 
source("https://raw.githubusercontent.com/kris-nader/R/predict_compounds.R"); 

```
### Step 0: Process the data
Although most tools in this analysis require the raw count matrix, it is beneficial to visualize the data at each step of the process.

First, let's load an example PBMC scRNA-seq dataset, consisting of ~3000 cells obtained from a human AML patient. Next, we normalize and cluster our raw count matrix using Seurat (see <a href="https://satijalab.org/seurat/articles/pbmc3k_tutorial.html">Seurat tutorial for more details</a>). The raw data can be found <a href='https://raw.githubusercontent.com/kris-nader/TBD/main/sample_x_exp.rawdata.txt.zip'>here</a>.


```R
# Load example dataset or upload your own expression matrix (rows - genes, column - cells) as e.g.: data = read.table("exp.rawdata.txt", header = TRUE, row.names = 1, sep = "\t").
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

# check number of PC 
ElbowPlot(patient_sample)

# clustering 
patient_sample = FindNeighbors(patient_sample, dims = 1:10)
patient_sample = FindClusters(patient_sample, resolution = 0.8)
patient_sample = RunUMAP(patient_sample, dims = 1:10)

```
### Step 1: Automated Cell type annotation with ScType
In this step, we use a standard ScType workflow, which requires the single cell RNAseq object and the tissue type of interest as input. We refer users to the <a href="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx">ScTypeDB</a> for more information on the defined cell markers and to the original <a href="https://www.nature.com/articles/s41467-022-28803-w">publication</a>.

In this tutorial, the sample was derived from a patient with Acute Myeloid Leukemia (AML), and we will identify cell types using `known_tissue_type=Immune system` parameter. The resulting cell types can be visualized on the UMAP using `Seurat::DimPlot`.

```R
patient_sample=run_sctype(patient_sample,known_tissue_type="Immune system",plot=FALSE)
```
### Step 2: Identification of malignant/normal clusters
In this step, we use multiple tools to generate a confident ensemble prediction. To improve the accuracy of the predictions, we recommend using the normal cells identified in step 1 as input for copyKat and SCEVAN. The `runEnsemble` function executes each tool (copyKat, scType+new markers, SCEVAN) and computes the ensemble prediction. 

```R
norm_cells=get_normal_cells(patient_sample,c("Memory CD4+ T cells","CD8+ NKT-like cells"))
patient_sample=run_ensemble(patient_sample,disease="AML",known_normal_cells=norm_cells,plot=FALSE)
visualize_ensemble_step(patient_sample)
```

<p align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/example_ensemble.png">
</p>

### Step 3: Identification of genetically distinct subclones
This step uses healthy/reference cells identified by step 2(ensemble model) to identify genetically distinct sublcones. Note that this step may be computationally intensive. 
```R
patient_sample=run_infercnv(patient_sample)
```
<p align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/example_infercnv.png">
</p>

### Step 4: Extract subclone specific DEG
We will focus on broad levels subclones in this tutorial, but more specific subclones can be used in this step for more specific analysis. For subclones A and B:
```R
subcloneA=subclone_DEG(patient_sample,"A","healthy")
subcloneB=subclone_DEG(patient_sample,"B","healthy")
```

### Step 5: Use subclone specific DEG as input to the pre-trained LightGBM model.
First, we begin by filtering the differentially expressed genes based on 2 critera:
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

