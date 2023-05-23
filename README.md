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
<img src="https://github.com/kris-nader/TBD/blob/main/workflow.png">
</span>
Prediction of subclone-specific and cancer-selective compounds is performed in two major steps using only the expression count matrix.

(a) An automated cell type annotation tool, [ScType](https://github.com/IanevskiAleksandr/sc-type) is used  to accurately identify the cell types and  determine which clusters can be used as reference/healthy cells for the next step. Providing healthy cells at this stage is optional, but it greatly improves the prediction of the next step. Then an ensemble prediction is done using three different approaches to ensure confident calling of healthy cells.

1. [CopyKAT](https://github.com/navinlabcode/copykat), a Bayesian segmentation method to identify clusters of aneuploid vs diploid cells.
2. [ScType](https://github.com/IanevskiAleksandr/sc-type) + custom marker dataset derived from [CellMarker2.0](http://117.50.127.228/CellMarker/CellMarker_download.html) to develop a marker-based approach to distinguish healthy from malignant clusters. 
3. [SCEVAN](https://github.com/AntonioDeFalco/SCEVAN), employs a segmentation method that utilizes a Mumford and Shah energy model to call cell states:tumor vs normal. 

<p align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/ensemble_pred.png">
</p>

Finally, [inferCNV](https://github.com/broadinstitute/infercnv) is applied to identify genetically distinct subclones from malignant cells. 

(b) Subclone-specific differentially-expressed genes are identified through differential expression analysis. These identified genes, along with drug information such as molecular fingerprints and drug doses , are used as inputs for the trained LightGBM model. This model then predicts the most active compounds and their effective doses for each subclone, based on the provided inputs.

To train the LightGBM model, we created a comprehensive dataset that combines transcriptional changes observed in small molecule perturbation experiments ([LINCS L1000 dataset](https://clue.io/about)), drug chemical structures represented as fingerprints, and drug-dose response data ([PharmacoDB resource](http://pharmacodb.ca/)). By matching doses from the LINCS L1000 dataset with dose-response curves from PharmacoDB, we obtained interpolated cell viability data as the outcome variable for our prediction model.


<br>


## Prediciting subclone specific drug combinations 
<br>
Although most tools in this analysis require the raw count matrix, it is beneficial to visualize the data at each step of the process. This involves implementing a pre-processing workflow for quality control, normalization, identification of highly variable genes, scaling, and performing principal component analysis (PCA), followed by Uniform Manifold Approximation and Projection (UMAP) or t-Distributed Stochastic Neighbor Embedding (t-SNE). For this we recommend users follow the 
<a href='https://satijalab.org/seurat/articles/pbmc3k_tutorial.html' >Seurat-Guided Clustering Tutorial</a>.

### Step -1: Load the data and the functions

### Step 0: Process the data
First, let's load an example PBMC scRNA-seq dataset, consisting of ~3000 cells obtained from a human AML patient. Next, we normalize and cluster our raw count matrix using Seurat (see Seurat tutorial for more details, <a href="https://satijalab.org/seurat/articles/pbmc3k_tutorial.html">https://satijalab.org/seurat/articles/pbmc3k_tutorial.html</a>). The raw data can be found <a href='https://raw.githubusercontent.com/kris-nader/TBD/main/sample_x_exp.rawdata.txt.zip'>here</a>.

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

#visualize
DimPlot(patient_sample, reduction = "umap")
```
### Step 1: Automated Cell type annotation with ScType
In this step, we utilize a standard ScType workflow, which requires the single cell RNAseq object and the tissue type of interest as input. By default, the function employs the predefined scType database, containing markers for various tissues such as the immune system, pancreas, liver, eye, kidney, brain, lung, adrenal gland, heart, intestine, muscle, placenta, spleen, stomach, and thymus. We refer users to the <a href="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx">ScTypeDB</a> for more information on the defined cell markers.

Users can easily customize the analysis by uploading their own custom marker database for their specific tissue of interest using the `custom_marker_db` parameter. In short, the custom marker database should resemble that of the ScTypeDB(xlsx format) with four columns (tissue type, cell name, geneSymbolmore 1-- positive markers , and geneSymbolmore2--negative markers). In this tutorial, the sample was derived from a patient with Acute Myeloid Leukemia (AML), and we will identify cell types using `known_tissue_type=Immune system` parameter. The resulting cell types can be visualized on the UMAP using `Seurat::DimPlot`.

```R
patient_sample=run_sctype(patient_sample,known_tissue_type="Immune system",plot=FALSE)
```
### Step 2: Identification of malignant/normal clusters
In this step, we use multiple tools to generate a confident ensemble prediction. To improve the accuracy of the predictions, we recommend using the normal cells identified in step 1 as input for copyKat and SCEVAN. Afterwards, an ensemble prediction is constructed based on the combined results of these tools, which takes advantage of the distinct approaches to confidently identify both healthy and malignant cell clusters. The function `runEnsemble` will execute each of these tools(copyKat,scType+new markers,SCEVAN) and then compute the ensemble prediciton. We can also visualize the results of each individual tool and of the ensemble prediction.
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

### Step 4: Comparative analysis of the subclone and normal cluster to extract subclone specific DEG
We will focus on broad levels subclones in this tutorial, but more specific subclones can be used in this step for more specific analysis. For subclones A and B:
```R
subcloneA=subclone_DEG(patient_sample,"A","healthy")
subcloneB=subclone_DEG(patient_sample,"B","healthy")
```

### Step 5: Use subclone specific DEG as input to the pre-trained LightGBM model.
First, we begin by filtering the differentially expressed genes based on 2 critera:

The first condition (p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) filters genes with adjusted p-values less than or equal to 0.05 and an average log2 fold change greater than 1 or less than -1. This condition selects genes that are significantly differentially expressed.

The second condition (avg_log2FC > -0.1 & avg_log2FC < 0.1) filters genes with an average log2 fold change greater than -0.1 and less than 0.1. This condition selects genes that are not significantly differentially expressed but have a small fold change.
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

