# Identification of clone-specific and cancer-selective multi-targeting therapeutic options based on single-cell transcriptomic profiles of cancer patients
<br>
  
We introduce NAME, a computational framework for predicting drug combinations based solely on scRNA-seq data.

This tool is based on identification of genetically distinct cancer cell populations(clones), and their transcriptomic differences in individual patient samples, compared to non-malignant healthy cells from the same sample, and then leveraging a reference database of large-scale phenotypic profiles (both transcriptomic and viability readouts) measured in cancer cell lines in response to single-drug perturbations to pre-train a gradient boosting ML model that predicts drug response differences across cell populations.

When applied to patient samples, the model outcome is a list of effective multi-targeting options (either targeted-agents, chemotherapies, or their combinations) that selectively co-inhibit the key cancer subpopulations in the given patient sample, in comparison to non-malignant cell populations. The **concentration-specific drug response predictions** come with confidence quantification to guide the translational applications. 

<br>

## TBD workflow
<p align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/workflow.png">
</p>

Prediction of subclone-specific and cancer-selective compounds is performed in two major steps a-b: 

(a) Raw sequencing data from selected tissue are processed and aligned to give a scRNA-seq expression count matrix which is then processed accordingly.An automated cell type annotation tool, [ScType](https://github.com/IanevskiAleksandr/sc-type) , is applied to accurately identify the cell types in the sample and determine which clusters can be used as reference/healthy cells for the next step. Providing healthy cells at this stage is optional, but it greatly improves the prediction of the next step. These reference cells are used to differentiate between healthy and malignant cells. In this step, three different approaches are employed to build an ensemble prediction that ensures confident calling of healthy cells. 
1. [CopyKAT](https://github.com/navinlabcode/copykat), a Bayesian segmentation method to identify clusters of aneuploid vs diploid cells.
2. [ScType](https://github.com/IanevskiAleksandr/sc-type) + custom marker dataset derived from [CellMarker2.0](http://117.50.127.228/CellMarker/CellMarker_download.html) to develop a marker-based approach to distinguish healthy from malignant clusters. 
3. [SCEVAN](https://github.com/AntonioDeFalco/SCEVAN), employs a segmentation method that utilizes a Mumford and Shah energy model to call cell states:tumor vs normal. 

Finally, a majority vote is taken based on all three predictions to form a confident prediction. This ensemble prediction allows researchers to obtain an accurate overview of the cell types present in the sample and to differentiate between healthy and malignant cells. Thus, it is a valuable tool for various research applications, as shown in the schematic below.

<p align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/ensemble_pred.png">
</p>

Finally, [inferCNV](https://github.com/broadinstitute/infercnv) is applied to infer large-scale copy number variations and estimate genetically distinct subclones from malignant cells. 

(b) Subsequently, subclone-specific differentially-expressed genes are identified through differential expression analysis. These identified genes, along with drug information such as molecular fingerprints and drug doses , are used as inputs for the trained LightGBM model. This model then predicts the most active compounds and their effective doses for each subclone, based on the provided inputs. 

To train the LightGBM model, a comprehensive dataset was compiled with the objective of integrating transcriptional changes observed in small molecule perturbation experiments ([LINCS L1000 dataset](https://clue.io/about)) with drug chemical structures represented as fingerprints and drug-dose response data collected from various sources ([PharmacoDB resource](http://pharmacodb.ca/)). Doses from the LINCS L1000 dataset were matched with dose-response curves obtained from the PharmacoDB resource, and the interpolated cell viability data was used as the outcome variable for prediction model.

<br>


## Prediciting subclone specific drug combinations 
<br>
Although most tools in this analysis require the raw count matrix, it is beneficial to visualize the data at each step of the process. This involves implementing a pre-processing workflow for quality control, normalization, identification of highly variable genes, scaling, and performing principal component analysis (PCA), followed by Uniform Manifold Approximation and Projection (UMAP) or t-Distributed Stochastic Neighbor Embedding (t-SNE). For this we recommend users follow the 
<a href='https://satijalab.org/seurat/articles/pbmc3k_tutorial.html' >Seurat-Guided Clustering Tutorial</a>.

### Step 0: Load the data and the functions
```R
# load functions for identification of healthy and malignant clusters 
source("https://raw.githubusercontent.com/kris-nader/TBD/master/R/identify_mal_norm.R")
# load functions for identification of genetically distinct subclones
source("https://raw.githubusercontent.com/kris-nader/TBD/master/R/identify_subclones.R")
# load functions for predicting subclone specific therapeutic options

patient_sample=readRDS("./example_data.RDS")
```
### Step 1: Automated Cell type annotation with ScType
In this step, we utilize a standard ScType workflow, which requires the single cell RNAseq object and the tissue type of interest as input. By default, the function employs the predefined scType database, containing markers for various tissues such as the immune system, pancreas, liver, eye, kidney, brain, lung, adrenal gland, heart, intestine, muscle, placenta, spleen, stomach, and thymus. We refer users to the <a href="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx">ScTypeDB</a> for more information on the defined cell markers.

However, users can easily customize the analysis by uploading their own custom marker database for their specific tissue of interest using the `custom_marker_db` parameter. In short, the custom marker database should resemble that of the ScTypeDB(xlsx format) with four columns (tissue type, cell name, geneSymbolmore 1-- positive markers , and geneSymbolmore2--negative markers). In this tutorial, the sample was derived from a patient with Acute Myeloid Leukemia (AML), and we will identify cell types using `tissue=Immune System` parameter. The resulting cell types can be visualized on the UMAP using `Seurat::DimPlot`.


```R
sctype_source()
patient_sample=run_sctype(patient_sample,tissue="Immune System",plot=TRUE)
```
### Step 2: Identification of malignant/normal clusters
In this step, we use multiple tools to generate a confident ensemble prediction. To improve the accuracy of the predictions, we recommend using the normal cells identified in step 1 as input for copyKat and SCEVAN. Afterwards, an ensemble prediction is constructed based on the combined results of these tools, which takes advantage of the distinct approaches to confidently identify both healthy and malignant cell clusters. The function `runEnsemble` will execute each of these tools(copyKat,scType+new markers,SCEVAN) and then compute the ensemble prediciton. We can also visualize the results of each individual tool and of the ensemble prediction.
```R
normal_cells=c("NKT-like cells","CD4+ T cells", "Naive B cells")
patient_sample=runEnsemble(patient_sample,known_tissue_type="AML",known_normal_cells=normal_cells)
visualize_ensembl_step(patient_sample)
```

### Step 3: Identification of genetically distinct subclones
This step uses healthy/reference cells identified by step 2(ensemble model) to identify genetically distinct sublcones. Note that this step may be computationally intensive. 
```R
patient_sample=runinferCNV(patient_sample,"healthy")
```


### Step 4: Comparative analysis of the subclone and normal cluster to extract subclone specific DEG
We will focus on more broad levels subclones in this tutorrial, but more specific subclones can be used in this step. For example, for subclones A and B:
```R
subcloneA=subclone_DEG(patient_sample,"A","healthy")
subcloneB=subclone_DEG(patient_sample,"B","healthy")
```


### Step 5: Use subclone specific DEG as input to the pre-trained LightGBM model.
For each run of `run_drug_combo_pred`, the model predicts drug:dose:%inhibtition based on a predefined set of drug:dose:response integrated from LINCS L1000 and PharmacoDB. To predict response for a drug not included in the database, refer to our next section on predicting response of new drugs.
```R
subcloneA_drugs=run_drug_combo_pred(subcloneA)
subcloneB_drugs=run_drug_combo_pred(subcloneB)
```

## Predicting response of new drug:dose
<br>


## Contact information
For any questions please contact **Aleksandr Ianevski** [aleksandr.ianevski@helsinki.fi] and  **Kristen Nader** [kristen.nader@helsinki.fi]

## Reference Papers
1. Ianevski, A., Giri, A. K. &amp; Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Nature Communications 13, (2022). 
2. Gao, R. et al. Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes. Nature Biotechnology 39, 599–608 (2021). 
3. De Falco, A., Caruso, F., Su, X.-D., Iavarone, A. &amp; Ceccarelli, M. A variational algorithm to detect the clonal copy number substructure of tumors from scrna-seq data. Nature Communications 14, (2023). 
4. Hu, C. et al. CellMarker 2.0: An updated database of manually curated cell markers in human/mouse and web tools based on scRNA-Seq Data. Nucleic Acids Research 51, (2022). 
5. 
