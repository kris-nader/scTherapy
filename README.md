# Identification of clone-specific and cancer-selective multi-targeting therapeutic options based on single-cell transcriptomic profiles of cancer patients
<br>
**Article**: yabba dabba do

<p style="text-align:justify;"> <b>name of the tool</b> description</p>

We introduce NAME to adress the limitations of various machine learning (ML) methods have been developed to predict effective anticancer drug combinations using large-scale multi-omics data. 

This tool is based on identification of genetically distinct cancer cell populations(clones), and their transcriptomic differences in individual patient samples, compared to non-malignant healthy cells from the same sample, and then leveraging a reference database of large-scale phenotypic profiles (both transcriptomic and viability readouts) measured in cancer cell lines in response to single-drug perturbations to pre-train a gradient boosting ML model that predicts drug response differences across cell populations.

When applied to patient samples, the model outcome is a list of effective multi-targeting options (either targeted-agents, chemotherapies, or their combinations) that selectively co-inhibit the key cancer subpopulations in the given patient sample, in comparison to non-malignant cell populations. The **concentration-specific drug response predictions** come with confidence quantification to guide the translational applications. 

<br>

## TBD workflow
<p align="center"> 
<img src="https://github.com/kris-nader/TBD/blob/main/workflow.png">
</p>

Prediction of subclone-specific and cancer-selective compounds is performed in two major steps a-b: 

(a) Raw sequencing data from selected tissue are processed and aligned to give a scRNA-seq expression count matrix. Unsupervised clustering of cells is then performed to separate malignant and healthy cell clusters using two analytical tools, [ScType](https://github.com/IanevskiAleksandr/sc-type) and [CopyKAT](https://github.com/navinlabcode/copykat). Finally, [inferCNV](https://github.com/broadinstitute/infercnv) is applied to infer large-scale copy number variations and estimate genetically distinct subclones from malignant cells. 

(b) Subsequently, subclone-specific differentially-expressed genes are identified through differential expression analysis. These identified genes, along with drug information such as molecular fingerprints and drug doses , are used as inputs for the pre-trained LightGBM model. This model then predicts the most active compounds and their effective doses for each subclone, based on the provided inputs. 

To pre-train the LightGBM model, a comprehensive dataset was compiled with the objective of integrating transcriptional changes observed in small molecule perturbation experiments ([LINCS L1000 dataset](https://clue.io/about)) with drug chemical structures represented as fingerprints and drug-dose response data collected from various sources ([PharmacoDB resource](http://pharmacodb.ca/)). Doses from the LINCS L1000 dataset were matched with dose-response curves obtained from the PharmacoDB resource, and the interpolated cell viability data was used as the outcome variable for prediction model.

<br>

### Step 1: Identification of malignant/normal clusters
First you will need to identify clusters of normal and malignant cells. We propose an ensembl of 3 tools: ScType,SCEVAN, and scATOMIC which will be merged with the input seurat object.

```R
patient_sample=runEnsembl(patient_sample)
```
If everything is successful, you should observe an output analogous to the following:
```
###################################################
## Running ensembl tool: ScType,SCEVAN, scATOMIC ##
###################################################

Success: ScType time=0.07 
Success: SCEVAN time=0.07 
Success: scATOMIC time=0.07 
...
Ensembl model: time=0.07 

Done!
```

### Step 2: Identification of genetically distinct subclones
This step is computationally intensive. 

```R
patient_sample=runinferCNV(patient_sample)
```
Once again, If everything is successful, you should observe an output analougous to the following: 
If everything is successful, you should observe an output analogous to the following:
```
#####################################
## Running subclone tool: inferCNV ##
#####################################

Success: inferCNV time=0.07 

Done!
```

```R
R code here

```

## Contact information
For any questions please contact **Aleksandr Ianevski** [aleksandr.ianevski@helsinki.fi] and  **Kristen Nader** [kristen.nader@helsinki.fi]

