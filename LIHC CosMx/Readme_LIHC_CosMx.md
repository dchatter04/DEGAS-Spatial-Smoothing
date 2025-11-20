# Analysis for LIHC samples from CosMx

# Step 1 Bulk RNA-seq data cleanup
# Step 2 Preprocessing
# Step 3 DEGAS model training and predictions
# Step 4 Postprocessing smoothing
# Step 5 Other relevant plots


## **Step 1 Bulk RNA-seq data cleanup**

* Download the bulkRNA-seq gene expression and clinical data from https://gdac.broadinstitute.org/ your local Drive.
  - Source of the gene-exp data is https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LIHC/20160128/gdac.broadinstitute.org_LIHC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz.
  - Clinical covariates are in https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LIHC/20160128/gdac.broadinstitute.org_LIHC.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz
* Run the lines in script **BulkRNA-seq_datatcleanup.R** 
* Preprocessed patient-level gene expression data will be stored as **patDat.csv** and the labels of each patient
  whether they are classified as Type 2 Diabetic (t2d) or non-diabetic (no-diasese) is provided in **patLab.csv** file.

## **Step 2 Preprocessing of the CosMx data (single cell spatially resolved transcriptomics (scSRT))**

* Download the Type 2 Diabetes Xenium data for 2 tissue samples (2 type 2 diabetic, and 2 non-diabetic) from <a href="https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/">[https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/](https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/)</a> into your local drive
* Run the lines in script **Preprocess_CosMx.R**
* Preprocessed combined gene expression data is stored as **stDat.rds** file.

## **Step 3 DEGAS model training and predictions**

* Run the code **Model_training_and_prediction.R**
* Align the bulkRNA-seq and scSRT data with common set of highly variable genes and train the model. Save the trained model in **model1.RDS** file.
* Obtain the saved patient -level predictions in **preds1.csv** file.

## **Step 4 Postprocessing smoothing**

* Run the code **SwiS_smoothing.R**
* Sliding Window (**swSmoothAtlas**) function is provided in the list of functions document.
* Smoothed predictions for each sample is individually stored in **swSmoothed_ND1.rds**, **swSmoothed_T2D1.rds**, **swSmoothed_ND2.rds**, and **swSmoothed_T2D2.rds** respectively.

## **Step 5 Other relevant plots**

* Run the script **Plots.R** for the relavant boxplots and the pairwise Wilcoxon test results for DEGAS prediction scores for non-diabetic vs type 2 diabetic samples.
