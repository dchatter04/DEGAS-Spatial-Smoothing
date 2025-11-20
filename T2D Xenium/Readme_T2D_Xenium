# **Analysis of T2D Xenium data**

## **Step 1 Bulk RNA-seq data cleanup**

* Download the bulkRNA_seq_Diabetes.csv file (<a href="https://github.com/dchatter04/DEGAS-Spatial-Smoothing/blob/main/T2D%20Xenium/bulkRNA_seq_Diabetes.csv">[https://github.com/dchatter04/DEGAS-Spatial-Smoothing/blob/main/T2D%20Xenium/bulkRNA_seq_Diabetes.csv](https://github.com/dchatter04/DEGAS-Spatial-Smoothing/blob/main/T2D%20Xenium/bulkRNA_seq_Diabetes.csv)</a>)
  into your local Drive.
* Run the lines in script BulkRNA-seq_datatcleanup.R (<a href="https://github.com/dchatter04/DEGAS-Spatial-Smoothing/blob/main/T2D%20Xenium/BulkRNA-seq_datacleanup.R">[https://github.com/dchatter04/DEGAS-Spatial-Smoothing/blob/main/T2D%20Xenium/BulkRNA-seq_datacleanup.R](https://github.com/dchatter04/DEGAS-Spatial-Smoothing/blob/main/T2D%20Xenium/BulkRNA-seq_datacleanup.R)</a>)
* Preprocessed patient-level gene expression data will be stored as **patDat.csv** and the labels of each patient
  whether they are classified as Type 2 Diabetic (t2d) or non-diabetic (no-diasese) is provided in **patLab.csv** file.

## **Step 2 Preprocessing of the Xenium data (single cell spatially resolved transcriptomics (scSRT))**

* Download the Type 2 Diabetes Xenium data for 4 tissue samples (2 type 2 diabetic, and 2 non-diabetic) from <a href="https://www.synapse.org/Synapse:syn68699752/files/">[https://www.synapse.org/Synapse:syn68699752/files/](https://www.synapse.org/Synapse:syn68699752/files/)</a> into your local drive
* Run the lines in script **Preprocess_Xenium.R**
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
