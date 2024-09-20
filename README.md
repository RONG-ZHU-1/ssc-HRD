This GitHub repository contains the data files and analysis code used in **"Detecting homologous recombination deficiency for breast cancer through integrative analysis of genomic data"** published in _..._: 

Link: https://www...

The files are organised into three folders:

* _R_: The R code to reproduce all analyses.

_parameter_function_: The parameters and custom functions.

_training_model_: Develop the feature selection, Leave-One-Out Cross-Validation (LOOCV) training, and threshold selection process, and assess the performance of LOOCV.

_LOOCV_: The Leave-One-Out Cross-Validation (LOOCV) process. 

_ssc_model_: Construct the final five models using a random forest-based self-training method applied to three cohorts, incorporating specifically chosen features and thresholds.

_validation_: Validation on other cohorts, such as TransNEO, NEWTON and CCLE datasets; and show performance. 


* _data_: The features, models and list of genes required to perform the analyses described in the paper. 

_model.self.rf.RData_: the final five random forest-based self-training models;

_loc.feature0.name_: the list of the names of all features in each model;

_loc.feature.name_: the list of the names of selected features in each model;

_gene781_: a curated list of 781 genes related to cancer and the homologous recombination pathway.


* _validation_SCANB.R_: is an example of using SCAN-B data to illustrate how to use our model for prediction.

* _output_: The predicted HRD probabilities in five models for all cohorts.

