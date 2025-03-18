This GitHub repository contains the data files and analysis code used in **"Detecting homologous recombination deficiency for breast cancer through integrative analysis of genomic data"** published in _..._: 

Link: https://www...

The files are organised into three folders:

* _R_: The R code to reproduce all analyses.

_parameter_function_: The parameters and custom functions.

_training_model_: Develop the feature selection, Leave-One-Out Cross-Validation (LOOCV) training, and threshold selection process, and assess the performance of LOOCV.

_LOOCV_: The Leave-One-Out Cross-Validation (LOOCV) process. 

_ssc_model_: Construct the final five models using a random forest-based self-training method applied to two cohorts, incorporating specifically chosen features and thresholds.

_validation_: Validation on other cohorts, such as TransNEO, NEWTON and CCLE datasets; and show performance. 


* _data_: The features, models and list of genes required to perform the analyses described in the paper. 

_model.self.rf.RData_: the final five random forest-based self-training models;

_loc.feature0.name_: the list of the names of all features in each model;

_loc.feature.name_: the list of the names of selected features in each model;

_gene781_: a curated list of 781 genes related to cancer and the homologous recombination pathway.


* _output_: The predicted HRD probabilities in five models for all cohorts.


* _validation_SCANB.R_: is an example of using SCAN-B data to illustrate how to use our model for prediction.


```
#--load R packages, parameters, features' names and our five models.
set.seed(1) 
library(ssc) 
library(randomForest) 
source('R/parameter_function.R')  
load('data/loc.feature.name.RData') 
load('data/model.self.rf.RData') 

#--load  and SCANB dataset, i.e. Data.SCANB includes the features of 231 samples, where the first column is the sample name, the second column is the HRD label with factor levels = c(1,0), and the 3:53 columns are the considered 51 features.
load('data.SCANB.Rdata') #231 53

#--data conversion by scale(log(x+1)): 
matlog.SCANB=data.SCANB
matlog.SCANB[,c(3:53)]=scale(log(data.SCANB[,c(3:53)]+1))

#--establish the dataset for 5 models:
mat.model.SCANB=list(
model1=data.frame(matlog.SCANB[,match(c('id','is.hrd',loc.feature.name$model1),colnames(matlog.SCANB))]),
model2=data.frame(matlog.SCANB[,match(c('id','is.hrd',loc.feature.name$model2),colnames(matlog.SCANB))]),
model3=data.frame(matlog.SCANB[,match(c('id','is.hrd',loc.feature.name$model3),colnames(matlog.SCANB))]),
model4=data.frame(matlog.SCANB[,match(c('id','is.hrd',loc.feature.name$model4),colnames(matlog.SCANB))]),
model5=data.frame(matlog.SCANB[,match(c('id','is.hrd',loc.feature.name$model5),colnames(matlog.SCANB))])
                     )

#--sequentially input the 5 models, and obtain the prediction probabilities and classifications for each model, where Prob.SCANB [m] is the prediction probability of the m-th model.
pred.SCANB=prob.SCANB=list() 
for(m in 1:5){   
  prob.SCANB[[m]]=ssc:::predProb(model.self.rf[[m]]$model,data.matrix(mat.model.SCANB[[m]][,-c(1,2)]),model.self.rf[[m]]$pred,model.self.rf[[m]]$pred.pars)[,1]
  pred.SCANB[[m]]=ifelse(prob.SCANB[[m]]>=mat.cut.f1$cut[m],1,0) 
}
