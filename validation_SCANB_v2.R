#--load R packages, parameters, features' name and our five models.
set.seed(1) 
library(ssc) 
library(randomForest) 
source('R/parameter_function.R')  
load('data/loc.feature.name.RData') 
load('data/model.self.rf.RData') 

#--load  and SCANB dataset, i.e. Data.SCANB includes the features of 231 samples, 
#where the first column is the sample name, 
#the second column is the HRD label,factor values with levels = c(1,0), 
#and the 3:53 columns are the considered 51 features.
load('data.SCANB.Rdata') #231 53

#--data conversion by scale(log(x+1)) 
matlog.SCANB=data.SCANB
matlog.SCANB[,c(3:53)]=scale(log(data.SCANB[,c(3:53)]+1))

#--establish the dataset for 5 models
mat.model.SCANB=list(model1=data.frame(matlog.SCANB[,match(c('id','is.hrd',loc.feature.name$model1),colnames(matlog.SCANB))]),
                     model2=data.frame(matlog.SCANB[,match(c('id','is.hrd',loc.feature.name$model2),colnames(matlog.SCANB))]),
                     model3=data.frame(matlog.SCANB[,match(c('id','is.hrd',loc.feature.name$model3),colnames(matlog.SCANB))]),
                     model4=data.frame(matlog.SCANB[,match(c('id','is.hrd',loc.feature.name$model4),colnames(matlog.SCANB))]),
                     model5=data.frame(matlog.SCANB[,match(c('id','is.hrd',loc.feature.name$model5),colnames(matlog.SCANB))]))

#--sequentially input the 5 models, and obtain the prediction probabilities and classifications for each model,
#where Prob.SCANB [m] is the prediction probability of the m-th model.
pred.SCANB=prob.SCANB=list() 
for(m in 1:5){ 
  prob.SCANB[[m]]=ssc:::predProb(model.self.rf[[m]]$model,data.matrix(mat.model.SCANB[[m]][,-c(1,2)]),
                                 model.self.rf[[m]]$pred,model.self.rf[[m]]$pred.pars)[,1]
  pred.SCANB[[m]]=ifelse(prob.SCANB[[m]]>=mat.cut.f1$cut[m],1,0) 
}
