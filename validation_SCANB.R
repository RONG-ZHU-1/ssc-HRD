
#--load R packages, parameters, features' name and our five models.
set.seed(1) 
library(ssc) 
library(randomForest) 
library(reshape2)
library(ggplot2) 
library(ggpubr) 

source('parameter_function.R')  
load('data/loc.feature.name.RData') 
load('data/model.self.rf.RData') 

#--load  and SCANB dataset, i.e. Data.SCANB includes the features of 231 samples, 
#where the first column is the sample name, 
#the second column is the HRD label,factor values with levels = c(1,0), 
#and the 3:53 columns are the considered 51 features.
load('example/data.SCANB.Rdata') #231 53

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


#--boxplot
mat.prob.SCANB=cbind(prob.SCANB[[1]],prob.SCANB[[2]],prob.SCANB[[3]],prob.SCANB[[4]],prob.SCANB[[5]]) 
mat.prob.SCANB=data.frame(mat.prob.SCANB)
colnames(mat.prob.SCANB)=name.model 
mat.prob.SCANB$is.hrd=matlog.SCANB$is.hrd
mat.prob.SCANB=mat.prob.SCANB[!is.na(mat.prob.SCANB$is.hrd),]
mat.prob.SCANB$is.hrd=factor(mat.prob.SCANB$is.hrd,levels = c(1,0),labels= c('HRD+','HRD-')) 
mat.prob.SCANB=melt(mat.prob.SCANB)
dim(mat.prob.SCANB) 
colnames(mat.prob.SCANB)[c(1,3)]=c('HRD','prob')  

#--
gg.box.SCANB=ggplot(mat.prob.SCANB, aes(x=HRD, y=prob, fill=HRD)) +
  geom_boxplot()+stat_compare_means(method = "wilcox.test",label = "p.signif",color="blue")+
  theme_minimal()+
  xlab('')+ ylab('Predicted probabilities')+
  facet_wrap(~variable,ncol=5)+ 
  geom_hline(data = mat.cut.f1, aes(yintercept = cut), linetype="dashed", color = "grey") 
gg.box.SCANB
