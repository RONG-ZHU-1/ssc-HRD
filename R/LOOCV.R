
#################################
#### LOOCV process ##############
#################################

rm(list=ls()) 
setwd('..')  
set.seed(1) #


# library ####
library(doParallel)
library(foreach)
library(ssc)#
library(proxy)#dist 
library(randomForest)


# input #### 
#-- combined log-scale transformed data and list of selected features
load("data/matlog1.all.RData")
load("data/loc.feature.name.RData")
n.model=5



# doParallel ####
pred.loo=prob.loo=matrix(,dim(matlog1.all)[1],n.model) 

No.cpu=80

#-- parallel processing for LOOCV 
cl <- makeCluster(No.cpu) 
registerDoParallel(cl)

# 
LOO=foreach(ii=1:dim(matlog1.all)[1], .packages=c('ssc','randomForest')) %dopar%{ 
  ## data ####
  mat.model=list(model1=data.frame(matlog1.all[-ii,match(c('id','is.hrd',loc.feature.name$model1),colnames(matlog1.all))]),
                 model2=data.frame(matlog1.all[-ii,match(c('id','is.hrd',loc.feature.name$model2),colnames(matlog1.all))]),
                 model3=data.frame(matlog1.all[-ii,match(c('id','is.hrd',loc.feature.name$model3),colnames(matlog1.all))]),
                 model4=data.frame(matlog1.all[-ii,match(c('id','is.hrd',loc.feature.name$model4),colnames(matlog1.all))]),
                 model5=data.frame(matlog1.all[-ii,match(c('id','is.hrd',loc.feature.name$model5),colnames(matlog1.all))]))
  # 
  for(m in 1:n.model){
    mat.model[[m]][,2]=factor(mat.model[[m]][,2],levels = c(1,0))
  }
  #--
  MAT.model=list()  
  for(m in 1:n.model){  
    MAT.model[[m]]=list(model=mat.model[[m]],  
                        loc_l=which(!is.na(mat.model[[m]]$is.hrd)), #label  
                        loc_u=which(is.na(mat.model[[m]]$is.hrd)), #unlabel  
                        x=mat.model[[m]][,-c(1,2)], #1=id,2=is.hrd
                        y=mat.model[[m]][,2]) #2=is.hrd
    MAT.model[[m]]$x_l=MAT.model[[m]]$x[-MAT.model[[m]]$loc_u,] #l=label
    MAT.model[[m]]$y_l=MAT.model[[m]]$y[-MAT.model[[m]]$loc_u]  #label
    MAT.model[[m]]$x_u=MAT.model[[m]]$x[MAT.model[[m]]$loc_u,]  #u=unlabel 
    MAT.model[[m]]$dist=as.matrix(proxy::dist(x = MAT.model[[m]]$x, method = "euclidean", by_rows = TRUE)) #distance matrices
    MAT.model[[m]]$kernel=as.matrix(exp(- 0.048 * MAT.model[[m]]$dist^2)) #kernel matrices
  }
  
  ## model.self.rf.loo #### 
  model.self.rf.loo= list()
  for(m in 1:n.model){
    print(m)
    model.self.rf.loo[[m]]   =selfTraining(x = MAT.model[[m]]$x, y = MAT.model[[m]]$y, learner = randomForest, 
                                           pred = "predict",pred.pars =list(type = "prob"))
  }
  
  ## predict #### 
  mat.model.test=list(model1=data.frame(matlog1.all[ii,match(c('id','is.hrd',loc.feature.name$model1),colnames(matlog1.all))]),
                      model2=data.frame(matlog1.all[ii,match(c('id','is.hrd',loc.feature.name$model2),colnames(matlog1.all))]),
                      model3=data.frame(matlog1.all[ii,match(c('id','is.hrd',loc.feature.name$model3),colnames(matlog1.all))]),
                      model4=data.frame(matlog1.all[ii,match(c('id','is.hrd',loc.feature.name$model4),colnames(matlog1.all))]),
                      model5=data.frame(matlog1.all[ii,match(c('id','is.hrd',loc.feature.name$model5),colnames(matlog1.all))]))
  for(m in 1:n.model){
    mat.model.test[[m]][,2]=factor(mat.model.test[[m]][,2],levels = c(1,0))
  } 
  #
  for(m in 1:n.model){
    print(m)
    pred.loo[ii,m]=predict(model.self.rf.loo[[m]],data.matrix(mat.model.test[[m]][,-c(1,2)])) #Levels: 1 0
    prob.loo[ii,m]=ssc:::predProb(model.self.rf.loo[[m]]$model,
                                  data.matrix(mat.model.test[[m]][,-c(1,2)]),
                                  model.self.rf.loo[[m]]$pred,
                                  model.self.rf.loo[[m]]$pred.pars)[,1]
  }
  #--save
  list(pred.loo[ii,],prob.loo[ii,]) 
}  

# stop parallel programming 
stopCluster(cl) 

## output ####
save(LOO,file='output/LOOCV.RData')

