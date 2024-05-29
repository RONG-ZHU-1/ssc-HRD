################################################################
# R version 4.2.3 (2023-03-15 ucrt)#############################
# Copyright (C) 2023 The R Foundation for Statistical Computing#
# Platform: x86_64-w64-mingw32/x64 (64-bit) ####################


################################################################
# construct the final ssc model ################################ 
################################################################

rm(list=ls()) 
setwd('...') 
set.seed(1) #

# library ####
library(proxy)#dist 
library(ssc)# 
library(randomForest)
library(pheatmap) #pheatmap
library(ggplot2) 
library(ggpubr)#stat_compare_means ,ggarrange


# SEC1. para,function ####
source('parameter_function.R')


# SEC2. data ####  
load("data/matlog1.all.RData")
load("data/loc.feature.name.RData")
load("data/loc.feature0.name.RData")


# SEC3. model on matlog1.all #### 
mat.model.all=list(model1=data.frame(matlog1.all[,match(c('id','is.hrd',loc.feature.name$model1),colnames(matlog1.all))]),
                   model2=data.frame(matlog1.all[,match(c('id','is.hrd',loc.feature.name$model2),colnames(matlog1.all))]),
                   model3=data.frame(matlog1.all[,match(c('id','is.hrd',loc.feature.name$model3),colnames(matlog1.all))]),
                   model4=data.frame(matlog1.all[,match(c('id','is.hrd',loc.feature.name$model4),colnames(matlog1.all))]),
                   model5=data.frame(matlog1.all[,match(c('id','is.hrd',loc.feature.name$model5),colnames(matlog1.all))]))
#
for(m in 1:n.model){
  mat.model.all[[m]][,2]=factor(mat.model.all[[m]][,2],levels = c(1,0))
}

#--
MAT.model.all=list()  
for(m in 1:n.model){  
  MAT.model.all[[m]]=list(model=mat.model.all[[m]],  
                          loc_l=which(!is.na(mat.model.all[[m]]$is.hrd)), #label  
                          loc_u=which(is.na(mat.model.all[[m]]$is.hrd)), #unlabel    
                          x=mat.model.all[[m]][,-c(1,2)], #1=id,2=paper.HRDetect is.hrd
                          y=mat.model.all[[m]][,2]) #2=is.hrd
  MAT.model.all[[m]]$x_l=MAT.model.all[[m]]$x[-MAT.model.all[[m]]$loc_u,] #l=label  
  MAT.model.all[[m]]$y_l=MAT.model.all[[m]]$y[-MAT.model.all[[m]]$loc_u]  #label
  MAT.model.all[[m]]$x_u=MAT.model.all[[m]]$x[MAT.model.all[[m]]$loc_u,]  #u=unlabel  
  MAT.model.all[[m]]$dist=as.matrix(proxy::dist(x = MAT.model.all[[m]]$x, method = "euclidean", by_rows = TRUE)) #distance matrices
  MAT.model.all[[m]]$kernel=as.matrix(exp(- 0.048 * MAT.model.all[[m]]$dist^2)) #kernel matrices
}


## ssc model: model.self.rf ####  
model.self.rf= list()
for(m in 1:n.model){
  print(m)
  model.self.rf[[m]]   =selfTraining(x = MAT.model.all[[m]]$x, y = MAT.model.all[[m]]$y, learner = randomForest, 
                                     pred = "predict",pred.pars =list(type = "prob"))
}  
#save(model.self.rf,file='data/model.self.rf.RData')
#load('data/model.self.rf.RData')


## importance of feature #### 
# Get Gini importance scores 
importance_scores=list()
for(m in 1:n.model){
  importance_scores[[m]]=importance(model.self.rf[[m]]$model, type = 2)
  importance_scores[[m]]=data.frame(feature=rownames(importance_scores[[m]]),
                                        Gini=importance_scores[[m]]) 
}  

#-- heatmap 
mat.imp.all=data.frame(feature=c(loc.feature0.name$cn1,
                                 loc.feature0.name$cn2,
                                 loc.feature0.name$snv), #30,11,9
                       block=c(rep(1,12),rep(2,9),rep(3,30))#2,3,1 for cna,ascn,snv,
)
for(m in 1:n.model){
  mat.imp.all=merge(mat.imp.all,importance_scores[[m]][,1:2],by.x='feature',by.y='feature',all.x=TRUE)
}
colnames(mat.imp.all)[3:7]=name.model
dim(mat.imp.all) #51  7
mat.imp.all=mat.imp.all[match(c(loc.feature0.name$cn1,
                                loc.feature0.name$cn2,
                                loc.feature0.name$snv),mat.imp.all$feature),]
rownames(mat.imp.all)=mat.imp.all$feature
mat.imp.all[is.na(mat.imp.all)]=0 #set NA to 0 for pheatmap 
annotation_row=data.frame(block=factor(mat.imp.all$block,levels = c(1:3),labels = c('CNA','ASCN','SNV')))
rownames(annotation_row)=mat.imp.all$feature
pheatmap(mat.imp.all[,3:7],cluster_row = FALSE, cluster_cols= FALSE, main ='Feature Importance',
         display_numbers = matrix(ifelse(mat.imp.all[,3:7]>0, "*", ""), nrow(mat.imp.all)),
         annotation_row=annotation_row,angle_col = "0",gaps_row = c(12,21),
         filename = "fig4_feature_importance_heatmap_all.pdf", 
         width =8, height =10, 
         cellwidth =50, cellheight =10) 



## bar plot of predicted probabilities ####
prob.all=list()   
for(m in 1:n.model){
  print(m)
  prob.all[[m]]=ssc:::predProb(model.self.rf[[m]]$model,data.matrix(MAT.model.all[[m]]$x),
                               model.self.rf[[m]]$pred,model.self.rf[[m]]$pred.pars)[,1] 
} 
#--
gg.prob.all.f1=plot.bar.5.given(prob.all,matlog1.all[,c('id','color1')],n.model=5,name.col=name.col,
                                   cut=mat.cut.f1$cut)  
#--
gg.bar.legend=get_legend(gg.prob.all.f1[[1]]+labs(fill = ""))
#ggarrange(gg.bar.legend)
#
pdf(file = paste0("figs4_bar_all_given.pdf"), width =8, height =8) #
ggarrange(gg.prob.all.f1[[1]]+ theme(legend.position = "none"),
          gg.prob.all.f1[[2]]+ theme(legend.position = "none"), 
          gg.prob.all.f1[[3]]+ theme(legend.position = "none"), 
          gg.prob.all.f1[[4]]+ theme(legend.position = "none"), 
          gg.prob.all.f1[[5]]+ theme(legend.position = "none"),
          gg.bar.legend,
          ncol=1,heights = c(rep(1,n.model),0.5) 
          )
dev.off() 



# save ####
save.image('output/ssc_model.RData') 
