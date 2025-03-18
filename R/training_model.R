################################################################
# R version 4.2.3 (2023-03-15 ucrt)#############################
# Copyright (C) 2023 The R Foundation for Statistical Computing#
# Platform: x86_64-w64-mingw32/x64 (64-bit) ####################


################################################################
# feature selection, training model and performance of LOOCV ### 
################################################################

rm(list=ls()) 
setwd('...') 
set.seed(1) #

# library ####
library(reshape2)#melt
library(ggplot2) 
library(ggpubr)#stat_compare_means 
library(Hmisc) #rcorr
library(corrplot)#corrplot  
library(doParallel)
library(foreach)
library(ssc)#
library(proxy)#dist 
library(randomForest)
library(caret) #Classification and Regression Training ; confusionMatrix()
library(mltools) #calculate Matthews correlation coefficient, mcc(preds, actual)
library(PRROC) #pr.curve,roc.curve
library(ggpubr) #ggarrange  
library(ROCR) #Visualizing the Performance of Scoring Classifiers;performance(), prediction() 



# SEC1. parameter,function ####
source('parameter_function.R')

# SEC2. data ####   
load('data/mat.TCGA.RData')  #1003 85
load('data/mat.ICGC.RData')    #401 90 

## location #### 
#-- 
loc.TCGA=list(hrd=match('is.hrd',colnames(mat.TCGA)),
              SBS.v3.3.brca=match('SBS1.brca',colnames(mat.TCGA)):match('SBS30.brca',colnames(mat.TCGA)),
              ID.v3.3=match('ID1',colnames(mat.TCGA)):match('ID18',colnames(mat.TCGA)),
              cna.6score=match(c('cnaBurden','cnaLoad',"MACN",'TDP','TDP_size','Chromoth_state'),colnames(mat.TCGA)),
              CX.brca=c(match('CX1.brca',colnames(mat.TCGA)):match('CX9.brca',colnames(mat.TCGA))), 
              CN.brca=match('CN1.brca',colnames(mat.TCGA)):match('CN17.brca',colnames(mat.TCGA)),
              scar=match(c('scar.loh','scar.lst','scar.tai','scar.sum'),colnames(mat.TCGA)),
              color1=match('color1',colnames(mat.TCGA)),
              cohort=match('cohort',colnames(mat.TCGA)) 
)  
#-- 
loc.ICGC=list(hrd=match('is.hrd',colnames(mat.ICGC)),
              SBS.v3.3.brca=match('SBS1.brca',colnames(mat.ICGC)):match('SBS30.brca',colnames(mat.ICGC)),
              ID.v3.3=match('ID1',colnames(mat.ICGC)):match('ID18',colnames(mat.ICGC)),
              cna.6score=match(c('cnaBurden','cnaLoad',"MACN",'TDP','TDP_size','Chromoth_state'),colnames(mat.ICGC)),
              CX.brca=c(match('CX1.brca',colnames(mat.ICGC)):match('CX9.brca',colnames(mat.ICGC))),
              CN.brca=match('CN1.brca',colnames(mat.ICGC)):match('CN17.brca',colnames(mat.ICGC)),
              scar=match(c('scar.loh','scar.lst','scar.tai','scar.sum'),colnames(mat.ICGC)),
              color1=match('color1',colnames(mat.ICGC)),
              cohort=match('cohort',colnames(mat.ICGC))
)


## matlog.x=scale(log(x+1)) ####   
matlog.TCGA=mat.TCGA 
for(i in unique(as.numeric(unlist(loc.TCGA[c(2:7)])))){ 
  if(sum(na.omit(matlog.TCGA[,i]))==0|sum(na.omit(matlog.TCGA[,i])<0)>0){
    matlog.TCGA[,i]=NA
  }else{
    matlog.TCGA[,i]=scale(log(mat.TCGA[,i]+1))
  }
}
#
matlog.ICGC=mat.ICGC 
for(i in unique(as.numeric(unlist(loc.ICGC[c(2:7)])))){ 
  if(sum(na.omit(matlog.ICGC[,i]))==0|sum(na.omit(matlog.ICGC[,i])<0)>0){
    matlog.ICGC[,i]=NA
  }else{
    matlog.ICGC[,i]=scale(log(mat.ICGC[,i]+1))
  }
} 

### matlog1.TCGA, n*55, matlog1.all=MTI #### 
matlog1.TCGA=matlog.TCGA[,c(1,loc.TCGA$hrd,loc.TCGA$SBS.v3.3.brca,loc.TCGA$ID.v3.3,
                            loc.TCGA$cna.6score,loc.TCGA$CX.brca,
                            loc.TCGA$CN.brca,loc.TCGA$scar[4],loc.TCGA$color1,loc.TCGA$cohort)]
matlog1.ICGC=matlog.ICGC[,c(1,loc.ICGC$hrd,loc.ICGC$SBS.v3.3.brca,loc.ICGC$ID.v3.3,
                            loc.ICGC$cna.6score,loc.ICGC$CX.brca,
                            loc.ICGC$CN.brca,loc.ICGC$scar[4],loc.ICGC$color1,loc.ICGC$cohort)]
 
#--combine
matlog1.all=rbind(matlog1.TCGA,
                  matlog1.ICGC)
dim(matlog1.all)
#save(matlog1.all,file='data/matlog1.all.RData')

 
## location of all 51 features;loc.feature0,12,9,30=51 #### 
loc.feature0=list(cn1=33:44, 
                  cn2=45:53,
                  snv=3:32) #
loc.feature0$all=unlist(loc.feature0) #
loc.feature0$model1=loc.feature0$cn1
loc.feature0$model2=c(loc.feature0$cn1,loc.feature0$cn2)
loc.feature0$model3=loc.feature0$snv
loc.feature0$model4=c(loc.feature0$snv,loc.feature0$cn1)
loc.feature0$model5=loc.feature0$all 
#--
loc.feature0.name=list(colnames(matlog1.all)[loc.feature0[[1]]],
                       colnames(matlog1.all)[loc.feature0[[2]]],
                       colnames(matlog1.all)[loc.feature0[[3]]],
                       colnames(matlog1.all)[loc.feature0[[4]]],
                       colnames(matlog1.all)[loc.feature0[[5]]],
                       colnames(matlog1.all)[loc.feature0[[6]]],
                       colnames(matlog1.all)[loc.feature0[[7]]],
                       colnames(matlog1.all)[loc.feature0[[8]]],
                       colnames(matlog1.all)[loc.feature0[[9]]])
names(loc.feature0.name)=names(loc.feature0) 
#save(loc.feature0.name,file='data/loc.feature0.name.RData')



### gg.box.feature, feature~hrd on matlog1.all ####
#--
df.box.feature.cn1=matlog1.all[!is.na(matlog1.all$is.hrd),c(2,loc.feature0$cn1)] #
for(i in 2:dim(df.box.feature.cn1)[2]){
  df.box.feature.cn1[,i]=as.numeric(df.box.feature.cn1[,i])
} 
df.box.feature.cn1=melt(df.box.feature.cn1,1)#
dim(df.box.feature.cn1)#52230     3
colnames(df.box.feature.cn1)=c('HRD','feature','value')
df.box.feature.cn1$HRD=ifelse(df.box.feature.cn1$HRD==1,'HRD+','HRD-')
df.box.feature.cn1$block='CNA'
#
df.box.feature.cn2=matlog1.all[!is.na(matlog1.all$is.hrd),c(2,loc.feature0$cn2)] #
for(i in 2:dim(df.box.feature.cn2)[2]){
  df.box.feature.cn2[,i]=as.numeric(df.box.feature.cn2[,i])
} 
df.box.feature.cn2=melt(df.box.feature.cn2,1)#
dim(df.box.feature.cn2)#52230     3
colnames(df.box.feature.cn2)=c('HRD','feature','value')
df.box.feature.cn2$HRD=ifelse(df.box.feature.cn2$HRD==1,'HRD+','HRD-')
df.box.feature.cn2$block='ASCN'
#
df.box.feature.snv=matlog1.all[!is.na(matlog1.all$is.hrd),c(2,loc.feature0$snv)] #
for(i in 2:dim(df.box.feature.snv)[2]){
  df.box.feature.snv[,i]=as.numeric(df.box.feature.snv[,i])
} 
df.box.feature.snv=melt(df.box.feature.snv,1)#
dim(df.box.feature.snv)#52230     3
colnames(df.box.feature.snv)=c('HRD','feature','value')
df.box.feature.snv$HRD=ifelse(df.box.feature.snv$HRD==1,'HRD+','HRD-')
df.box.feature.snv$block='SNV'
#
df.box.feature.cn1$HRD=factor(df.box.feature.cn1$HRD,levels = c('HRD+','HRD-'))
df.box.feature.cn2$HRD=factor(df.box.feature.cn2$HRD,levels = c('HRD+','HRD-'))
df.box.feature.snv$HRD=factor(df.box.feature.snv$HRD,levels = c('HRD+','HRD-'))
#--
gg.box.feature.cn1=ggplot(df.box.feature.cn1,aes(x=HRD,y=value,fill=HRD))+  
  facet_wrap(~feature,scales = "free",ncol=6)+  
  geom_boxplot()+ylab('')+ xlab('')+ 
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.y.npc = 0.8,col='blue')+  
  guides(fill="none") 
#
gg.box.feature.cn2=ggplot(df.box.feature.cn2,aes(x=HRD,y=value,fill=HRD))+  
  facet_wrap(~feature,scales = "free",ncol=6)+ #,ncol=1 
  geom_boxplot()+ylab('')+ xlab('')+ 
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.y.npc = 0.8,col='blue')+  
  guides(fill="none") 
#
gg.box.feature.snv=ggplot(df.box.feature.snv,aes(x=HRD,y=value,fill=HRD))+  
  facet_wrap(~feature,scales = "free",ncol=6)+ 
  geom_boxplot()+ylab('')+ xlab('')+ 
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.y.npc = 0.8,col='blue')+  
  guides(fill="none") 

#--Supplementary Figure 1
pdf('fig_box_feature_3block.pdf', width = 10, height = 11) 
ggarrange(NULL,NULL,gg.box.feature.cn1,
          NULL,gg.box.feature.cn2,
          NULL,gg.box.feature.snv, 
          labels=c(NA,name.model[1],NA,name.model[2],NA,name.model[3]),
          ncol=1,nrow=7, heights = c(0.1,0.2,2,0.2,2,0.2,5))  
dev.off()



### wilcox.test for feature~hrd on matlog1.all #### 
feature.wilcox.pval.all=c()
for (i in 1:length(loc.feature0$all)) {
  feature.wilcox.pval.all=c(feature.wilcox.pval.all,wilcox.test(matlog1.all[,loc.feature0$all[i]]~matlog1.all$is.hrd)$p.value)
}
names(feature.wilcox.pval.all)=colnames(matlog1.all)[loc.feature0$all]



#### delete p>0.05, loc.feature,10,7,21  ####
names(feature.wilcox.pval.all)[which(feature.wilcox.pval.all<=0.05)] #38  
loc.feature.keepp=match(names(feature.wilcox.pval.all)[which(feature.wilcox.pval.all<=0.05)],colnames(matlog1.all)) #38

loc.feature=list(
  cn1=match(intersect(colnames(matlog1.all)[loc.feature0$cn1],colnames(matlog1.all)[loc.feature.keepp]),colnames(matlog1.all)),
  cn2=match(intersect(colnames(matlog1.all)[loc.feature0$cn2],colnames(matlog1.all)[loc.feature.keepp]),colnames(matlog1.all)),
  snv=match(intersect(colnames(matlog1.all)[loc.feature0$snv],colnames(matlog1.all)[loc.feature.keepp]),colnames(matlog1.all))
) #10,7,21


#### pearson cor>0.8,delete TDP_size,keep cnaLoad/CX5.brca;loc.feature=37=9,7,21 #### 
#--cn1
mat=na.omit(matlog1.all[,loc.feature$cn1])  
mat.cor.cn1=rcorr(as.matrix(mat)) 
mat.cor.cn1$loc0.8=which(abs(mat.cor.cn1$r)>=0.8,arr.ind = T)
mat.cor.cn1$loc0.8=mat.cor.cn1$loc0.8[which(mat.cor.cn1$loc0.8[,1]!=mat.cor.cn1$loc0.8[,2]),] 
mat.cor.cn1$loc0.8 

#--cn2
mat=na.omit(matlog1.all[,loc.feature$cn2])  
mat.cor.cn2=rcorr(as.matrix(mat)) 
mat.cor.cn2$loc0.8=which(abs(mat.cor.cn2$r)>=0.8,arr.ind = T)
mat.cor.cn2$loc0.8=mat.cor.cn2$loc0.8[which(mat.cor.cn2$loc0.8[,1]!=mat.cor.cn2$loc0.8[,2]),] 
mat.cor.cn2$loc0.8 #0


#--snv
mat=na.omit(matlog1.all[,loc.feature$snv])  
mat.cor.snv=rcorr(as.matrix(mat)) 
mat.cor.snv$loc0.8=which(abs(mat.cor.snv$r)>=0.8,arr.ind = T)
mat.cor.snv$loc0.8=mat.cor.snv$loc0.8[which(mat.cor.snv$loc0.8[,1]!=mat.cor.snv$loc0.8[,2]),] 
mat.cor.snv$loc0.8 

#--delete TDP_size,SBS2.brca
loc.feature$cn1=loc.feature$cn1[-4]
loc.feature$snv=loc.feature$snv[-1]
loc.feature$all=unlist(loc.feature)#37,30=9,5,16  
loc.feature$model1=loc.feature$cn1
loc.feature$model2=c(loc.feature$cn1,loc.feature$cn2)
loc.feature$model3=loc.feature$snv
loc.feature$model4=c(loc.feature$snv,loc.feature$cn1)
loc.feature$model5=loc.feature$all 

#--save
loc.feature.name=list(colnames(matlog1.all)[loc.feature[[1]]],
                       colnames(matlog1.all)[loc.feature[[2]]],
                       colnames(matlog1.all)[loc.feature[[3]]],
                       colnames(matlog1.all)[loc.feature[[4]]],
                       colnames(matlog1.all)[loc.feature[[5]]],
                       colnames(matlog1.all)[loc.feature[[6]]],
                       colnames(matlog1.all)[loc.feature[[7]]],
                       colnames(matlog1.all)[loc.feature[[8]]],
                       colnames(matlog1.all)[loc.feature[[9]]])
names(loc.feature.name)=names(loc.feature)
#save(loc.feature.name,file='data/loc.feature.name.RData')


### cor plot #### 
mat.select.cor.all=cor(matlog1.all[,loc.feature$model5])  
mat.select.pval.all <- cor.mtest(matlog1.all[,loc.feature$model5])
#
i1=9;i2=14;i3=30
#--Supplementary Figure 6 
pdf(file = paste0("fig_cor_select_3block_between12.pdf"), width=5, height =8) 
corrplot(mat.select.cor.all[1:i1,(i1+1):i2], method="color",   
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = mat.select.pval.all[1:i1,(i1+1):i2], sig.level = 0.05, insig = "blank",
         col=colorRampPalette(c("blue","white","red"))(100)
)
dev.off()
#
pdf(file = paste0("fig_cor_select_3block_between13.pdf"), width=11, height =8)#
corrplot(mat.select.cor.all[1:i1,(i2+1):i3], method="color",   
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = mat.select.pval.all[1:i1,(i2+1):i3], sig.level = 0.05, insig = "blank",
         col=colorRampPalette(c("blue","white","red"))(100) 
)
dev.off()
#
pdf(file = paste0("fig_cor_select_3block_between23.pdf"), width=11, height =5)#
corrplot(mat.select.cor.all[(i1+1):i2,(i2+1):i3], method="color", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = mat.select.pval.all[(i1+1):i2,(i2+1):i3], sig.level = 0.05, insig = "blank",
         col=colorRampPalette(c("blue","white","red"))(100) 
)#
dev.off()



#--
pdf(file = paste0("fig_cor_select_cn1.pdf"), width =8, height =8) 
corrplot(mat.select.cor.all[1:i1,1:i1], method="color",  
         type="upper",  
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = mat.select.pval.all[1:i1,1:i1], sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE, 
         col=colorRampPalette(c("blue","white","red"))(100) 
) 
dev.off()
#
pdf(file = paste0("fig_cor_select_cn2.pdf"), width =5, height =5)#
corrplot(mat.select.cor.all[(i1+1):i2,(i1+1):i2], method="color", 
         type="upper",  
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = mat.select.pval.all[(i1+1):i2,(i1+1):i2], sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,  
         col=colorRampPalette(c("blue","white","red"))(100))
dev.off()
#
pdf(file = paste0("fig_cor_select_snv.pdf"), width =10, height =10)#
corrplot(mat.select.cor.all[(i2+1):i3,(i2+1):i3], method="color", 
         type="upper",  
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = mat.select.pval.all[(i2+1):i3,(i2+1):i3], sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE, 
         col=colorRampPalette(c("blue","white","red"))(100)
)
dev.off()





### cor(HRDetect,feature) in ICGC #### 
#--from HRDetect
HRDetect.features_mean <- c(0.218,2.096,1.260,1.935,2.195,4.390)
HRDetect.features_sd <- c(0.090,3.555,1.657,1.483,0.750,3.179)
mat.ICGC.HRDetect.feature=mat.ICGC[,c('paper.del.mh.prop','paper.e.3','paper.SV3','paper.SV5','paper.hrd','paper.e.8','scar.loh')]
which(is.na(mat.ICGC.HRDetect.feature),arr.ind = TRUE) #4 hrd is NA, match scar.loh
mat.ICGC.HRDetect.feature$paper.hrd[is.na(mat.ICGC.HRDetect.feature$paper.hrd)]=mat.ICGC.HRDetect.feature$scar.loh[is.na(mat.ICGC.HRDetect.feature$paper.hrd)]
mat.ICGC.HRDetect.feature=mat.ICGC.HRDetect.feature[,1:6]
colnames(mat.ICGC.HRDetect.feature)=substr(colnames(mat.ICGC.HRDetect.feature),7,20)
for(i in 1:dim(mat.ICGC.HRDetect.feature)[2]){
  mat.ICGC.HRDetect.feature[,i]=(log(mat.ICGC.HRDetect.feature[,i]+1)-HRDetect.features_mean[i])/HRDetect.features_sd[i]
}

#--
mat.cor.ICGC=cor(cbind(matlog1.ICGC[,loc.feature$model5],
                       mat.ICGC.HRDetect.feature))  
mat.cor.ICGC.pval <- cor.mtest(cbind(matlog1.ICGC[,loc.feature$model5],
                                     mat.ICGC.HRDetect.feature))
#
i1=9;i2=14;i3=30;i4=36
#
pdf(file = paste0("fig_cor_ICGC_HRDetect_v1.pdf"), width=12, height =8)#
par(mfrow=c(2,1)) 
corrplot(mat.cor.ICGC[(i3+1):i4,1:i2], method="color",  
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = mat.cor.ICGC.pval[(i3+1):i4,1:i2], sig.level = 0.05, insig = "blank",
         col=colorRampPalette(c("blue","white","red"))(100), 
         mar=c(0,0,3,0),
         main ='ASCN'
) 
corrplot(mat.cor.ICGC[(i3+1):i4,(i2+1):i3], method="color",  
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = mat.cor.ICGC.pval[(i3+1):i4,(i2+1):i3], sig.level = 0.05, insig = "blank",
         col=colorRampPalette(c("blue","white","red"))(100), 
         mar=c(1,0,2,0),
         main ='SNV'
)#
dev.off()


## mat.v1.clinical ####
#--
mat.v1.clinical.TCGA=mat.TCGA[,match(c('id','is.hrd',name.clinical,'cohort'),colnames(mat.TCGA))]
mat.v1.clinical.TCGA$is.hrd=as.character(mat.v1.clinical.TCGA$is.hrd)
mat.v1.clinical.TCGA$is.hrd[is.na(mat.v1.clinical.TCGA$is.hrd)]='NA'
mat.v1.clinical.TCGA$is.hrd=factor(mat.v1.clinical.TCGA$is.hrd,
                                   levels =c('1','0','NA'),
                                   labels =c('HRD+','HRD-','HRD unknown'))
#--
mat.v1.clinical.ICGC=mat.ICGC[,match(c('id','is.hrd',name.clinical,'cohort'),colnames(mat.ICGC))]
mat.v1.clinical.ICGC$is.hrd=as.character(mat.v1.clinical.ICGC$is.hrd)
mat.v1.clinical.ICGC$is.hrd[is.na(mat.v1.clinical.ICGC$is.hrd)]='NA'
mat.v1.clinical.ICGC$is.hrd=factor(mat.v1.clinical.ICGC$is.hrd,
                                   levels =c('1','0','NA'),
                                   labels =c('HRD+','HRD-','HRD unknown'))




### plot clinical ####    
name.type.v1.n=c(2,2,2,2,3,2,5,10)
name.type.v1.label=c('<50','>=50',
                     'ER+','ER-',
                     'Her2+','Her2-',
                     'PR+','PR-',
                     "I","II","III",
                     'LN+','LN-',
                     'Ba','H','A','B','N',
                     1:10)
mat.type.v1=data.frame(type=c(rep(name.clinical,name.type.v1.n),
                              rep(name.clinical,name.type.v1.n),
                              rep(name.clinical,name.type.v1.n)),  
                       label=c(name.type.v1.label,name.type.v1.label,name.type.v1.label),#80
                       HRD=c(rep('HRD+',length(name.type.v1.label)),
                             rep('HRD-',length(name.type.v1.label)),
                             rep('HRD unknown',length(name.type.v1.label)))) 
#add prob; 
mat.type.v1$TCGA=c(prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD+',3])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD+',4])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD+',5])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD+',6])), 
                   NA,NA,NA, 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD+',8])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD+',9])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD+',10])), #HRD+ 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD-',3])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD-',4])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD-',5])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD-',6])), 
                   NA,NA,NA, 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD-',8])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD-',9])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD-',10])), #HRD-
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD unknown',3])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD unknown',4])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD unknown',5])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD unknown',6])), 
                   NA,NA,NA, 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD unknown',8])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD unknown',9])), 
                   prop.table(table(mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd=='HRD unknown',10])) #HRD unknown
) 
mat.type.v1$ICGC=c(prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD+',3])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD+',4])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD+',5])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD+',6])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD+',7])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD+',8])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD+',9])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD+',10])),  #HRD+ 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD-',3])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD-',4])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD-',5])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD-',6])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD-',7])),
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD-',8])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD-',9])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD-',10])), #HRD-
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD unknown',3])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD unknown',4])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD unknown',5])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD unknown',6])),
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD unknown',7])),
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD unknown',8])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD unknown',9])), 
                   prop.table(table(mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd=='HRD unknown',10]))#HRD unknown
) 

#--- 
mat.type.v1$type=factor(mat.type.v1$type,levels = name.clinical) #8 
mat.type.v1$label=factor(mat.type.v1$label,levels =name.type.v1.label)
mat.type.v1$HRD=factor(mat.type.v1$HRD,levels =c('HRD+','HRD-','HRD unknown'))
# 
mat.type.v1$color=rep(c("#9ac9db",'#ff8884',
                        "#66C2A5",'#f8ac8c',
                        "#82B0D2",'#ff8884',
                        "#66C2A5",'#f8ac8c',
                        '#CFEAF1','#B883D4','#FA7F6F',
                        "#66C2A5",'#f8ac8c', 
                        colipam,
                        coliClust),3)
#
dim(mat.type.v1)


baseSize=13
#--TCGA 
HRD.labs=paste0(c('HRD+ (N=','HRD- (N=','HRD unknown (N='),
                as.numeric(table(mat.v1.clinical.TCGA$is.hrd)),')')
names(HRD.labs)=c('HRD+','HRD-','HRD unknown')
#
gg.sub.v1.TCGA=ggplot(mat.type.v1[!is.na(mat.type.v1$TCGA),],aes(y=TCGA,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) +  
  scale_fill_manual(values=mat.type.v1[!is.na(mat.type.v1$TCGA),'color'])+ 
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)) +
  facet_grid(~HRD,labeller = labeller(HRD=HRD.labs)) + 
  labs(y="Proportion of cases")+
  scale_x_discrete(expand = c(0, 0), breaks=name.clinical,
                   labels=name.clinical)+ 
  theme_grey(base_size = baseSize) + 
  theme(
    axis.text=element_text(size = baseSize*0.7, face="bold",colour = "black"),
    axis.title=element_blank(),
    axis.ticks=element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    panel.background = element_blank(),
    strip.text.x = element_text(size = baseSize*0.8, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.00001, "lines"),
    panel.spacing.y = unit(0.8, "lines"))+
  coord_flip() 

#--ICGC  
HRD.labs=paste0(c('HRD+ (N=','HRD- (N=','HRD unknown (N='),
                as.numeric(table(mat.v1.clinical.ICGC$is.hrd)),')')
names(HRD.labs)=c('HRD+','HRD-','HRD unknown')
#
gg.sub.v1.ICGC=ggplot(mat.type.v1[!is.na(mat.type.v1$ICGC),],aes(y=ICGC,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) +  
  scale_fill_manual(values=mat.type.v1[!is.na(mat.type.v1$ICGC),'color'])+ 
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)) +
  facet_grid(~HRD,labeller = labeller(HRD=HRD.labs)) +  
  labs(y="Proportion of cases")+
  scale_x_discrete(expand = c(0, 0), breaks=name.clinical,
                   labels=name.clinical)+ 
  theme_grey(base_size = baseSize) + 
  theme(
    axis.text=element_text(size = baseSize*0.7, face="bold",colour = "black"),
    axis.title=element_blank(),
    axis.ticks=element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    panel.background = element_blank(),
    strip.text.x = element_text(size = baseSize*0.8, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.00001, "lines"),
    panel.spacing.y = unit(0.8, "lines"))+
  coord_flip()
gg.sub.v1.ICGC

#--
pdf(file = "fig_clinical.pdf")
ggarrange(gg.sub.v1.TCGA,gg.sub.v1.ICGC, 
          labels=c('TCGA','ICGC'),ncol=1,nrow=2) 
dev.off()



### chisq.test ####  
mat.v2.clinical.TCGA=mat.v1.clinical.TCGA[mat.v1.clinical.TCGA$is.hrd%in%c('HRD+','HRD-'),]
mat.v2.clinical.TCGA$is.hrd=factor(mat.v2.clinical.TCGA$is.hrd,levels =c('HRD+','HRD-'))
#
mat.v2.clinical.ICGC=mat.v1.clinical.ICGC[mat.v1.clinical.ICGC$is.hrd%in%c('HRD+','HRD-'),]
mat.v2.clinical.ICGC$is.hrd=factor(mat.v2.clinical.ICGC$is.hrd,levels =c('HRD+','HRD-'))

#--
mat.v2.clinical.test=data.frame(type=name.clinical)
mat.v2.clinical.test$TCGA=round(c(chisq.test(table(mat.v2.clinical.TCGA[,c(3,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.TCGA[,c(4,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.TCGA[,c(5,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.TCGA[,c(6,2)]))$p.value,
                                  NA,
                                  chisq.test(table(mat.v2.clinical.TCGA[,c(8,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.TCGA[,c(9,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.TCGA[,c(10,2)]))$p.value),3) 
#
mat.v2.clinical.test$ICGC=round(c(chisq.test(table(mat.v2.clinical.ICGC[,c(3,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.ICGC[,c(4,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.ICGC[,c(5,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.ICGC[,c(6,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.ICGC[,c(7,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.ICGC[,c(8,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.ICGC[,c(9,2)]))$p.value,
                                  chisq.test(table(mat.v2.clinical.ICGC[,c(10,2)]))$p.value),3)

# SEC3. see LOOCV.R to obtain LOOCV.RData ####


# SEC4. performance ####
## load('LOOCV.RData') ####
load('output/LOOCV.RData')
#
loo.pred=t(array(unlist(lapply(LOO,"[[",1)),dim=c(5,1404))) 
loo.pred=data.frame(loo.pred)
colnames(loo.pred)=name.model 
loo.pred=2-loo.pred
#
loo.prob=t(array(unlist(lapply(LOO,"[[",2)),dim=c(5,1404)))
loo.prob=data.frame(loo.prob)
colnames(loo.prob)=name.model 


## select f1 cut ####
perfm.loo=list()  
for(m in 1:n.model){
  print(m)
  perfm.loo[[m]]=perf.prob.f1(loo.prob[which(!is.na(matlog1.all$is.hrd)),m],
                              factor(matlog1.all$is.hrd[which(!is.na(matlog1.all$is.hrd))],levels =c(1,0)))#, should levels=c(1,0) 
} #m


#--result
Result.loo=matrix(,length(name.perf.v1),n.model)
for(m in 1:n.model){ 
  Result.loo[,m]=perfm.loo[[m]]$perf 
}
Result.loo=data.frame(Result.loo) 
Result.loo$index=name.perf.v1  
dim(Result.loo)  
colnames(Result.loo)[1:n.model]=name.model


#--optimal threshold
mat.cut.f1=data.frame(variable=name.model,cut=as.numeric(Result.loo[7,1:5]))
mat.cut.f1$variable=factor(mat.cut.f1$variable,levels=name.model)
mat.cut.f1


## roc/pr #### 
loo.roc=loo.pr=list() 
  for(m in 1:n.model){ 
    loo.roc[[m]]=roc.curve(loo.prob[which(matlog1.all$is.hrd==1),m],
                             loo.prob[which(matlog1.all$is.hrd==0),m], curve = TRUE,rand.compute=T)
    loo.pr[[m]]=pr.curve(loo.prob[which(matlog1.all$is.hrd==1),m],
                          loo.prob[which(matlog1.all$is.hrd==0),m], curve = TRUE,rand.compute=T)
    }
names(loo.roc)=name.model#paste0('model',1:n.model) 
names(loo.pr)=name.model#paste0('model',1:n.model) 
#
loo.roc.TCGA=loo.pr.TCGA=list() 
for(m in 1:n.model){ 
  loo.roc.TCGA[[m]]=roc.curve(loo.prob[which(matlog1.all$is.hrd==1&matlog1.all$cohort=='TCGA'),m],
                              loo.prob[which(matlog1.all$is.hrd==0&matlog1.all$cohort=='TCGA'),m], curve = TRUE,rand.compute=T)
  loo.pr.TCGA[[m]]=pr.curve(loo.prob[which(matlog1.all$is.hrd==1&matlog1.all$cohort=='TCGA'),m],
                            loo.prob[which(matlog1.all$is.hrd==0&matlog1.all$cohort=='TCGA'),m], curve = TRUE,rand.compute=T)
}
names(loo.roc.TCGA)=name.model#paste0('model',1:n.model) 
names(loo.pr.TCGA)=name.model#paste0('model',1:n.model) 
#
loo.roc.ICGC=loo.pr.ICGC=list() 
for(m in 1:n.model){ 
  loo.roc.ICGC[[m]]=roc.curve(loo.prob[which(matlog1.all$is.hrd==1&matlog1.all$cohort=='ICGC'),m],
                              loo.prob[which(matlog1.all$is.hrd==0&matlog1.all$cohort=='ICGC'),m], curve = TRUE,rand.compute=T)
  loo.pr.ICGC[[m]]=pr.curve(loo.prob[which(matlog1.all$is.hrd==1&matlog1.all$cohort=='ICGC'),m],
                            loo.prob[which(matlog1.all$is.hrd==0&matlog1.all$cohort=='ICGC'),m], curve = TRUE,rand.compute=T)
}
names(loo.roc.ICGC)=name.model#paste0('model',1:n.model) 
names(loo.pr.ICGC)=name.model#paste0('model',1:n.model) 

 
#--
pdf(file = paste0("fig_loo_roc.pdf"))
par(mfrow=c(1,3))
plot.roc.5(loo.roc,pos="bottomright",cohort='Combined')  
plot.roc.5(loo.roc.TCGA,pos="bottomright",cohort='TCGA')
plot.roc.5(loo.roc.ICGC,pos="bottomright",cohort='ICGC') 
dev.off()

#--
pdf(file = paste0("fig_loo_pr.pdf"))
par(mfrow=c(1,3))
plot.pr.5(loo.pr,pos="bottomleft",cohort='Combined')  
plot.pr.5(loo.pr.TCGA,pos="bottomleft",cohort='TCGA') 
plot.pr.5(loo.pr.ICGC,pos="bottomleft",cohort='ICGC') 
dev.off()



## boxplot ####
#
mat.prob.loo=loo.prob
mat.prob.loo$is.hrd=matlog1.all$is.hrd
mat.prob.loo=mat.prob.loo[!is.na(mat.prob.loo$is.hrd),]
mat.prob.loo$is.hrd=factor(mat.prob.loo$is.hrd,levels = c(1,0),labels= c('HRD+','HRD-'))
mat.prob.loo=melt(mat.prob.loo) 
colnames(mat.prob.loo)[c(1,3)]=c('HRD','prob')  
#
mat.prob.loo.TCGA=loo.prob[which(matlog1.all$cohort=='TCGA'),]
mat.prob.loo.TCGA$is.hrd=matlog1.TCGA$is.hrd
mat.prob.loo.TCGA=mat.prob.loo.TCGA[!is.na(mat.prob.loo.TCGA$is.hrd),]
mat.prob.loo.TCGA$is.hrd=factor(mat.prob.loo.TCGA$is.hrd,levels = c(1,0),labels= c('HRD+','HRD-'))
mat.prob.loo.TCGA=melt(mat.prob.loo.TCGA) 
colnames(mat.prob.loo.TCGA)[c(1,3)]=c('HRD','prob')  
#
mat.prob.loo.ICGC=loo.prob[which(matlog1.all$cohort=='ICGC'),]
mat.prob.loo.ICGC$is.hrd=matlog1.ICGC$is.hrd
mat.prob.loo.ICGC=mat.prob.loo.ICGC[!is.na(mat.prob.loo.ICGC$is.hrd),]
mat.prob.loo.ICGC$is.hrd=factor(mat.prob.loo.ICGC$is.hrd,levels = c(1,0),labels= c('HRD+','HRD-'))
mat.prob.loo.ICGC=melt(mat.prob.loo.ICGC) 
colnames(mat.prob.loo.ICGC)[c(1,3)]=c('HRD','prob') 

#--
gg.box.loo=ggplot(mat.prob.loo, aes(x=HRD, y=prob, fill=HRD)) + 
  geom_boxplot()+stat_compare_means(method = "wilcox.test",label = "p.signif",color="blue",
                                    label.y.npc = 0.9,hide.ns = T)+  
  xlab('')+ylab('')+  
  facet_wrap(~variable,ncol=5)+ 
  geom_hline(data = mat.cut.f1, aes(yintercept = cut), linetype="dashed", color = "grey")
#
gg.box.loo.TCGA=ggplot(mat.prob.loo.TCGA, aes(x=HRD, y=prob, fill=HRD)) +  
  geom_boxplot()+stat_compare_means(method = "wilcox.test",label = "p.signif",color="blue",
                                    label.y.npc = 0.9)+   
  xlab('')+ylab('')+ 
  facet_wrap(~variable,ncol=5)+ 
  geom_hline(data = mat.cut.f1, aes(yintercept = cut), linetype="dashed", color = "grey") 
#
gg.box.loo.ICGC=ggplot(mat.prob.loo.ICGC, aes(x=HRD, y=prob, fill=HRD)) +  
  geom_boxplot()+stat_compare_means(method = "wilcox.test",label = "p.signif",color="blue",
                                    label.y.npc = 0.9)+ 
  xlab('')+ylab('')+ 
  facet_wrap(~variable,ncol=5)+  
  geom_hline(data = mat.cut.f1, aes(yintercept = cut), linetype="dashed", color = "grey")

#--
pdf(file ="fig_box_loo.pdf")
ggarrange(NULL,NULL,gg.box.loo+theme(legend.position='none'),
          NULL,gg.box.loo.TCGA+theme(legend.position='none'),
          NULL,gg.box.loo.ICGC+theme(legend.position='none'),ncol =1,
          labels=c(NA,'Combined',NA,"   TCGA",NA,"   ICGC"), 
          heights = c(0.1,0.12,1,0.12,1,0.12,1)) 
dev.off()




# save #### 
save.image('output/training_model.RData') # 
