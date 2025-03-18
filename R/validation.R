################################################################
# R version 4.2.3 (2023-03-15 ucrt)#############################
# Copyright (C) 2023 The R Foundation for Statistical Computing#
# Platform: x86_64-w64-mingw32/x64 (64-bit) ####################


################################################################
# validation and further analysis ############################## 
################################################################

rm(list=ls()) 
setwd('...')
set.seed(1) #

# library ####
library(reshape2)#melt
library(ggplot2) 
library(ggpubr)#stat_compare_means 
#library(Hmisc) #rcorr
#library(corrplot)#corrplot  
library(ssc)#
library(proxy)#dist 
library(randomForest)
library(caret) #Classification and Regression Training ; confusionMatrix()
library(mltools) #calculate Matthews correlation coefficient, mcc(preds, actual)
library(PRROC) #pr.curve,roc.curve
library(ggpubr) #ggarrange  
library(ROCR) #Visualizing the Performance of Scoring Classifiers;performance(), prediction() 
library(survival)  #coxph
library("survminer") #ggsurvplot 
library(pheatmap) #pheatmap
library(gridExtra) #grid.arrange
library(cowplot) #plot_grid
library(grid) #grid
library(MASS) #polr
library(corrplot)#corrplot  



# SEC1. para,function  #### 
source('parameter_function.R')

# SEC2. data ####  
#--feature
load('data/loc.feature.name.RData')
#--model
load('data/model.self.rf.RData')
#--data
load('data/mat.METABRIC.RData') #2009 86 
load('data/mat.TCGA.RData')  #1003 85
load('data/mat.ICGC.RData')    #401 90
load('data/mat.TransNEO.RData')  #164 89 
load('data/mat.NEWTON.RData') #525 64
load('data/mat.CCLE.RData') #853 84
#--mutation data
load("data/mat.ismut.TCGA.RData")#1026 782
load("data/mat.ismut.ICGC.RData")#479 782


## location #### 
#--
loc.METABRIC=list(hrd=match('is.hrd',colnames(mat.METABRIC)),
                  SBS.v3.3.brca=match('SBS1.brca',colnames(mat.METABRIC)):match('SBS30.brca',colnames(mat.METABRIC)),
                  ID.v3.3=match('ID1',colnames(mat.METABRIC)):match('ID18',colnames(mat.METABRIC)),
                  cna.6score=match(c('cnaBurden','cnaLoad',"MACN",'TDP','TDP_size','Chromoth_state'),colnames(mat.METABRIC)),
                  CX.brca=c(match('CX1.brca',colnames(mat.METABRIC)):match('CX9.brca',colnames(mat.METABRIC))),
                  CN.brca=match('CN1.brca',colnames(mat.METABRIC)):match('CN17.brca',colnames(mat.METABRIC)),
                  scar=match(c('scar.loh','scar.lst','scar.tai','scar.sum'),colnames(mat.METABRIC)),
                  color1=match('color1',colnames(mat.METABRIC)),
                  cohort=match('cohort',colnames(mat.METABRIC)) 
) 
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

loc.TransNEO=list(hrd=match('is.hrd',colnames(mat.TransNEO)),
                  SBS.v3.3.brca=match('SBS1.brca',colnames(mat.TransNEO)):match('SBS30.brca',colnames(mat.TransNEO)),
                  ID.v3.3=match('ID1',colnames(mat.TransNEO)):match('ID18',colnames(mat.TransNEO)),
                  cna.6score=match(c('cnaBurden','cnaLoad',"MACN",'TDP','TDP_size','Chromoth_state'),colnames(mat.TransNEO)),
                  CX.brca=c(match('CX1.brca',colnames(mat.TransNEO)):match('CX9.brca',colnames(mat.TransNEO))),
                  CN.brca=match('CN1.brca',colnames(mat.TransNEO)):match('CN17.brca',colnames(mat.TransNEO)),
                  scar=match(c('scar.loh','scar.lst','scar.tai','scar.sum'),colnames(mat.TransNEO)),
                  color1=match('color1',colnames(mat.TransNEO)),
                  cohort=match('cohort',colnames(mat.TransNEO))
)  

loc.NEWTON=list(SBS.v3.3.brca=match('SBS1.brca',colnames(mat.NEWTON)):match('SBS30.brca',colnames(mat.NEWTON)),
                ID.v3.3=match('ID1',colnames(mat.NEWTON)):match('ID18',colnames(mat.NEWTON)),
                cna.6score=match(c('cnaBurden','cnaLoad',"MACN",'TDP','TDP_size','Chromoth_state'),colnames(mat.NEWTON)),
                CX.brca=c(match('CX1.brca',colnames(mat.NEWTON)):match('CX9.brca',colnames(mat.NEWTON))), 
                cohort=match('cohort',colnames(mat.NEWTON))
)  

loc.CCLE=list(hrd=match('is.hrd',colnames(mat.CCLE)),
              SBS.v3.3.brca=match('SBS1.brca',colnames(mat.CCLE)):match('SBS30.brca',colnames(mat.CCLE)),
              ID.v3.3=match('ID1',colnames(mat.CCLE)):match('ID18',colnames(mat.CCLE)),
              cna.6score=match(c('cnaBurden','cnaLoad',"MACN",'TDP','TDP_size','Chromoth_state'),colnames(mat.CCLE)),
              CX.brca=c(match('CX1.brca',colnames(mat.CCLE)):match('CX9.brca',colnames(mat.CCLE))),
              CN.brca=match('CN1.brca',colnames(mat.CCLE)):match('CN17.brca',colnames(mat.CCLE)),
              scar=match(c('scar.loh','scar.lst','scar.tai','scar.sum'),colnames(mat.CCLE)),
              color1=match('color1',colnames(mat.CCLE)),
              cohort=match('cohort',colnames(mat.CCLE))
)  


## matlog.x=scale(log(x+1)) ####  
matlog.METABRIC=mat.METABRIC 
for(i in unique(as.numeric(unlist(loc.METABRIC[c(2:7)])))){ 
  if(sum(na.omit(matlog.METABRIC[,i]))==0|sum(na.omit(matlog.METABRIC[,i])<0)>0){
    matlog.METABRIC[,i]=NA  #all=0 or some<0
  }else{
    matlog.METABRIC[,i]=scale(log(mat.METABRIC[,i]+1))
  }
} 
#
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

#-- 
matlog.TransNEO=mat.TransNEO 
for(i in unique(as.numeric(unlist(loc.TransNEO[c(2:7)])))){ 
  if(sum(na.omit(matlog.TransNEO[,i]))==0|sum(na.omit(matlog.TransNEO[,i])<0)>0){
    matlog.TransNEO[,i]=NA 
  }else{
    matlog.TransNEO[,i]=scale(log(mat.TransNEO[,i]+1))
  }
} 

#-- 
matlog.NEWTON=mat.NEWTON 
for(i in unique(as.numeric(unlist(loc.NEWTON[c(1:4)])))){ 
  if(sum(na.omit(matlog.NEWTON[,i]))==0|sum(na.omit(matlog.NEWTON[,i])<0)>0){
    matlog.NEWTON[,i]=NA 
  }else{
    matlog.NEWTON[,i]=scale(log(mat.NEWTON[,i]+1))
  }
} 

#-- 
matlog.CCLE=mat.CCLE 
for(i in unique(as.numeric(unlist(loc.CCLE[c(2:7)])))){  
  if(sum(na.omit(matlog.CCLE[,i]))==0|sum(na.omit(matlog.CCLE[,i])<0)>0){
    matlog.CCLE[,i]=NA 
  }else{
    matlog.CCLE[,i]=scale(log(mat.CCLE[,i]+1))
  }
} 


#--matlog1 
matlog1.METABRIC=matlog.METABRIC[,c(1,loc.METABRIC$hrd,loc.METABRIC$SBS.v3.3.brca,loc.METABRIC$ID.v3.3,
                                    loc.METABRIC$cna.6score,loc.METABRIC$CX.brca,
                                    loc.METABRIC$CN.brca,loc.METABRIC$scar[4],loc.METABRIC$color1,loc.METABRIC$cohort)]
matlog1.TCGA=matlog.TCGA[,c(1,loc.TCGA$hrd,loc.TCGA$SBS.v3.3.brca,loc.TCGA$ID.v3.3,
                            loc.TCGA$cna.6score,loc.TCGA$CX.brca,
                            loc.TCGA$CN.brca,loc.TCGA$scar[4],loc.TCGA$color1,loc.TCGA$cohort)]
matlog1.ICGC=matlog.ICGC[,c(1,loc.ICGC$hrd,loc.ICGC$SBS.v3.3.brca,loc.ICGC$ID.v3.3,
                            loc.ICGC$cna.6score,loc.ICGC$CX.brca,
                            loc.ICGC$CN.brca,loc.ICGC$scar[4],loc.ICGC$color1,loc.ICGC$cohort)]

matlog1.TransNEO=matlog.TransNEO[,c(1,loc.TransNEO$hrd,loc.TransNEO$SBS.v3.3.brca,loc.TransNEO$ID.v3.3,
                                    loc.TransNEO$cna.6score,loc.TransNEO$CX.brca,
                                    loc.TransNEO$CN.brca,loc.TransNEO$scar[4],
                                    loc.TransNEO$color1,loc.TransNEO$cohort)] 

matlog1.NEWTON=matlog.NEWTON[,c(1,loc.NEWTON$SBS.v3.3.brca,loc.NEWTON$ID.v3.3,
                                loc.NEWTON$cna.6score,loc.NEWTON$CX.brca, loc.NEWTON$cohort)] 

matlog1.CCLE=matlog.CCLE[,c(1,loc.CCLE$hrd,loc.CCLE$SBS.v3.3.brca,loc.CCLE$ID.v3.3,
                            loc.CCLE$cna.6score,loc.CCLE$CX.brca,
                            loc.CCLE$CN.brca,loc.CCLE$scar[4],
                            loc.CCLE$color1,loc.CCLE$cohort)]



# SEC3. METABRIC ####
## prob/pred ####        
mat.model.METABRIC=list(model1=data.frame(matlog1.METABRIC[,match(c('id','is.hrd',loc.feature.name$model1),colnames(matlog1.METABRIC))]),
                        model2=data.frame(matlog1.METABRIC[,match(c('id','is.hrd',loc.feature.name$model2),colnames(matlog1.METABRIC))]),
                        model3=data.frame(matlog1.METABRIC[,match(c('id','is.hrd',loc.feature.name$model3),colnames(matlog1.METABRIC))]),
                        model4=data.frame(matlog1.METABRIC[,match(c('id','is.hrd',loc.feature.name$model4),colnames(matlog1.METABRIC))]),
                        model5=data.frame(matlog1.METABRIC[,match(c('id','is.hrd',loc.feature.name$model5),colnames(matlog1.METABRIC))]))

for(m in 1:n.model){
  mat.model.METABRIC[[m]][,2]=factor(mat.model.METABRIC[[m]][,2],levels = c(1,0))
}

#--
pred.METABRIC=prob.METABRIC=list() 
for(m in 1:n.model){
  print(m) 
  prob.METABRIC[[m]]=ssc:::predProb(model.self.rf[[m]]$model,data.matrix(mat.model.METABRIC[[m]][,-c(1,2)]),
                                    model.self.rf[[m]]$pred,model.self.rf[[m]]$pred.pars)[,1]
  pred.METABRIC[[m]]=ifelse(prob.METABRIC[[m]]>=mat.cut.f1$cut[m],1,0) 
} 

#--add
prob.METABRIC.mat=cbind(prob.METABRIC[[1]],prob.METABRIC[[2]],prob.METABRIC[[3]],
                        prob.METABRIC[[4]],prob.METABRIC[[5]])
prob.METABRIC.mat=data.frame(prob.METABRIC.mat)
colnames(prob.METABRIC.mat)=paste0('prob.',1:n.model)
#
pred.METABRIC.mat=cbind(pred.METABRIC[[1]],pred.METABRIC[[2]],pred.METABRIC[[3]],
                        pred.METABRIC[[4]],pred.METABRIC[[5]])
pred.METABRIC.mat=data.frame(pred.METABRIC.mat)
colnames(pred.METABRIC.mat)=paste0('pred.',1:n.model)
#
mat.METABRIC.add=cbind(mat.METABRIC,prob.METABRIC.mat,pred.METABRIC.mat) 


## RFS, relapse free survival #### 
loc.na.RFS=unique(c(which(is.na(mat.METABRIC.add$RFS.time)),which(mat.METABRIC.add$RFS.time<1),which(is.na(mat.METABRIC.add$RFS.event))))#186 1216
METABRIC.SA=data.frame(id=mat.METABRIC.add$id[-loc.na.RFS],
                       time=mat.METABRIC.add$RFS.time[-loc.na.RFS],
                       event=mat.METABRIC.add$RFS.event[-loc.na.RFS],
                       mat.METABRIC.add[-loc.na.RFS,])


### hazard ratio in CT+&ER- ####
loc.formula=grep('prob.',colnames(mat.METABRIC.add)) 
formula.test=sapply(colnames(mat.METABRIC.add)[loc.formula],
                    function(x) as.formula(paste('Surv(time, event)~', x)))
#--coxph(Surv())  
#--ER-
METABRIC.SA.CTER.result=matrix(,length(formula.test),10) 
for(i in 1:length(formula.test)){METABRIC.SA.CTER.result[i,]= unlist(c(hazard_table(coxph(formula.test[[i]], 
                                                                                          data = METABRIC.SA[which(METABRIC.SA$CT=='CT+'&METABRIC.SA$ER=='ER-'),])),'RFS'))}

#
METABRIC.SA.CTER.result=data.frame(METABRIC.SA.CTER.result)
colnames(METABRIC.SA.CTER.result)=c('feature','HR','Lower','Higher','P','formulas','n','nevent',
                                    'p.sig','type')  
#--chr to num
METABRIC.SA.CTER.result$cohort='METABRIC'
METABRIC.SA.CTER.result$HR=as.numeric(METABRIC.SA.CTER.result$HR)
METABRIC.SA.CTER.result$Lower=as.numeric(METABRIC.SA.CTER.result$Lower)
METABRIC.SA.CTER.result$Higher=as.numeric(METABRIC.SA.CTER.result$Higher)
METABRIC.SA.CTER.result$P=as.numeric(METABRIC.SA.CTER.result$P)
METABRIC.SA.CTER.result$index=c(1:dim(METABRIC.SA.CTER.result)[1])
METABRIC.SA.CTER.result$type1=name.model 
METABRIC.SA.CTER.result$ER='ER-' 


#--ER+
METABRIC.SA.CTERp.result=matrix(,length(formula.test),10) 
for(i in 1:length(formula.test)){METABRIC.SA.CTERp.result[i,]= unlist(c(hazard_table(coxph(formula.test[[i]], data = METABRIC.SA[which(METABRIC.SA$CT=='CT+'&METABRIC.SA$ER=='ER+'),])),'RFS'))}

#
METABRIC.SA.CTERp.result=data.frame(METABRIC.SA.CTERp.result)
colnames(METABRIC.SA.CTERp.result)=c('feature','HR','Lower','Higher','P','formulas','n','nevent',
                                     'p.sig','type')  
METABRIC.SA.CTERp.result$cohort='METABRIC'
METABRIC.SA.CTERp.result$HR=as.numeric(METABRIC.SA.CTERp.result$HR)
METABRIC.SA.CTERp.result$Lower=as.numeric(METABRIC.SA.CTERp.result$Lower)
METABRIC.SA.CTERp.result$Higher=as.numeric(METABRIC.SA.CTERp.result$Higher)
METABRIC.SA.CTERp.result$P=as.numeric(METABRIC.SA.CTERp.result$P)
METABRIC.SA.CTERp.result$index=c(1:dim(METABRIC.SA.CTERp.result)[1])
METABRIC.SA.CTERp.result$type1=name.model 
METABRIC.SA.CTERp.result$ER='ER+' 


#--fig
mat.HR.CTER.METABRIC.RFS=rbind(METABRIC.SA.CTER.result,METABRIC.SA.CTERp.result)
#
gg.HR.CT_ER.METABRIC.RFS.prob=ggplot(data=mat.HR.CTER.METABRIC.RFS, 
                                     aes(y=index, x=HR, xmin=Lower, xmax=Higher)) +  
  geom_point(aes(color=p.sig),size=2)+ 
  geom_errorbarh(aes(xmax =Higher, xmin =Lower,color=p.sig), size = .5, height = .2) + 
  scale_y_continuous(breaks=1:dim(mat.HR.CTER.METABRIC.RFS)[1], 
                     labels=mat.HR.CTER.METABRIC.RFS$type1, 
                     trans = "reverse") + 
  labs( x='Hazard ratio (95% CI)', y = '') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  scale_color_manual(values=c("#375E97","gray70","#FB6542"),breaks=c("good","ns","poor"))+
  coord_cartesian(xlim=c(0,3))+
  facet_wrap(~ER, ncol=2)+
  theme(legend.position = "none",
        axis.line = element_blank(),
        strip.text = element_text(size = 15),
        strip.background =element_blank(), 
        plot.margin = unit(c(1,2,1,0.5), "lines"),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15)) 
#--
pdf('fig_SA_HR_METABRIC_CT_ER.pdf')#
gg.HR.CT_ER.METABRIC.RFS.prob
dev.off() 


# SEC4. TCGA ####
## prob/pred ####
mat.model.TCGA=list(model1=data.frame(matlog1.TCGA[,match(c('id','is.hrd',loc.feature.name$model1),colnames(matlog1.TCGA))]),
                    model2=data.frame(matlog1.TCGA[,match(c('id','is.hrd',loc.feature.name$model2),colnames(matlog1.TCGA))]),
                    model3=data.frame(matlog1.TCGA[,match(c('id','is.hrd',loc.feature.name$model3),colnames(matlog1.TCGA))]),
                    model4=data.frame(matlog1.TCGA[,match(c('id','is.hrd',loc.feature.name$model4),colnames(matlog1.TCGA))]),
                    model5=data.frame(matlog1.TCGA[,match(c('id','is.hrd',loc.feature.name$model5),colnames(matlog1.TCGA))]))

for(m in 1:n.model){
  mat.model.TCGA[[m]][,2]=factor(mat.model.TCGA[[m]][,2],levels = c(1,0))
}

#--
pred.TCGA=prob.TCGA=list() 
for(m in 1:n.model){ 
  prob.TCGA[[m]]=ssc:::predProb(model.self.rf[[m]]$model,data.matrix(mat.model.TCGA[[m]][,-c(1,2)]),
                                model.self.rf[[m]]$pred,model.self.rf[[m]]$pred.pars)[,1]
  pred.TCGA[[m]]=ifelse(prob.TCGA[[m]]>=mat.cut.f1$cut[m],1,0) 
} 

#--add  
prob.TCGA.mat=cbind(prob.TCGA[[1]],prob.TCGA[[2]],prob.TCGA[[3]],
                    prob.TCGA[[4]],prob.TCGA[[5]])
prob.TCGA.mat=data.frame(prob.TCGA.mat)
colnames(prob.TCGA.mat)=paste0('prob.',1:n.model)
#
pred.TCGA.mat=cbind(pred.TCGA[[1]],pred.TCGA[[2]],pred.TCGA[[3]],
                    pred.TCGA[[4]],pred.TCGA[[5]])
pred.TCGA.mat=data.frame(pred.TCGA.mat)
colnames(pred.TCGA.mat)=paste0('pred.',1:n.model)
#
mat.TCGA.add=cbind(mat.TCGA,prob.TCGA.mat,pred.TCGA.mat) 


# SEC5. ICGC ####
## prob/pred ####
mat.model.ICGC=list(model1=data.frame(matlog1.ICGC[,match(c('id','is.hrd',loc.feature.name$model1),colnames(matlog1.ICGC))]),
                    model2=data.frame(matlog1.ICGC[,match(c('id','is.hrd',loc.feature.name$model2),colnames(matlog1.ICGC))]),
                    model3=data.frame(matlog1.ICGC[,match(c('id','is.hrd',loc.feature.name$model3),colnames(matlog1.ICGC))]),
                    model4=data.frame(matlog1.ICGC[,match(c('id','is.hrd',loc.feature.name$model4),colnames(matlog1.ICGC))]),
                    model5=data.frame(matlog1.ICGC[,match(c('id','is.hrd',loc.feature.name$model5),colnames(matlog1.ICGC))]))

for(m in 1:n.model){
  mat.model.ICGC[[m]][,2]=factor(mat.model.ICGC[[m]][,2],levels = c(1,0))
}

#--
pred.ICGC=prob.ICGC=list() 
for(m in 1:n.model){ 
  prob.ICGC[[m]]=ssc:::predProb(model.self.rf[[m]]$model,data.matrix(mat.model.ICGC[[m]][,-c(1,2)]),
                                model.self.rf[[m]]$pred,model.self.rf[[m]]$pred.pars)[,1]
  pred.ICGC[[m]]=ifelse(prob.ICGC[[m]]>=mat.cut.f1$cut[m],1,0)
}

#-- add 
prob.ICGC.mat=cbind(prob.ICGC[[1]],prob.ICGC[[2]],prob.ICGC[[3]],
                    prob.ICGC[[4]],prob.ICGC[[5]])
prob.ICGC.mat=data.frame(prob.ICGC.mat)
colnames(prob.ICGC.mat)=paste0('prob.',1:n.model)
#
pred.ICGC.mat=cbind(pred.ICGC[[1]],pred.ICGC[[2]],pred.ICGC[[3]],
                    pred.ICGC[[4]],pred.ICGC[[5]])
pred.ICGC.mat=data.frame(pred.ICGC.mat)
colnames(pred.ICGC.mat)=paste0('pred.',1:n.model)
#
mat.ICGC.add=cbind(mat.ICGC,prob.ICGC.mat,pred.ICGC.mat) 



# SEC6. mutation analysis #### 
#781 genes from CHORD
load('data/gene781.RData') 


## fisher.test ####
loc.pred=grep('pred.',colnames(mat.TCGA.add)) 
mat.mut.TCGA.pred=merge(mat.mut.TCGA,mat.TCGA.add[,c(1,loc.pred)],by='id') 
mat.mut.TCGA.pred=mat.mut.TCGA.pred[match(mat.TCGA.add$id,mat.mut.TCGA.pred$id),]
#--ICGC
loc.pred=grep('pred.',colnames(mat.ICGC.add)) 
mat.mut.ICGC.pred=merge(mat.mut.ICGC,mat.ICGC.add[,c(1,loc.pred)],by='id') 
mat.mut.ICGC.pred=mat.mut.ICGC.pred[match(mat.ICGC.add$id,mat.mut.ICGC.pred$id),]
#--
mat.mut.2cohort.pred=rbind(mat.mut.TCGA.pred,mat.mut.ICGC.pred)#1404 787
dim(mat.mut.2cohort.pred)

fisher.2cohort=list()
fisher.2cohort.pval=fisher.2cohort.OR=matrix(,length(gene781),n.model)
for(g in 1:length(gene781)){ 
  print(g)
  fisher.2cohort[[g]]=list()
  if(sum(mat.mut.2cohort.pred[,g+1])==0){ #this gth gene no mut
    fisher.2cohort[[g]]=NA
    fisher.2cohort.pval[g,]=fisher.2cohort.OR[g,]=NA
  }else{
    for(m in 1:n.model){ 
      fisher.2cohort[[g]][[m]]=fisher.test(table(mat.mut.2cohort.pred[,g+1],
                                                 mat.mut.2cohort.pred[,782+m]))
      fisher.2cohort.pval[g,m]=fisher.2cohort[[g]][[m]]$p.value 
      fisher.2cohort.OR[g,m]=fisher.2cohort[[g]][[m]]$estimate
    }
  }
} 
#--p value
fisher.2cohort.pval=data.frame(fisher.2cohort.pval)
rownames(fisher.2cohort.pval)=gene781
colnames(fisher.2cohort.pval)=name.model 

#--keep significant cases
fisher.2cohort.pval.sig=fisher.2cohort.pval
for(m in 1:n.model){
  fisher.2cohort.pval.sig[,m]=ifelse(fisher.2cohort.pval[,m]<=0.05,1,0)
}
fisher.2cohort.pval.sig.v1=fisher.2cohort.pval.sig[which(rowSums(fisher.2cohort.pval.sig)!=0),]
fisher.2cohort.pval.sig.v1=data.frame(fisher.2cohort.pval.sig.v1)

#--p.adjust
fisher.2cohort.adj=fisher.2cohort.pval
for(m in 1:n.model){ 
  fisher.2cohort.adj[,m]=p.adjust(fisher.2cohort.pval[,m],'hochberg')
}
fisher.2cohort.adj.sig=ifelse(fisher.2cohort.adj<=0.05,1,0) #1=sig 
fisher.2cohort.adj.sig.v1=fisher.2cohort.adj.sig[which(rowSums(fisher.2cohort.adj.sig)!=0),]
fisher.2cohort.adj.sig.v1=data.frame(fisher.2cohort.adj.sig.v1)  

#--OR
fisher.2cohort.OR=data.frame(fisher.2cohort.OR)
rownames(fisher.2cohort.OR)=gene781
colnames(fisher.2cohort.OR)=name.model 

#--OR with  significant cases
fisher.2cohort.OR.v1=fisher.2cohort.OR[rownames(fisher.2cohort.OR)%in%rownames(fisher.2cohort.adj.sig.v1),] 

#--change order with OR
rank(fisher.2cohort.OR.v1$`SNV+ASCN`)
fisher.2cohort.OR.v2=fisher.2cohort.OR.v1[order(fisher.2cohort.OR.v1$`SNV+ASCN`,decreasing =T),]
fisher.2cohort.adj.sig.v2=fisher.2cohort.adj.sig.v1[order(fisher.2cohort.OR.v1$`SNV+ASCN`,decreasing =T),]

#--genes with OR>1 
fisher.2cohort.sig.gene=rownames(fisher.2cohort.OR.v2)
fisher.2cohort.sig.gene.OR1=rownames(fisher.2cohort.OR.v2)[1:7]
#"DPYD"   "BRCA2"  "BRCA1"  "SETDB1" "DCC"    "EPHA3"  "TP53" 
fisher.2cohort.sig.gene.OR0=rownames(fisher.2cohort.OR.v2)[8:11]
#"GATA3"   "PIK3CA" "MAP3K1" "CDH1" 




## mutational type vs HRD prob #### 
variant.all.sort=c("Splice_Site","Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation",
                   "Missense_Mutation","In_Frame_Del","In_Frame_Ins","Nonstop_Mutation",
                   "Translation_Start_Site","Readthrough_Mutation",
                   "RNA","DNV" )
mat.muttype.2cohort=rbind(mat.muttype.TCGA,mat.muttype.ICGC)
loc.gene.somatic=match(fisher.2cohort.sig.gene.OR1,colnames(mat.muttype.2cohort)) 
df.muttype.2cohort=list()
for(i in 1:length(loc.gene.somatic)){
  print(i)
  tmp=data.frame(gene=mat.muttype.2cohort[,loc.gene.somatic[i]], 
                 mat.muttype.2cohort$prob.1,
                 mat.muttype.2cohort$prob.2,
                 mat.muttype.2cohort$prob.3,
                 mat.muttype.2cohort$prob.4,
                 mat.muttype.2cohort$prob.5)
  colnames(tmp)=c('gene',name.model)
  tmp=na.omit(tmp)  
  tmp=melt(tmp) 
  tmp$gene=factor(tmp$gene,levels = unique(tmp$gene)[na.omit(match(variant.all.sort,unique(tmp$gene)))])
  colnames(tmp)=c('gene','model','prob') 
  df.muttype.2cohort[[i]]=tmp
} 
names(df.muttype.2cohort)=fisher.2cohort.sig.gene.OR1


## OR heatmap ####
#--set unsig as -1 in white
my_palette <- c('white',colorRampPalette(colors = c("darkblue", "lightblue"))(n =length(seq(0,0.9,by=0.1))),
                "gray38", "gray38",
                c(colorRampPalette(colors = c("tomato1","darkred"))(n =length(seq(1.1,25,by=1)))))

#--
tmp=which(fisher.2cohort.adj.sig.v2==0,arr.ind = TRUE) 
tmp1=fisher.2cohort.OR.v2 
if(dim(tmp)[1]!=0){
  for(i in 1:dim(tmp)[1]){
    tmp1[tmp[i,1],tmp[i,2]]=-1 #set unsig as -1 in white
  }
}
tmp1=data.frame(tmp1) 
tmp1$gene=rownames(tmp1) 
df.heatmap.2cohort=tmp1 
colnames(df.heatmap.2cohort)[1:5]=name.model

#--fig
newnames <- lapply(rownames(df.heatmap.2cohort),function(x) bquote(italic(.(x))))
heatmap.2cohort=pheatmap(df.heatmap.2cohort[,1:5],cluster_row = FALSE, cluster_cols= FALSE, 
                         color = my_palette, breaks =c(-1,seq(0,0.9,by=0.1),0.9999,1.0001,seq(1.1,25,by=1)), 
                         scale = "none", #legend=FALSE,
                         main ='combined TCGA and ICGC cohorts',
                         labels_row = as.expression(newnames))


#--
gg.muttype.2cohort=gg.muttype.2cohort.v2=gg.muttype.2cohort.n=list()
for(i in 1:length(loc.gene.somatic)){   
  gg.muttype.2cohort[[i]]=ggplot(df.muttype.2cohort[[i]], aes(gene, prob)) + geom_boxplot() +ylim(c(0,1))+
    facet_grid(~model)+labs(x=fisher.2cohort.sig.gene.OR1[i])+
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
    coord_flip()+ylab('Predicted probabilities')+ xlab(" ")
  
  ## Create the table-base pallete
  table_base <- ggplot(data.frame(table(df.muttype.2cohort[[i]]$gene)/n.model), aes(y=Var1)) +
    ylab(NULL) + xlab("  ") + 
    theme(plot.title = element_text(hjust = 0.5, size=10), 
          axis.text.x = element_text(color="white", hjust =-3, size = 3), # This is used to help with alignment
          axis.line = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.y = element_blank(), 
          legend.position = "none",
          panel.background = element_blank(), 
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()  
    )
  #table_base
  ## HR table
  tab1 <- table_base + 
    labs(title = "space") +
    geom_text(aes(y =Var1, x = 1, label = Freq), size =3) +  
    ggtitle("n")
  gg.muttype.2cohort.n[[i]]=tab1
  #
  lay <-  matrix(c(rep(1,12),2), nrow = 1)
  gg.muttype.2cohort.v2[[i]]=ggarrange(grid.arrange(gg.muttype.2cohort[[i]], tab1,layout_matrix = lay))
}

#--
tmp=ggarrange(gg.muttype.2cohort.v2[[1]]+xlab(''),#+ylab('Predicted probabilities'),
              gg.muttype.2cohort.v2[[2]]+xlab(''),#+ylab('Predicted probabilities'),
              gg.muttype.2cohort.v2[[3]]+xlab(''),#+ylab('Predicted probabilities'),
              gg.muttype.2cohort.v2[[4]]+xlab(''),#+ylab('Predicted probabilities'),
              gg.muttype.2cohort.v2[[5]]+xlab(''),#+ylab('Predicted probabilities'), 
              gg.muttype.2cohort.v2[[6]]+xlab(''),
              gg.muttype.2cohort.v2[[7]]+xlab(''),
              labels =fisher.2cohort.sig.gene.OR1,ncol=4,nrow =2,
              font.label = list(face = "bold.italic"))


# SEC7. TransNEO ####
## prob/pred ####
mat.model.TransNEO=list(model1=data.frame(matlog1.TransNEO[,match(c('id','is.hrd',loc.feature.name$model1),colnames(matlog1.TransNEO))]),
                        model2=data.frame(matlog1.TransNEO[,match(c('id','is.hrd',loc.feature.name$model2),colnames(matlog1.TransNEO))]),
                        model3=data.frame(matlog1.TransNEO[,match(c('id','is.hrd',loc.feature.name$model3),colnames(matlog1.TransNEO))]),
                        model4=data.frame(matlog1.TransNEO[,match(c('id','is.hrd',loc.feature.name$model4),colnames(matlog1.TransNEO))]),
                        model5=data.frame(matlog1.TransNEO[,match(c('id','is.hrd',loc.feature.name$model5),colnames(matlog1.TransNEO))]))

for(m in 1:n.model){
  mat.model.TransNEO[[m]][,2]=factor(mat.model.TransNEO[[m]][,2],levels = c(1,0))
}

#--
pred.TransNEO=prob.TransNEO=list() 
for(m in 1:n.model){ 
  prob.TransNEO[[m]]=ssc:::predProb(model.self.rf[[m]]$model,data.matrix(mat.model.TransNEO[[m]][,-c(1,2)]),
                                    model.self.rf[[m]]$pred,model.self.rf[[m]]$pred.pars)[,1]
  pred.TransNEO[[m]]=ifelse(prob.TransNEO[[m]]>=mat.cut.f1$cut[m],1,0)
}

#-- add pred/prob 
prob.TransNEO.mat=cbind(prob.TransNEO[[1]],prob.TransNEO[[2]],prob.TransNEO[[3]],
                        prob.TransNEO[[4]],prob.TransNEO[[5]])
prob.TransNEO.mat=data.frame(prob.TransNEO.mat)
colnames(prob.TransNEO.mat)=paste0('prob.',1:n.model)
#
pred.TransNEO.mat=cbind(pred.TransNEO[[1]],pred.TransNEO[[2]],pred.TransNEO[[3]],
                        pred.TransNEO[[4]],pred.TransNEO[[5]])
pred.TransNEO.mat=data.frame(pred.TransNEO.mat)
colnames(pred.TransNEO.mat)=paste0('pred.',1:n.model)
#
mat.TransNEO.add=cbind(mat.TransNEO,prob.TransNEO.mat,pred.TransNEO.mat) 


### boxplot HRD #### 
mat.prob.TransNEO=cbind(prob.TransNEO[[1]],prob.TransNEO[[2]],prob.TransNEO[[3]],prob.TransNEO[[4]],prob.TransNEO[[5]]) 
mat.prob.TransNEO=data.frame(mat.prob.TransNEO)
colnames(mat.prob.TransNEO)=name.model 
mat.prob.TransNEO$is.hrd=matlog1.TransNEO$is.hrd
mat.prob.TransNEO=mat.prob.TransNEO[!is.na(mat.prob.TransNEO$is.hrd),]
mat.prob.TransNEO$is.hrd=factor(mat.prob.TransNEO$is.hrd,levels = c(1,0),labels= c('HRD+','HRD-')) 
mat.prob.TransNEO=melt(mat.prob.TransNEO) 
colnames(mat.prob.TransNEO)[c(1,3)]=c('HRD','prob') 
#
gg.box.TransNEO=ggplot(mat.prob.TransNEO, aes(x=HRD, y=prob, fill=HRD)) +  
  geom_boxplot()+stat_compare_means(method = "wilcox.test",label = "p.signif",color="blue")+ 
  theme_minimal()+
  xlab('')+ ylab('Predicted probabilities')+  
  facet_wrap(~variable,ncol=5)+ 
  geom_hline(data = mat.cut.f1, aes(yintercept = cut), linetype="dashed", color = "grey")  
#
pdf(file ="fig_box_TransNEO.pdf",width=8, height=4)
gg.box.TransNEO+theme(legend.position='none')
dev.off()


### boxplot RCB #### 
mat.prob.TransNEO.RCB=cbind(prob.TransNEO[[1]],prob.TransNEO[[2]],prob.TransNEO[[3]],prob.TransNEO[[4]],prob.TransNEO[[5]]) 
mat.prob.TransNEO.RCB=data.frame(mat.prob.TransNEO.RCB)
colnames(mat.prob.TransNEO.RCB)=name.model  
mat.prob.TransNEO.RCB$RCB.category=mat.TransNEO.add$RCB.category
mat.prob.TransNEO.RCB=mat.prob.TransNEO.RCB[!is.na(mat.prob.TransNEO.RCB$RCB.category),] 
mat.prob.TransNEO.RCB=melt(mat.prob.TransNEO.RCB)

gg.TransNEO.box.RCB=ggplot(mat.prob.TransNEO.RCB,aes(x=RCB.category,y=value,fill=RCB.category))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  facet_wrap(~variable,ncol=5)+ geom_hline(yintercept=0.5, linetype="dashed", color = "grey")+  
  ylab('Predicted probabilities')+
  stat_compare_means(ref.group="pCR",label = "p.signif",color="#FB6542",
                     label.y.npc = 0.8,hide.ns = T)+ 
  theme(panel.grid.major.x = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,1), "lines"))+guides(fill="none")
#
pdf('fig_TransNEO_box_RCB.pdf') 
gg.TransNEO.box.RCB 
dev.off()


## monotonically associated with RCB,polr() #### 
loc.prob=grep('prob.',colnames(mat.TransNEO.add))  
formula.TransNEO=colnames(mat.TransNEO.add)[loc.prob] 

#Ordered Logistic or Probit Regression, check valid pred is monotonically associated with RCB class 
TransNEO.polr.RCB=list()
TransNEO.polr.RCB.pval=c() 
for(i in 1:length(formula.TransNEO)){
  TransNEO.polr.RCB[[i]]=polr(as.formula(paste('RCB.category~',formula.TransNEO[i])), 
                              mat.TransNEO.add[!is.na(mat.TransNEO.add$RCB.category),],  Hess=TRUE)
  TransNEO.polr.RCB.pval=c(TransNEO.polr.RCB.pval,pnorm(abs(coef(summary(TransNEO.polr.RCB[[i]]))[1, "t value"]), lower.tail = FALSE) * 2)
}
TransNEO.polr.RCB.pval= round(TransNEO.polr.RCB.pval,3)
names(TransNEO.polr.RCB.pval)=formula.TransNEO 


## OR, Associations with pCR ####
#--simple glm 
TransNEO.glm.pCR=list()
TransNEO.glm.pCR.pval=c() 
TransNEO.glm.pCR[[1]]=glm(as.formula(paste('pCR.RD~',formula.TransNEO[1])), 
                          mat.TransNEO.add[!is.na(mat.TransNEO.add$pCR.RD),],family = "binomial")
TransNEO.glm.pCR.pval=c(TransNEO.glm.pCR.pval,summary(TransNEO.glm.pCR[[1]])$coefficients[2,'Pr(>|z|)'])

for(i in 2:length(formula.TransNEO)){
  TransNEO.glm.pCR[[i]]=glm(as.formula(paste('pCR.RD~',formula.TransNEO[i])), 
                            mat.TransNEO.add[!is.na(mat.TransNEO.add$pCR.RD),],family = "binomial")
  TransNEO.glm.pCR.pval=c(TransNEO.glm.pCR.pval,summary(TransNEO.glm.pCR[[i]])$coefficients[2,'Pr(>|z|)'])
}
TransNEO.glm.pCR.pval= round(TransNEO.glm.pCR.pval,3)
names(TransNEO.glm.pCR.pval)=formula.TransNEO  

#--odds.ratio
TransNEO.OR.pCR = odds.ratio_table(TransNEO.glm.pCR[[1]])  
TransNEO.OR.pCR$method=formula.TransNEO[1]
for(i in 2:length(formula.TransNEO)){
  if(!is.na(TransNEO.glm.pCR)[i]){
    tmp=odds.ratio_table(TransNEO.glm.pCR[[i]])
    tmp$method=formula.TransNEO[i]
    TransNEO.OR.pCR=rbind(TransNEO.OR.pCR,tmp)
  }
} 
TransNEO.OR.pCR$class <- "pCR"
TransNEO.OR.pCR$type <- "Simple logist." 
TransNEO.OR.pCR$index=c(1:dim(TransNEO.OR.pCR)[1]) 
TransNEO.OR.pCR$feature=name.model
#
gg.TransNEO.OR.pCR.5prob=ggplot(data=TransNEO.OR.pCR, aes(y=index, x=odds.ratio, xmin=ci.low, xmax=ci.hi)) + #
  geom_point(aes(color=p.sig),size=2)+
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dotted") + 
  geom_errorbarh(aes(xmax =ci.hi, xmin =ci.low,color=p.sig), height = .2) +  
  scale_y_continuous(breaks=1:dim(TransNEO.OR.pCR)[1], 
                     labels=TransNEO.OR.pCR$feature, 
                     trans = "reverse") + 
  labs(  x='Odds ratio (95% CI)', y = '') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  scale_color_manual(values=c("#375E97","gray70","#FB6542"),breaks=c("down","ns","up"))+
  coord_cartesian(xlim=c(0,2))+
  theme(legend.position = "none",
        axis.line = element_blank(),
        strip.text = element_text(size = 10),
        strip.background =element_blank(), 
        axis.text=element_text(size=10),
        axis.title=element_text(size=10)) 
#
pdf('fig_TransNEO_box_RCB_OR.pdf', width =13, height =4) #
ggarrange(gg.TransNEO.box.RCB,gg.TransNEO.OR.pCR.5prob, 
          nrow =1,ncol=2,widths =c(3.7,1),labels = c('a)','b)')) 
dev.off()



# SCAN-B, see validation_SCANB.R ####



# SEC8. NEWTON,no mcn,model1,3,4;CT ####  
matlog1.NEWTON$CX3.brca=0
mat.model.NEWTON=list(model1=data.frame(matlog1.NEWTON[,c(1,match(loc.feature.name$model1,colnames(matlog1.NEWTON)))]),
                      model2=NA,
                      model3=data.frame(matlog1.NEWTON[,c(1,match(loc.feature.name$model3,colnames(matlog1.NEWTON)))]), 
                      model4=data.frame(matlog1.NEWTON[,c(1,match(loc.feature.name$model4,colnames(matlog1.NEWTON)))]),
                      model5=NA) 
#--
pred.NEWTON=prob.NEWTON=list() 
for(m in 1:n.model){
  print(m)
  if(m==2|m==5){
    pred.NEWTON[[m]]=NA
    prob.NEWTON[[m]]=NA
  }else{  
    prob.NEWTON[[m]]=ssc:::predProb(model.self.rf[[m]]$model,data.matrix(mat.model.NEWTON[[m]][,-1]),
                                    model.self.rf[[m]]$pred,model.self.rf[[m]]$pred.pars)[,1]
    pred.NEWTON[[m]]=factor(ifelse(prob.NEWTON[[m]]>=mat.cut.f1$cut[m],1,0),levels=c(1,0)) 
  }
}  


#-- add pred/prob 
mat.NEWTON.add=mat.NEWTON
mat.NEWTON.add$prob.1=prob.NEWTON[[1]]
mat.NEWTON.add$prob.3=prob.NEWTON[[3]]
mat.NEWTON.add$prob.4=prob.NEWTON[[4]]
mat.NEWTON.add$pred.1=pred.NEWTON[[1]]
mat.NEWTON.add$pred.3=pred.NEWTON[[3]]
mat.NEWTON.add$pred.4=pred.NEWTON[[4]]




# SEC9. CCLE ####
## prob/pred ####
mat.model.CCLE=list(model1=data.frame(matlog1.CCLE[,match(c('id','is.hrd',loc.feature.name$model1),colnames(matlog1.CCLE))]),
                    model2=data.frame(matlog1.CCLE[,match(c('id','is.hrd',loc.feature.name$model2),colnames(matlog1.CCLE))]),
                    model3=data.frame(matlog1.CCLE[,match(c('id','is.hrd',loc.feature.name$model3),colnames(matlog1.CCLE))]),
                    model4=data.frame(matlog1.CCLE[,match(c('id','is.hrd',loc.feature.name$model4),colnames(matlog1.CCLE))]),
                    model5=data.frame(matlog1.CCLE[,match(c('id','is.hrd',loc.feature.name$model5),colnames(matlog1.CCLE))]))

for(m in 1:n.model){
  mat.model.CCLE[[m]][,2]=factor(mat.model.CCLE[[m]][,2],levels = c(1,0))
}

#--
pred.CCLE=prob.CCLE=list() 
for(m in 1:n.model){ 
  prob.CCLE[[m]]=ssc:::predProb(model.self.rf[[m]]$model,data.matrix(mat.model.CCLE[[m]][,-c(1,2)]),
                                model.self.rf[[m]]$pred,model.self.rf[[m]]$pred.pars)[,1]
  pred.CCLE[[m]]=ifelse(prob.CCLE[[m]]>=mat.cut.f1$cut[m],1,0)
}
#-- add 
prob.CCLE.mat=cbind(prob.CCLE[[1]],prob.CCLE[[2]],prob.CCLE[[3]],
                    prob.CCLE[[4]],prob.CCLE[[5]])
prob.CCLE.mat=data.frame(prob.CCLE.mat)
colnames(prob.CCLE.mat)=paste0('prob.',1:n.model)
#
pred.CCLE.mat=cbind(pred.CCLE[[1]],pred.CCLE[[2]],pred.CCLE[[3]],
                    pred.CCLE[[4]],pred.CCLE[[5]])
pred.CCLE.mat=data.frame(pred.CCLE.mat)
colnames(pred.CCLE.mat)=paste0('pred.',1:n.model)
#
mat.CCLE.add=cbind(mat.CCLE,prob.CCLE.mat,pred.CCLE.mat) 



## corplot with PRISM #### 
#--Correlation between the HRD and sensitivity to PARP inhibitors 
#evaluated using the PRISM metric.
#--only brca
gg.CCLE.cor.PRISM_v1_brca=list()
tmp=mat.CCLE.add[which(mat.CCLE.add$CANCER_TYPE=='Breast Cancer'),c(paste0(c("OLAPARIB","TALAZOPARIB" ,"RUCAPARIB","NIRAPARIB"),'.PRISM'),
                                                                    paste0('prob.',1:n.model))]
name.drug=c("OLAPARIB","TALAZOPARIB" ,"RUCAPARIB","NIRAPARIB")
colnames(tmp)[1:4]=name.drug
colnames(tmp)[5:9]=name.model  
tmp=melt(tmp,id=1:4)  
colnames(tmp)[5:6]=c('model','prob') #
for(m in 1:4){
  tmp1=tmp[,c(m,5,6)]
  colnames(tmp1)[1]='PRISM'
  gg.CCLE.cor.PRISM_v1_brca[[m]]=ggscatter(tmp1, x = "prob", y = "PRISM",size = 2,
                                           add = "reg.line",  
                                           add.params = list(color = "blue", fill = "lightgray"), 
                                           conf.int = TRUE # Add confidence interval
  )+facet_wrap(~model,ncol=n.model)+ylab(name.drug[m])+xlab('Predicted probabilities')+ 
    stat_cor(method = "pearson", label.x =0.1, label.y =0.4,size=5)
} 
#
pdf(file ="fig_CCLE_cor_PRISM_v1_brca.pdf", width = 12, height =10)
ggarrange(gg.CCLE.cor.PRISM_v1_brca[[1]]+xlab(''),gg.CCLE.cor.PRISM_v1_brca[[2]]+xlab(''),
          gg.CCLE.cor.PRISM_v1_brca[[3]]+xlab(''),
          gg.CCLE.cor.PRISM_v1_brca[[4]],
          ncol=1,nrow=4) 
dev.off()



# SEC10. clinical #### 
## chisq.test #### 
mat.clinical.test.pred=data.frame(type=name.clinical)
mat.clinical.test.pred$METABRIC.pred.1=round(c(chisq.test(table(mat.METABRIC.add[,c('pred.1',name.clinical[1])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.1',name.clinical[2])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.1',name.clinical[3])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.1',name.clinical[4])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.1',name.clinical[5])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.1',name.clinical[6])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.1',name.clinical[7])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.1',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$METABRIC.pred.2=round(c(chisq.test(table(mat.METABRIC.add[,c('pred.2',name.clinical[1])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.2',name.clinical[2])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.2',name.clinical[3])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.2',name.clinical[4])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.2',name.clinical[5])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.2',name.clinical[6])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.2',name.clinical[7])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.2',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$METABRIC.pred.3=round(c(chisq.test(table(mat.METABRIC.add[,c('pred.3',name.clinical[1])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.3',name.clinical[2])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.3',name.clinical[3])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.3',name.clinical[4])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.3',name.clinical[5])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.3',name.clinical[6])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.3',name.clinical[7])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.3',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$METABRIC.pred.4=round(c(chisq.test(table(mat.METABRIC.add[,c('pred.4',name.clinical[1])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.4',name.clinical[2])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.4',name.clinical[3])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.4',name.clinical[4])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.4',name.clinical[5])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.4',name.clinical[6])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.4',name.clinical[7])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.4',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$METABRIC.pred.5=round(c(chisq.test(table(mat.METABRIC.add[,c('pred.5',name.clinical[1])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.5',name.clinical[2])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.5',name.clinical[3])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.5',name.clinical[4])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.5',name.clinical[5])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.5',name.clinical[6])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.5',name.clinical[7])]))$p.value, 
                                               chisq.test(table(mat.METABRIC.add[,c('pred.5',name.clinical[8])]))$p.value),3)
#
mat.clinical.test.pred$TCGA.pred.1=round(c(chisq.test(table(mat.TCGA.add[,c('pred.1',name.clinical[1])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.1',name.clinical[2])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.1',name.clinical[3])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.1',name.clinical[4])]))$p.value, NA, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.1',name.clinical[6])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.1',name.clinical[7])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.1',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$TCGA.pred.2=round(c(chisq.test(table(mat.TCGA.add[,c('pred.2',name.clinical[1])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.2',name.clinical[2])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.2',name.clinical[3])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.2',name.clinical[4])]))$p.value, NA,  
                                           chisq.test(table(mat.TCGA.add[,c('pred.2',name.clinical[6])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.2',name.clinical[7])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.2',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$TCGA.pred.3=round(c(chisq.test(table(mat.TCGA.add[,c('pred.3',name.clinical[1])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.3',name.clinical[2])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.3',name.clinical[3])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.3',name.clinical[4])]))$p.value,NA,  
                                           chisq.test(table(mat.TCGA.add[,c('pred.3',name.clinical[6])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.3',name.clinical[7])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.3',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$TCGA.pred.4=round(c(chisq.test(table(mat.TCGA.add[,c('pred.4',name.clinical[1])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.4',name.clinical[2])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.4',name.clinical[3])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.4',name.clinical[4])]))$p.value, NA, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.4',name.clinical[6])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.4',name.clinical[7])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.4',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$TCGA.pred.5=round(c(chisq.test(table(mat.TCGA.add[,c('pred.5',name.clinical[1])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.5',name.clinical[2])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.5',name.clinical[3])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.5',name.clinical[4])]))$p.value, NA,  
                                           chisq.test(table(mat.TCGA.add[,c('pred.5',name.clinical[6])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.5',name.clinical[7])]))$p.value, 
                                           chisq.test(table(mat.TCGA.add[,c('pred.5',name.clinical[8])]))$p.value),3)
#
mat.clinical.test.pred$ICGC.pred.1=round(c(chisq.test(table(mat.ICGC.add[,c('pred.1',name.clinical[1])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.1',name.clinical[2])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.1',name.clinical[3])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.1',name.clinical[4])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.1',name.clinical[5])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.1',name.clinical[6])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.1',name.clinical[7])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.1',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$ICGC.pred.2=round(c(chisq.test(table(mat.ICGC.add[,c('pred.2',name.clinical[1])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.2',name.clinical[2])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.2',name.clinical[3])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.2',name.clinical[4])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.2',name.clinical[5])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.2',name.clinical[6])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.2',name.clinical[7])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.2',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$ICGC.pred.3=round(c(chisq.test(table(mat.ICGC.add[,c('pred.3',name.clinical[1])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.3',name.clinical[2])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.3',name.clinical[3])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.3',name.clinical[4])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.3',name.clinical[5])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.3',name.clinical[6])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.3',name.clinical[7])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.3',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$ICGC.pred.4=round(c(chisq.test(table(mat.ICGC.add[,c('pred.4',name.clinical[1])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.4',name.clinical[2])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.4',name.clinical[3])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.4',name.clinical[4])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.4',name.clinical[5])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.4',name.clinical[6])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.4',name.clinical[7])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.4',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$ICGC.pred.5=round(c(chisq.test(table(mat.ICGC.add[,c('pred.5',name.clinical[1])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.5',name.clinical[2])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.5',name.clinical[3])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.5',name.clinical[4])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.5',name.clinical[5])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.5',name.clinical[6])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.5',name.clinical[7])]))$p.value, 
                                           chisq.test(table(mat.ICGC.add[,c('pred.5',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$TransNEO.pred.1=round(c(chisq.test(table(mat.TransNEO.add[,c('pred.1',name.clinical[1])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.1',name.clinical[2])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.1',name.clinical[3])]))$p.value,NA,
                                               chisq.test(table(mat.TransNEO.add[,c('pred.1',name.clinical[5])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.1',name.clinical[6])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.1',name.clinical[7])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.1',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$TransNEO.pred.2=round(c(chisq.test(table(mat.TransNEO.add[,c('pred.2',name.clinical[1])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.2',name.clinical[2])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.2',name.clinical[3])]))$p.value,NA, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.2',name.clinical[5])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.2',name.clinical[6])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.2',name.clinical[7])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.2',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$TransNEO.pred.3=round(c(chisq.test(table(mat.TransNEO.add[,c('pred.3',name.clinical[1])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.3',name.clinical[2])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.3',name.clinical[3])]))$p.value,NA,
                                               chisq.test(table(mat.TransNEO.add[,c('pred.3',name.clinical[5])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.3',name.clinical[6])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.3',name.clinical[7])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.3',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$TransNEO.pred.4=round(c(chisq.test(table(mat.TransNEO.add[,c('pred.4',name.clinical[1])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.4',name.clinical[2])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.4',name.clinical[3])]))$p.value,NA, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.4',name.clinical[5])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.4',name.clinical[6])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.4',name.clinical[7])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.4',name.clinical[8])]))$p.value),3)
mat.clinical.test.pred$TransNEO.pred.5=round(c(chisq.test(table(mat.TransNEO.add[,c('pred.5',name.clinical[1])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.5',name.clinical[2])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.5',name.clinical[3])]))$p.value,NA, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.5',name.clinical[5])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.5',name.clinical[6])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.5',name.clinical[7])]))$p.value, 
                                               chisq.test(table(mat.TransNEO.add[,c('pred.5',name.clinical[8])]))$p.value),3)
#
mat.clinical.test.pred$NEWTON.pred.1=round(c(chisq.test(table(mat.NEWTON.add[,c('pred.1',name.clinical[1])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.1',name.clinical[2])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.1',name.clinical[3])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.1',name.clinical[4])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.1',name.clinical[5])]))$p.value,NA,NA,NA),3)
mat.clinical.test.pred$NEWTON.pred.3=round(c(chisq.test(table(mat.NEWTON.add[,c('pred.3',name.clinical[1])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.3',name.clinical[2])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.3',name.clinical[3])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.3',name.clinical[4])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.3',name.clinical[5])]))$p.value,NA,NA,NA),3)
mat.clinical.test.pred$NEWTON.pred.4=round(c(chisq.test(table(mat.NEWTON.add[,c('pred.4',name.clinical[1])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.4',name.clinical[2])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.4',name.clinical[3])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.4',name.clinical[4])]))$p.value, 
                                             chisq.test(table(mat.NEWTON.add[,c('pred.4',name.clinical[5])]))$p.value,NA,NA,NA),3)


## mat.type.clinical.pred ####
name.clinical.n=c(2,2,2,2,3,2,5,10)
name.clinical.label=c('<50','>=50',
                      'ER+','ER-',
                      'Her2+','Her2-',
                      'PR+','PR-',
                      "I","II","III",
                      'LN+','LN-',
                      'Ba','H','A','B','N',
                      1:10)
mat.type.clinical.pred=data.frame(type=c(rep(name.clinical,name.clinical.n),
                                         rep(name.clinical,name.clinical.n)),
                                  label=c(name.clinical.label,name.clinical.label),
                                  HRD=c(rep('HRD+',length(name.clinical.label)),
                                        rep('HRD-',length(name.clinical.label))),
                                  color=rep(c("#9ac9db",'#ff8884',
                                              "#66C2A5",'#f8ac8c',
                                              "#82B0D2",'#ff8884',
                                              "#66C2A5",'#f8ac8c',
                                              '#CFEAF1','#B883D4','#FA7F6F',
                                              "#66C2A5",'#f8ac8c', 
                                              colipam,
                                              coliClust),2)) 
mat.type.clinical.pred$HRD=factor(mat.type.clinical.pred$HRD,levels =c('HRD+','HRD-')) 

#add prob;
mat.type.clinical.pred$METABRIC.pred.1=c(prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==1),name.clinical[1]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==1),name.clinical[2]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==1),name.clinical[3]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==1),name.clinical[4]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==1),name.clinical[5]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==1),name.clinical[6]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==1),name.clinical[7]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==1),name.clinical[8]])), #HRD+
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==0),name.clinical[1]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==0),name.clinical[2]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==0),name.clinical[3]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==0),name.clinical[4]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==0),name.clinical[5]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==0),name.clinical[6]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==0),name.clinical[7]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.1==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$METABRIC.pred.2=c(prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==1),name.clinical[1]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==1),name.clinical[2]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==1),name.clinical[3]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==1),name.clinical[4]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==1),name.clinical[5]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==1),name.clinical[6]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==1),name.clinical[7]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==1),name.clinical[8]])), #HRD+
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==0),name.clinical[1]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==0),name.clinical[2]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==0),name.clinical[3]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==0),name.clinical[4]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==0),name.clinical[5]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==0),name.clinical[6]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==0),name.clinical[7]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.2==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$METABRIC.pred.3=c(prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==1),name.clinical[1]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==1),name.clinical[2]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==1),name.clinical[3]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==1),name.clinical[4]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==1),name.clinical[5]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==1),name.clinical[6]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==1),name.clinical[7]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==1),name.clinical[8]])), #HRD+
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==0),name.clinical[1]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==0),name.clinical[2]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==0),name.clinical[3]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==0),name.clinical[4]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==0),name.clinical[5]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==0),name.clinical[6]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==0),name.clinical[7]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.3==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$METABRIC.pred.4=c(prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==1),name.clinical[1]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==1),name.clinical[2]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==1),name.clinical[3]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==1),name.clinical[4]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==1),name.clinical[5]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==1),name.clinical[6]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==1),name.clinical[7]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==1),name.clinical[8]])), #HRD+
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==0),name.clinical[1]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==0),name.clinical[2]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==0),name.clinical[3]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==0),name.clinical[4]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==0),name.clinical[5]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==0),name.clinical[6]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==0),name.clinical[7]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.4==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$METABRIC.pred.5=c(prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==1),name.clinical[1]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==1),name.clinical[2]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==1),name.clinical[3]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==1),name.clinical[4]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==1),name.clinical[5]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==1),name.clinical[6]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==1),name.clinical[7]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==1),name.clinical[8]])), #HRD+
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==0),name.clinical[1]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==0),name.clinical[2]])),
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==0),name.clinical[3]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==0),name.clinical[4]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==0),name.clinical[5]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==0),name.clinical[6]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==0),name.clinical[7]])), 
                                         prop.table(table(mat.METABRIC.add[which(mat.METABRIC.add$pred.5==0),name.clinical[8]])))#HRD-
#
mat.type.clinical.pred$TCGA.pred.1=c(prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==1),name.clinical[1]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==1),name.clinical[2]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==1),name.clinical[3]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==1),name.clinical[4]])), 
                                     NA,NA,NA,  
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==1),name.clinical[6]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==1),name.clinical[7]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==1),name.clinical[8]])), #HRD+
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==0),name.clinical[1]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==0),name.clinical[2]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==0),name.clinical[3]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==0),name.clinical[4]])),
                                     NA,NA,NA, 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==0),name.clinical[6]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==0),name.clinical[7]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.1==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$TCGA.pred.2=c(prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==1),name.clinical[1]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==1),name.clinical[2]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==1),name.clinical[3]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==1),name.clinical[4]])),
                                     NA,NA,NA, 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==1),name.clinical[6]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==1),name.clinical[7]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==1),name.clinical[8]])), #HRD+
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==0),name.clinical[1]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==0),name.clinical[2]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==0),name.clinical[3]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==0),name.clinical[4]])),
                                     NA,NA,NA, 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==0),name.clinical[6]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==0),name.clinical[7]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.2==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$TCGA.pred.3=c(prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==1),name.clinical[1]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==1),name.clinical[2]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==1),name.clinical[3]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==1),name.clinical[4]])),
                                     NA,NA,NA, 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==1),name.clinical[6]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==1),name.clinical[7]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==1),name.clinical[8]])), #HRD+
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==0),name.clinical[1]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==0),name.clinical[2]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==0),name.clinical[3]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==0),name.clinical[4]])),
                                     NA,NA,NA, 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==0),name.clinical[6]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==0),name.clinical[7]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.3==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$TCGA.pred.4=c(prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==1),name.clinical[1]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==1),name.clinical[2]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==1),name.clinical[3]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==1),name.clinical[4]])),
                                     NA,NA,NA, 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==1),name.clinical[6]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==1),name.clinical[7]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==1),name.clinical[8]])), #HRD+
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==0),name.clinical[1]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==0),name.clinical[2]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==0),name.clinical[3]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==0),name.clinical[4]])), 
                                     NA,NA,NA, 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==0),name.clinical[6]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==0),name.clinical[7]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.4==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$TCGA.pred.5=c(prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==1),name.clinical[1]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==1),name.clinical[2]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==1),name.clinical[3]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==1),name.clinical[4]])),
                                     NA,NA,NA, 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==1),name.clinical[6]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==1),name.clinical[7]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==1),name.clinical[8]])), #HRD+
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==0),name.clinical[1]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==0),name.clinical[2]])),
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==0),name.clinical[3]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==0),name.clinical[4]])), 
                                     NA,NA,NA, 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==0),name.clinical[6]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==0),name.clinical[7]])), 
                                     prop.table(table(mat.TCGA.add[which(mat.TCGA.add$pred.5==0),name.clinical[8]])))#HRD-
#
mat.type.clinical.pred$ICGC.pred.1=c(prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==1),name.clinical[1]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==1),name.clinical[2]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==1),name.clinical[3]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==1),name.clinical[4]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==1),name.clinical[5]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==1),name.clinical[6]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==1),name.clinical[7]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==1),name.clinical[8]])), #HRD+
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==0),name.clinical[1]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==0),name.clinical[2]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==0),name.clinical[3]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==0),name.clinical[4]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==0),name.clinical[5]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==0),name.clinical[6]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==0),name.clinical[7]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.1==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$ICGC.pred.2=c(prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==1),name.clinical[1]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==1),name.clinical[2]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==1),name.clinical[3]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==1),name.clinical[4]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==1),name.clinical[5]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==1),name.clinical[6]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==1),name.clinical[7]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==1),name.clinical[8]])), #HRD+
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==0),name.clinical[1]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==0),name.clinical[2]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==0),name.clinical[3]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==0),name.clinical[4]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==0),name.clinical[5]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==0),name.clinical[6]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==0),name.clinical[7]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.2==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$ICGC.pred.3=c(prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==1),name.clinical[1]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==1),name.clinical[2]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==1),name.clinical[3]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==1),name.clinical[4]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==1),name.clinical[5]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==1),name.clinical[6]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==1),name.clinical[7]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==1),name.clinical[8]])), #HRD+
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==0),name.clinical[1]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==0),name.clinical[2]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==0),name.clinical[3]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==0),name.clinical[4]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==0),name.clinical[5]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==0),name.clinical[6]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==0),name.clinical[7]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.3==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$ICGC.pred.4=c(prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==1),name.clinical[1]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==1),name.clinical[2]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==1),name.clinical[3]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==1),name.clinical[4]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==1),name.clinical[5]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==1),name.clinical[6]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==1),name.clinical[7]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==1),name.clinical[8]])), #HRD+
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==0),name.clinical[1]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==0),name.clinical[2]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==0),name.clinical[3]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==0),name.clinical[4]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==0),name.clinical[5]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==0),name.clinical[6]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==0),name.clinical[7]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.4==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$ICGC.pred.5=c(prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==1),name.clinical[1]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==1),name.clinical[2]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==1),name.clinical[3]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==1),name.clinical[4]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==1),name.clinical[5]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==1),name.clinical[6]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==1),name.clinical[7]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==1),name.clinical[8]])), #HRD+
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==0),name.clinical[1]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==0),name.clinical[2]])),
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==0),name.clinical[3]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==0),name.clinical[4]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==0),name.clinical[5]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==0),name.clinical[6]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==0),name.clinical[7]])), 
                                     prop.table(table(mat.ICGC.add[which(mat.ICGC.add$pred.5==0),name.clinical[8]])))#HRD-
#
mat.type.clinical.pred$TransNEO.pred.1=c(prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==1),name.clinical[1]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==1),name.clinical[2]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==1),name.clinical[3]])), 
                                         NA,NA,
                                         NA,prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==1),name.clinical[5]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==1),name.clinical[6]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==1),name.clinical[7]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==1),name.clinical[8]])), #HRD+
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==0),name.clinical[1]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==0),name.clinical[2]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==0),name.clinical[3]])), 
                                         NA,NA,
                                         NA,prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==0),name.clinical[5]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==0),name.clinical[6]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==0),name.clinical[7]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.1==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$TransNEO.pred.2=c(prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==1),name.clinical[1]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==1),name.clinical[2]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==1),name.clinical[3]])),
                                         NA,NA,
                                         NA,prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==1),name.clinical[5]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==1),name.clinical[6]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==1),name.clinical[7]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==1),name.clinical[8]])), #HRD+
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==0),name.clinical[1]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==0),name.clinical[2]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==0),name.clinical[3]])), 
                                         NA,NA,
                                         NA,prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==0),name.clinical[5]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==0),name.clinical[6]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==0),name.clinical[7]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.2==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$TransNEO.pred.3=c(prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==1),name.clinical[1]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==1),name.clinical[2]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==1),name.clinical[3]])), 
                                         NA,NA,
                                         NA,prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==1),name.clinical[5]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==1),name.clinical[6]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==1),name.clinical[7]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==1),name.clinical[8]])), #HRD+
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==0),name.clinical[1]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==0),name.clinical[2]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==0),name.clinical[3]])), 
                                         NA,NA,
                                         NA,prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==0),name.clinical[5]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==0),name.clinical[6]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==0),name.clinical[7]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.3==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$TransNEO.pred.4=c(prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==1),name.clinical[1]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==1),name.clinical[2]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==1),name.clinical[3]])), 
                                         NA,NA,
                                         NA,prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==1),name.clinical[5]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==1),name.clinical[6]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==1),name.clinical[7]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==1),name.clinical[8]])), #HRD+
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==0),name.clinical[1]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==0),name.clinical[2]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==0),name.clinical[3]])), 
                                         NA,NA,
                                         NA,prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==0),name.clinical[5]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==0),name.clinical[6]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==0),name.clinical[7]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.4==0),name.clinical[8]])))#HRD-
mat.type.clinical.pred$TransNEO.pred.5=c(prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==1),name.clinical[1]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==1),name.clinical[2]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==1),name.clinical[3]])), 
                                         NA,NA,
                                         NA,prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==1),name.clinical[5]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==1),name.clinical[6]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==1),name.clinical[7]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==1),name.clinical[8]])), #HRD+
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==0),name.clinical[1]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==0),name.clinical[2]])),
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==0),name.clinical[3]])), 
                                         NA,NA,
                                         NA,prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==0),name.clinical[5]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==0),name.clinical[6]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==0),name.clinical[7]])), 
                                         prop.table(table(mat.TransNEO.add[which(mat.TransNEO.add$pred.5==0),name.clinical[8]])))#HRD-
#
mat.type.clinical.pred$NEWTON.pred.1=c(prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.1==1),name.clinical[1]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.1==1),name.clinical[2]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.1==1),name.clinical[3]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.1==1),name.clinical[4]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.1==1),name.clinical[5]])), 
                                       rep(NA,17), #HRD+
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.1==0),name.clinical[1]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.1==0),name.clinical[2]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.1==0),name.clinical[3]])),  
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.1==0),name.clinical[4]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.1==0),name.clinical[5]])), 
                                       rep(NA,17))#HRD-
mat.type.clinical.pred$NEWTON.pred.3=c(prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.3==1),name.clinical[1]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.3==1),name.clinical[2]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.3==1),name.clinical[3]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.3==1),name.clinical[4]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.3==1),name.clinical[5]])),
                                       rep(NA,17),#HRD+
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.3==0),name.clinical[1]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.3==0),name.clinical[2]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.3==0),name.clinical[3]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.3==0),name.clinical[4]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.3==0),name.clinical[5]])),
                                       rep(NA,17))#HRD-
mat.type.clinical.pred$NEWTON.pred.4=c(prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.4==1),name.clinical[1]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.4==1),name.clinical[2]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.4==1),name.clinical[3]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.4==1),name.clinical[4]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.4==1),name.clinical[5]])),
                                       rep(NA,17),#HRD+
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.4==0),name.clinical[1]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.4==0),name.clinical[2]])),
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.4==0),name.clinical[3]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.4==0),name.clinical[4]])), 
                                       prop.table(table(mat.NEWTON.add[which(mat.NEWTON.add$pred.4==0),name.clinical[5]])),
                                       rep(NA,17))#HRD-



### barplot ####
baseSize=13

#--METABRIC  
#--pred.1
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$METABRIC.pred.1),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.METABRIC.add$pred.1,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
#
gg.bar.clinical.METABRIC.pred.1=ggplot(tmp,
                                       aes(y=METABRIC.pred.1,x=type))+
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+ 
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) +  
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

#--pred.2
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.METABRIC.add$pred.2,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$METABRIC.pred.2),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.METABRIC.pred.2=ggplot(tmp,
                                       aes(y=METABRIC.pred.2,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) +  
  scale_fill_manual(values=tmp$color)+ 
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) +  
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

#--pred.3
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.METABRIC.add$pred.3,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$METABRIC.pred.3),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.METABRIC.pred.3=ggplot(tmp,
                                       aes(y=METABRIC.pred.3,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.4
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.METABRIC.add$pred.4,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$METABRIC.pred.4),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.METABRIC.pred.4=ggplot(tmp,aes(y=METABRIC.pred.4,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.5
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.METABRIC.add$pred.5,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$METABRIC.pred.5),]
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.METABRIC.pred.5=ggplot(tmp,aes(y=METABRIC.pred.5,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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


#--TCGA 
#--pred.1
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$TCGA.pred.1),]
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
#
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.TCGA.add$pred.1,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
#
gg.bar.clinical.TCGA.pred.1=ggplot(tmp,
                                   aes(y=TCGA.pred.1,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.2
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.TCGA.add$pred.2,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$TCGA.pred.2),]
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.TCGA.pred.2=ggplot(tmp,
                                   aes(y=TCGA.pred.2,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.3
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.TCGA.add$pred.3,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$TCGA.pred.3),]
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.TCGA.pred.3=ggplot(tmp,
                                   aes(y=TCGA.pred.3,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.4
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.TCGA.add$pred.4,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$TCGA.pred.4),]
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.TCGA.pred.4=ggplot(tmp,aes(y=TCGA.pred.4,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.5
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.TCGA.add$pred.5,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$TCGA.pred.5),]
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.TCGA.pred.5=ggplot(tmp,aes(y=TCGA.pred.5,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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
#--pred.1
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$ICGC.pred.1),]
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
#
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.ICGC.add$pred.1,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
#
gg.bar.clinical.ICGC.pred.1=ggplot(tmp,
                                   aes(y=ICGC.pred.1,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.2
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.ICGC.add$pred.2,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$ICGC.pred.2),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.ICGC.pred.2=ggplot(tmp,
                                   aes(y=ICGC.pred.2,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.3
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.ICGC.add$pred.3,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$ICGC.pred.3),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.ICGC.pred.3=ggplot(tmp,
                                   aes(y=ICGC.pred.3,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.4
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.ICGC.add$pred.4,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$ICGC.pred.4),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.ICGC.pred.4=ggplot(tmp,aes(y=ICGC.pred.4,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.5
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.ICGC.add$pred.5,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$ICGC.pred.5),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
gg.bar.clinical.ICGC.pred.5=ggplot(tmp,aes(y=ICGC.pred.5,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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


#--TransNEO 
#--pred.1
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$TransNEO.pred.1),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.TransNEO.add$pred.1,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
gg.bar.clinical.TransNEO.pred.1=ggplot(tmp,
                                       aes(y=TransNEO.pred.1,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.2
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.TransNEO.add$pred.2,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
gg.bar.clinical.TransNEO.pred.2=ggplot(tmp,
                                       aes(y=TransNEO.pred.2,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.3
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.TransNEO.add$pred.3,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
gg.bar.clinical.TransNEO.pred.3=ggplot(tmp,
                                       aes(y=TransNEO.pred.3,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.4
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.TransNEO.add$pred.4,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
gg.bar.clinical.TransNEO.pred.4=ggplot(tmp,
                                       aes(y=TransNEO.pred.4,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.5
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.TransNEO.add$pred.5,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
gg.bar.clinical.TransNEO.pred.5=ggplot(tmp,
                                       aes(y=TransNEO.pred.5,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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


#--NEWTON 
#--pred.1
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$NEWTON.pred.1),] 
tmp$type=factor(tmp$type,levels = name.clinical) 
tmp$label=factor(tmp$label,levels = tmp$label)
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.NEWTON.add$pred.1,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
gg.bar.clinical.NEWTON.pred.1=ggplot(tmp,
                                     aes(y=NEWTON.pred.1,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.3
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$NEWTON.pred.3),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
#
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.NEWTON.add$pred.3,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
gg.bar.clinical.NEWTON.pred.3=ggplot(tmp,
                                     aes(y=NEWTON.pred.3,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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

#--pred.4
tmp=mat.type.clinical.pred[!is.na(mat.type.clinical.pred$NEWTON.pred.4),] 
tmp$type=factor(tmp$type,levels = name.clinical)
tmp$label=factor(tmp$label,levels = tmp$label[1:(dim(tmp)[1]/2)])
#
pred.HRD.labs=paste0(c('HRD+ (N=','HRD- (N='), 
                     as.numeric(table(factor(mat.NEWTON.add$pred.4,levels=c(1,0)))),')')
names(pred.HRD.labs)=c('HRD+','HRD-')
gg.bar.clinical.NEWTON.pred.4=ggplot(tmp,
                                     aes(y=NEWTON.pred.4,x=type))+
  geom_bar(stat="identity",  
           position = position_stack(reverse = TRUE),
           aes(fill=label),color="gray30",size=0.2) + 
  scale_fill_manual(values=tmp$color)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust =0.5,reverse =T)
  ) +
  facet_grid(~HRD,labeller = labeller(HRD=pred.HRD.labs)) + 
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



### combine #### 
pdf(file = "fig_clinical_predHRD.pdf", width =18, height =20)
ggarrange(NULL,ggarrange(gg.bar.clinical.METABRIC.pred.1,gg.bar.clinical.METABRIC.pred.2,gg.bar.clinical.METABRIC.pred.3, 
                         gg.bar.clinical.METABRIC.pred.4,gg.bar.clinical.METABRIC.pred.5,
                         labels=name.model,ncol=n.model,nrow=1,label.y = 1.05),
          NULL,ggarrange(gg.bar.clinical.TCGA.pred.1,gg.bar.clinical.TCGA.pred.2,gg.bar.clinical.TCGA.pred.3, 
                         gg.bar.clinical.TCGA.pred.4,gg.bar.clinical.TCGA.pred.5,
                         labels=name.model,ncol=n.model,nrow=1,label.y = 1.05),
          NULL,ggarrange(gg.bar.clinical.ICGC.pred.1,gg.bar.clinical.ICGC.pred.2,gg.bar.clinical.ICGC.pred.3, 
                         gg.bar.clinical.ICGC.pred.4,gg.bar.clinical.ICGC.pred.5,
                         labels=name.model,ncol=n.model,nrow=1,label.y = 1.05), 
          NULL,
          ggarrange(gg.bar.clinical.TransNEO.pred.1,gg.bar.clinical.TransNEO.pred.2,gg.bar.clinical.TransNEO.pred.3, 
                    gg.bar.clinical.TransNEO.pred.4,gg.bar.clinical.TransNEO.pred.5,
                    labels=name.model,ncol=n.model,nrow=1,label.y = 1.05),
          NULL,
          ggarrange(gg.bar.clinical.NEWTON.pred.1,NULL,gg.bar.clinical.NEWTON.pred.3, 
                    gg.bar.clinical.NEWTON.pred.4,NULL,
                    labels=c("CNA",'',"SNV","SNV+CNA",''),ncol=n.model,nrow=1,label.y = 1.05), 
          labels=c("METABRIC",NA, "  TCGA",NA, "   ICGC",NA,"TransNEO",NA," NEWTON"),
          ncol=1,nrow=10,
          heights =rep(c(0.11,1),5)) #name.cohort[1:3],c(0.08,1,0.08,1,0.08,1)
dev.off()




# save #### 
mat.SCANB.add$cohort='SCANB'
mat.GEL.add$cohort='GEL'
mat.GEL.add$id=paste0('id.',1:dim(mat.GEL.add)[1])
#--all predicted prob, CCLE only brca; mat.CCLE_brca.add
mat.prob.cohort=rbind(mat.TCGA.add[,c('cohort','id','is.hrd','prob.1','prob.2','prob.3','prob.4','prob.5')],
                      mat.ICGC.add[,c('cohort','id','is.hrd','prob.1','prob.2','prob.3','prob.4','prob.5')],
                      mat.METABRIC.add[,c('cohort','id','is.hrd','prob.1','prob.2','prob.3','prob.4','prob.5')],
                      mat.GEL.add[,c('cohort','id','is.hrd','prob.1','prob.2','prob.3','prob.4','prob.5')],
                      mat.SCANB.add[,c('cohort','id','is.hrd','prob.1','prob.2','prob.3','prob.4','prob.5')],
                      mat.TransNEO.add[,c('cohort','id','is.hrd','prob.1','prob.2','prob.3','prob.4','prob.5')],
                      mat.CCLE.add[which(mat.CCLE.add$CANCER_TYPE=='Breast Cancer'),
                                   c('cohort','id','is.hrd','prob.1','prob.2','prob.3','prob.4','prob.5')])
#
tmp=mat.NEWTON.add[,c('cohort','id','prob.1','prob.3','prob.4')]
tmp$is.hrd=NA
tmp$prob.2=NA
tmp$prob.5=NA
tmp=tmp[,c(1,2,6,3,7,4,5,8)]
colnames(tmp)==colnames(mat.prob.cohort)
mat.prob.cohort=rbind(mat.prob.cohort,tmp)  
colnames(mat.prob.cohort)[3:8]=c('trueHRD',name.model)  
write.csv(mat.prob.cohort, file = "output/mat.prob.cohort_rf.csv", row.names = FALSE)    



#--
save.image('output/validation.RData') # 
