

# parameter ####
#--5 models
n.model=5
name.model=c('CNA','ASCN','SNV','SNV+CNA','SNV+ASCN')


#--index of performance
name.perf.v1=c('AUC','AUCPR',
               paste0(c('Sensitivity','Specificity','Precision','MCC','threshold'),'.f1'))


#--established colours for the Integrative clusters: 
coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', 
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')

#--established colours for Pam50:
colipam <- c("#E41A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")


#--clinical
name.clinical=c('Age','ER','Her2','PR','Grade','LN','Pam50','IntClust') #8


#--optimal threshold obtained from LOOCV and F1 score
mat.cut.f1=data.frame(variable=name.model,cut=c(0.382,0.450,0.344, 0.296, 0.374))
mat.cut.f1$variable=factor(mat.cut.f1$variable,levels=name.model)


#--bar plot for BRCA1/2 status
name.col=data.frame(color=c('No methylation&At least one BRCA1/2 with loh(=1)&No BRCA1/2 mutation detected', #hrd NA
                            'BRCA1/2 mutation without loh', #hrd NA
                            'No BRCA1/2 mutation detected or BRCA1/2 without loh',#hrd-
                            'Germline BRCA1 mutation with loh', #hrd+
                            'Germline BRCA2 mutation with loh',
                            'Somatic BRCA1 mutation with loh',
                            'Somatic BRCA2 mutation with loh',
                            'BRCA1 promoter methylation with loh'),
                    col=c('grey3','black',"grey",'green4','blue4','green','royalblue','purple'))
name.col$color=factor(name.col$color) 


# function ####
#--matrix of the p-value of the correlation
# mat: is a matrix of data 
cor.mtest <- function(mat) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j]) 
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


#--performance on threshold from f1 score
#prob: is a data.frame
#truelabel: factor with levels=c(1,0) 
perf.prob.f1=function(prob,truelabel){  
  predic=ROCR::prediction(prob,truelabel) 
  temp=list()
  temp$predic=predic 
  temp$prob=prob 
  temp=list(roc  =ROCR::performance(predic,"tpr","fpr"), #y,x
            pr     =ROCR::performance(predic, "prec", "rec"), # precision/recall curve (x-axis: recall, y-axis: precision)
            auc    =ROCR::performance(predic,"auc")@y.values[[1]], 
            aucpr  =ROCR::performance(predic,"aucpr")@y.values[[1]]
  )  
  temp$score.f1     =2/(1/temp$pr@x.values[[1]]+1/temp$pr@y.values[[1]]) #f1 score 
  threshold.f1=temp$pr@alpha.values[[1]][which.max(temp$score.f1)]
  temp$pred.f1=factor(ifelse(prob>=threshold.f1,1,0),levels=c(1,0))  
  temp$index.f1=c(as.numeric(confusionMatrix(temp$pred.f1,truelabel)$byClass[c(1,2,5)]),
                  mltools::mcc(preds =temp$pred.f1, actuals=truelabel),
                  threshold.f1) 
  temp$perf=c(temp$auc,temp$aucpr,temp$index.f1) 
  names(temp$perf)=name.perf.v1   
  return(temp)
}


#--roc plot for 5 models
#x: list of 5 from roc.curve()
#pos: position, e.g. "bottomright"
#cohort: cohort name, e.g. 'TCGA'
plot.roc.5=function(x,pos,cohort){ 
  plot(x[[1]],col=3,rand.plot=T,auc.main = FALSE,legend = FALSE,main=cohort)
  plot(x[[2]],col=4,add=T)
  plot(x[[3]],col=7,add=T)
  plot(x[[4]],col=2,add=T)
  plot(x[[5]],col=1,add=T) 
  legend(pos,col=c(3,4,7,2,1),lty=1, lwd = 2, 
         legend=c(paste0(name.model[1],'(AUC=',round(x[[1]]$auc,3),')'),
                  paste0(name.model[2],'(AUC=',round(x[[2]]$auc,3),')'),
                  paste0(name.model[3],'(AUC=',round(x[[3]]$auc,3),')'),
                  paste0(name.model[4],'(AUC=',round(x[[4]]$auc,3),')'),
                  paste0(name.model[5],'(AUC=',round(x[[5]]$auc,3),')') 
         ))
}



#--pr plot for 5 models
#x: list of 5 from pr.curve()
#pos: position, e.g. "bottomright"
#cohort: cohort name, e.g. 'TCGA'
plot.pr.5=function(x,pos,cohort){ 
  plot(x[[1]],col=3,rand.plot=T,auc.main = FALSE,legend = FALSE,
       main=paste0(cohort,' (baseline=',round(x[[1]]$rand$auc.integral,3),')'))
  plot(x[[2]],col=4,add=T)
  plot(x[[3]],col=7,add=T)
  plot(x[[4]],col=2,add=T)
  plot(x[[5]],col=1,add=T) 
  legend(pos,col=c(3,4,7,2,1),lty=1,lwd= 2,  
         legend=c(paste0(name.model[1],'(AUCPR=',round(x[[1]]$auc.integral,3),')'),
                  paste0(name.model[2],'(AUCPR=',round(x[[2]]$auc.integral,3),')'),
                  paste0(name.model[3],'(AUCPR=',round(x[[3]]$auc.integral,3),')'),
                  paste0(name.model[4],'(AUCPR=',round(x[[4]]$auc.integral,3),')'),
                  paste0(name.model[5],'(AUCPR=',round(x[[5]]$auc.integral,3),')') 
         )) 
}


#--5 bar plots of predicted probabilities with given threshold
#prob: list of predicted probabilities
#mat: data frame that contains column ('id','color1'), where 'color1' shows the BRCA1/2 status 
#cut: optimal threshold in mat.cut.f1$cut
plot.bar.5.given=function(prob,mat,n.model=5,name.col=name.col,cut){ 
  gg.prob=list()
  for(m in 1:n.model){
    print(m)
    #--train
    tmp=mat[,c('id','color1')]  
    plotdf=data.frame(prob=prob[[m]], 
                      type=tmp$color1) 
    plotdf=plotdf[order(plotdf$prob),]  
    plotdf$loc=1:dim(plotdf)[1]
    #
    gg.prob[[m]]=ggplot(plotdf, aes(x=loc, y=prob, fill=type)) +
      geom_bar(stat="identity")+#theme_minimal()+
      geom_hline(yintercept=cut[m], linetype="dashed", color = "grey")+
      ylab(name.model[m])+xlab('')+
      theme(axis.title.x = element_blank(), 
            legend.position="bottom",
            legend.key.size = unit(0.5, 'cm'), #change legend key size
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.title = element_text(size=5), #change legend title font size
            legend.text = element_text(size=5))  + #change legend text font size
      scale_fill_manual(values =name.col$col[match(levels(factor(mat$color1)),name.col$color)]) 
  }
  return(gg.prob) 
}


#--hazard ratio results from coxph()
#x: model from coxph()
hazard_table=function(x){ 
  HR <- round(exp(coef(x)), 3)
  CI <- round(exp(confint(x)), 3)
  P <- round(coef(summary(x))[,5], 3)
  colnames(CI) <- c("Lower", "Higher") 
  formulas=rep(paste(as.character(x$formula)[2],as.character(x$formula)[1],as.character(x$formula)[3]),
               length(HR))
  n=x$n
  nevent=x$nevent
  hazard_table <- as.data.frame(cbind(HR, CI, P,formulas,n,nevent))
  hazard_table<-cbind(feature=row.names(hazard_table),hazard_table) 
  hazard_table$p.sig='ns'
  hazard_table$p.sig[hazard_table$HR> 1 & hazard_table$P < 0.05]="poor"
  hazard_table$p.sig[hazard_table$HR< 1 & hazard_table$P < 0.05]="good"  
  
  return(hazard_table)
}


#--odds ratio results from glm() logistic regression
#x: model from glm()
odds.ratio_table=function(x){ 
  lrm        <- summary(x)
  pval       <- lrm$coefficients[c(2:nrow(lrm$coefficients)),"Pr(>|z|)"]  
  odds.ratio <- exp(lrm$coefficients[c(2:nrow(lrm$coefficients)),"Estimate"])
  conf       <- suppressMessages(confint(x))
  ci.low     <- exp(conf[c(2:nrow(conf)),1])
  ci.hi      <- exp(conf[c(2:nrow(conf)),2])
  OR_table <- data.frame(feature=rownames(lrm$coefficients)[2:nrow(lrm$coefficients)],pval,odds.ratio,ci.low,ci.hi)
  OR_table$p.sig='ns'
  OR_table$p.sig[OR_table$odds.ratio> 1 & OR_table$pval < 0.05]="up"
  OR_table$p.sig[OR_table$odds.ratio< 1 & OR_table$pval < 0.05]="down"  
  return(OR_table)
}


