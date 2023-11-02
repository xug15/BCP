rm(list=ls())
#define response variable


library(stringr)
#read bacteria abundance file
pathway=read.table('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/KEGG_pathways_MinPath_prunned.tsv',sep="\t",row.names = 1,header = T)
interaction=read.table('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/high_nutrition.interaction.tsv',sep="\t",header = T)
#
head(interaction)
tail(interaction)
dim(interaction)
head(pathway)
pathway=pathway[rowSums(pathway)>0,]
dim(pathway)
interaction$specise1=str_replace(interaction$specise1,"Species","s")
interaction$specise2=str_replace(interaction$specise2,"Species","s")
head(interaction)
#find optimal lambda value that minimizes test MSE
pathway


totalnumber=dim(interaction)[1]*((dim(pathway)[1]*2)+1)
print(totalnumber)

xdata=matrix(rep(1,totalnumber),nrow=dim(interaction)[1])
dim(xdata)
for (i in seq(dim(interaction)[1]) ){
  print(i)
  rowvalue=pathway[,interaction[i,1]]
  rowvalue=append(rowvalue,pathway[,interaction[i,2]])
  rowvalue=append(rowvalue,interaction[i,3])
  
  xdata[i,]=rowvalue
  
}


xdata
xdata2=xdata
colnames(xdata)=c(paste('befor',rownames(pathway)),paste('after',rownames(pathway)),"interaction")
colnames(xdata2)=paste('F',seq(dim(xdata)[2]),sep='')
ydata=interaction[,3]
head(xdata)
dim(xdata)

xdata=xdata[,colSums(xdata !=0) > 0]
xdata2=xdata2[,colSums(xdata2 !=0) > 0]

write.csv(xdata,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/kegg.train.fullname.csv',row.names = F)
write.csv(xdata2,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/kegg.train.csv',row.names = F)

write.csv(ydata,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/kegg.label.csv')

totalnumber=dim(interaction)[1]*((dim(pathway)[1])+1)
print(totalnumber)

xdata_sub=matrix(rep(1,totalnumber),nrow=dim(interaction)[1])
for (i in seq(dim(interaction)[1]) ){
  print(i)
  rowvalue=pathway[,interaction[i,1]]-pathway[,interaction[i,2]]
  rowvalue=append(rowvalue,interaction[i,3])
  
  xdata_sub[i,]=rowvalue
  
}
colnames(xdata_sub)=c(rownames(pathway),"interaction")
xdata_sub=xdata_sub[,colSums(xdata_sub !=0) > 0]
xdata_sub2=xdata_sub
colnames(xdata_sub2)=append(paste('F',seq(dim(xdata_sub2)[2]-1),sep=''),'interaction')
write.csv(xdata_sub,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/kegg.sub.train.fullname.csv',row.names = F)
write.csv(xdata_sub2,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/kegg.sub.train.csv',row.names = F)
##

##
#Run LASSO

#step2 fit the lasso regression model.
#install.packages('glmnet')
library(glmnet)
head(xdata)
xlast=dim(xdata)[2]
x=xdata[,-xlast]
head(x)
head(ydata)
y=ydata
#perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(x, y, alpha = 1)
#
best_lambda <- cv_model$lambda.min
best_lambda
#produce plot of test MSE by lambda value
plot(cv_model) 
#Step 3 Analyze Final Model
#find coefficients of best model
best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
coef(best_model)
#
#define new observation
new = matrix(c(24, 2.5, 3.5, 18.5), nrow=1, ncol=4) 
new = x[1,]
#use lasso regression model to predict response value
predict(best_model, s = best_lambda, newx = new)

y[1]

#use fitted best model to make predictions
y_predicted <- predict(best_model, s = best_lambda, newx = x)

#find SST and SSE
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq



