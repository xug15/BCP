rm(list=ls())
#define response variable


library(stringr)
#read bacteria abundance file
pathway=read.table('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/COG_metagenome_with_description.tsv',row.names = 1,sep="\t",header = T)
interaction=read.table('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/high_nutrition.interaction.tsv',sep="\t",header = T)
#
head(interaction)
tail(interaction)
dim(interaction)
head(pathway)
tail(pathway)
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

write.csv(xdata,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/cog.train.fullname.csv',row.names = F)
write.csv(xdata2,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/cog.train.csv',row.names = F)

write.csv(ydata,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/cog.label.csv')
##


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
write.csv(xdata_sub,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/cog.sub.train.fullname.csv',row.names = F)
write.csv(xdata_sub2,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/cog.sub.train.csv',row.names = F)
##
