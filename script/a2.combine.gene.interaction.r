print("Combine.gene.interaction.R pathway.tsv interaction.tsv output")

rm(list=ls())
#define response variable
args = commandArgs(trailingOnly=TRUE)

library(stringr)
#read bacteria abundance file
#pathway=read.table('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/CAZymes_metagenome_with_description.tsv',row.names = 1,sep="\t",header = T)
#interaction=read.table('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/high_nutrition.interaction.tsv',sep="\t",header = T)

#pathway=read.table('D:/a-document/github/BCP/data/annotation/CAZymes_metagenome.tsv',row.names = 1,sep="\t",header = T)
#interaction=read.table('D:/a-document/github/BCP/data/interaction/b4.hn.interaction.csv',sep=",",header = F,skip=1)

pathway=read.table(args[1],row.names = 1,sep="\t",header = T)
interaction=read.table(args[2],sep=",",header = F,skip=1)

colnames(interaction)=c('species1','species2','interaction')
#
#print("a")
#head(interaction)
#tail(interaction)
#dim(interaction)
#head(pathway)
#print("aa")
pathway=pathway[rowSums(pathway)>0,]
#print("ab")
#dim(pathway)
#head(interaction)
interaction$species1=str_replace(interaction$species1,"Species","s")
interaction$species2=str_replace(interaction$species2,"Species","s")
#head(interaction)
#find optimal lambda value that minimizes test MSE
#pathway

#interaction[1] means the rows
#dim(pathway)[1]*2+1 means the feature of pathway of two species and add a label. 
totalnumber=dim(interaction)[1] * ((dim(pathway)[1]*2)+1)
#
#print(totalnumber)
#print("b")
xdata=matrix(rep(1,totalnumber),nrow=dim(interaction)[1])
#dim(xdata)
for (i in seq(dim(interaction)[1]) ){
  #print(i)
  rowvalue=pathway[,interaction[i,1]]
  rowvalue=append(rowvalue,pathway[,interaction[i,2]])
  rowvalue=append(rowvalue,interaction[i,3])
  
  xdata[i,]=rowvalue
  
}

#print("c")
#xdata
#xdata2=xdata
colnames(xdata)=c(paste('befor',rownames(pathway)),paste('after',rownames(pathway)),"interaction")
#colnames(xdata2)=paste('F',seq(dim(xdata)[2]),sep='')
ydata=interaction[,3]
#head(xdata)
#dim(xdata)

xdata=xdata[,colSums(xdata !=0) > 0]
#xdata2=xdata2[,colSums(xdata2 !=0) > 0]

write.csv(xdata,paste0(args[3],"ml.full.csv"),row.names = F)
#write.csv(xdata2,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/caz.train.csv',row.names = F)
write.csv(ydata,paste0(args[3],'ml.label.csv'))
##
#print("d")

totalnumber=dim(interaction)[1]*((dim(pathway)[1])+1)
#print(totalnumber)

xdata_sub=matrix(rep(1,totalnumber),nrow=dim(interaction)[1])
for (i in seq(dim(interaction)[1]) ){
  #print(i)
  rowvalue=pathway[,interaction[i,1]]-pathway[,interaction[i,2]]
  rowvalue=append(rowvalue,interaction[i,3])
  
  xdata_sub[i,]=rowvalue
  
}
#print("f")
colnames(xdata_sub)=c(rownames(pathway),"interaction")
xdata_sub=xdata_sub[,colSums(xdata_sub !=0) > 0]
#xdata_sub2=xdata_sub
#colnames(xdata_sub2)=append(paste('F',seq(dim(xdata_sub2)[2]-1),sep=''),'interaction')
write.csv(xdata_sub,paste0(args[3],".ml.sub.csv"),row.names = F)
#write.csv(xdata_sub2,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/pathannotation/caz.sub.train.csv',row.names = F)
##