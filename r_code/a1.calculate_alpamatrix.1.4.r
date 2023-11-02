#v1.3 to 1.4

rm(list=ls())

#read bacteria abundance file
abundance=read.csv('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/r_code/b1.b3.bacteria.abundance.csv',row.names = 1)
abundance=read.csv('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/r_code/b2.b456.bacteria.abundance.csv',row.names = 1)
abundance=read.csv('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/r_code/b2.b789.bacteria.abundance.csv',row.names = 1)
abundance=read.csv('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/r_code/b2.b101112.bacteria.abundance.csv',row.names = 1)
abundance=read.csv('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/r_code/b6.b131415.bacteria.abundance.csv',row.names = 1)
abundance=read.csv('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/r_code/b6.b161718.bacteria.abundance.csv',row.names = 1)
abundance=read.csv('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/r_code/b6.b192021.bacteria.abundance.csv',row.names = 1)
abundance=read.csv('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/r_code/b6.b222324.bacteria.abundance.csv.',row.names = 1)
abundance=read.csv('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/r_code/b6.b123456.bacteria.abundance.csv',row.names = 1)


head(abundance)
relative_abundance=abundance
#2. Calculate the percentage.
for (i in 1:dim(abundance)[2]){
  print(i)
  relative_abundance[,i]=abundance[,i]/colSums(abundance)[i]
  
}

#relative_abundance=relative_abundance+1e^-6
#relative_abundance2=relative_abundance+2
relative_abundancet=t(relative_abundance)
relative_abundancet[-1,]

relative_abundancet_append=relative_abundancet[-dim(relative_abundancet)[1],]
paste(colnames(relative_abundancet_append),'i',sep="_")

colnames(relative_abundancet_append)=paste(colnames(relative_abundancet_append),'i',sep="_")

relative_abundancet_full=cbind(relative_abundancet_append,relative_abundancet[-1,])

relative_abundancet
relative_abundancet_full=as.data.frame(relative_abundancet_full)
relative_abundancet_full

num_species=dim(relative_abundancet_full)[2]/2

alpa=matrix(NA,dim(relative_abundance)[1],dim(relative_abundance)[1])
for(i in 1:num_species){
  print(i)
  #print(colnames(relative_abundancet_full)[i])
  #print(colnames(relative_abundancet_full)[i+num_species])
  formular=paste(colnames(relative_abundancet_full)[i+num_species],"~",colnames(relative_abundancet_full)[i],"+",colnames(relative_abundancet_full)[i],"*(1 - ",colnames(relative_abundancet_full)[i],sep=" ")
  #startc=c()
  startl=list()
  startc=c()
  for(j in 1:num_species)
    {
  if(i != j){
    #print(j)
    var=paste('a',j,sep='')
    #print(var)
    #assign(var,1)
    startc=append(startc,var)
    #startc=append(startc,starte)
    startl=append(startl,1)
    
    formular=paste(formular,"-",var,"*",colnames(relative_abundancet_full)[j])
          }    
    }
  formular=paste(formular,") * 0.5 + 1e-6", sep=' ')
  names(startl)=startc
  print(formular)
ps=summary(nls(formular,data=relative_abundancet_full,start=startl))
#print(ps)
#print(ps$coefficients[,1])
inter=ps$coefficients[,1]
inter=append(inter,1,i-1)
print(inter)
alpa[i,]=inter
}


#返回的结果是相关互作系数
#The returned result is the correlation interaction coefficient
print(alpa)
#
abundance
relative_abundance[,1]
alpa
specisNum=3
generation=8

relative_abundance
predict_abundance=function(alpa,relative_abundance,specisNum,generation){
  
  alpam=alpa
  #set day 0 the abundance of each bacteria
  #species is 3
  S=specisNum
  #geneartion or days is 8 days
  g=generation
  step=0.5
  # 根据本文假设，所有物种的增长率相同，我们设置为1，即所有每一代在没有其他因素的影响下增长一倍
  r=matrix(1,S, 1)
  # AA为物种间相互作用强度矩阵，原文写的其取值范围为U[0,2A],所以我们在此分布上进行取样
  AA=alpam
  
  # 生存矩阵NN，记录结果
  NN=matrix(0, g, S)
  
  NN[1,]=relative_abundance[,1]
  #设置一个统一的迁入率
  D=1e-6
  Migra=matrix(1,1,S)*D
  
  for (i in 2:g){
    #i=2
    for (j in 1:S){
      # 注意下式
      k=r[j] * NN[i-1,j] * (1-AA[j,] %*% NN[i-1,]) * step
      NN[i,j]=NN[i-1,j]+k
      # 如果物种丰度为负值，我们就将其赋为0 
      if(NN[i,j] < 0){
        NN[i,j] = 0
      }
    }
    NN[i,]=NN[i,]+Migra
  }
  
  #t(NN)
  par(mfrow=c(2,1)) 
  
  matplot(seq_len(dim(NN)[1]), NN, type='l', ylim = c(0,1),lwd=2,
          xlab = "Time", ylab = expression(Predict~Species~Abundance~N[i]),main="Predict")
  legend("topright", colnames( t(relative_abundance) ),col=seq_len( dim( NN )[2] ), lty=seq_len(dim( NN )[2]), lwd=2  )
  matplot(seq_len(dim( t(relative_abundance)[1:generation,] )[1]),t(relative_abundance)[1:generation,],type="l",ylim=c(0,1),lwd=2,
          xlab = "Time", ylab = expression(Real~Species~Abundance~N[i]),main="Real")
  legend("topright", colnames( t(relative_abundance) ),col=seq_len( dim( t(relative_abundance) )[2] ), lty=seq_len(dim( t(relative_abundance) )[2]), lwd=2  )
  
  
  return(NN)
}
print(alpa)
predict_abundance(alpa,relative_abundance,3,8)

rm(list=ls())








