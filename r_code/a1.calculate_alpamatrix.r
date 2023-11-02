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
relative_abundance

# Calculate the difference with time.
## set begin matrix
growth=matrix(data=NA,dim(relative_abundance)[1],1)

#the calculate the difference between adjust times data. 
for( i in 2:dim(relative_abundance)[2]){
  print(i)
  #i=2
  df=as.matrix(relative_abundance[i]-relative_abundance[i-1])
  growth=cbind(growth,df)
}
growth=growth[,-1]
growth
#The calculate the determinant multiplication production
#计算行列式相乘的结果,这个是方向计算,所以为后面倒退的相乘后的结果.
determinant=matrix(data=NA,dim(relative_abundance)[1],1)
determinant
for (i in 1:dim(growth)[2]){
  print(i)
  multiply=as.matrix((1-2*growth[,i]/relative_abundance[,i])-relative_abundance[,i])
  determinant=cbind(determinant,multiply)
}

determinant=determinant[,-1]
determinant
#remove na column.
determinant2=determinant[,colSums(is.na(determinant))==0]
# the determinant product = 
# for example s1=s2*alpa21+s3*alpa31
#use column
#选择使用的列,这里选择1,2列.
cal_colum=c(1,2,3,4,5)

alpa=matrix(NA,dim(relative_abundance)[1],dim(relative_abundance)[1])
#By finding the inverse of the determinant, the relevant coefficients can be calculated,
#求行列式的逆,可以计算出相关的系数,
for (i in 1:dim(relative_abundance)[1]){
  print(i)
  #i=1
  print(determinant[i,cal_colum])
  product=unlist(determinant[i,cal_colum])
  species_determinant=c()
  for (j in 1:dim(relative_abundance)[1]){
    if(i ==j){
      next;
    }else{
      #print(j)
      #print(relative_abundance[j,cal_colum])
      species_determinant=append(species_determinant,unlist(relative_abundance[j,cal_colum]))
    }
    
  }
  length(species_determinant)/5
  spcecies_matrix=matrix(species_determinant,length(species_determinant)/length(cal_colum),length(cal_colum))
  alpa[i,]=append(solve(spcecies_matrix,product),1,i-1)
  
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

alpa2=alpa
alpa2[1,2]= 0.5
print(alpa2)
predict_abundance(alpa2,relative_abundance,3,8)



