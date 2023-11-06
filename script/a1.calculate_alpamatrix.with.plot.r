#!/usr/bin/env Rscript
#v1.6 to 1.7
print("Calinter.R abundance.csv output.csv")
rm(list=ls())

#read bacteria abundance file
# abundance=read.csv('d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/r_code/b1.b3.bacteria.abundance.csv',row.names = 1)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}
abundance=read.csv(args[1],row.names=1)
head(abundance)

##########body being
#pass the data.
generate_relatvie_dataframe=function(abundance){
  ## genearate relative abundance dataframe begin.
  relative_abundance=abundance
  #2. Calculate the percentage.
  for (i in 1:dim(abundance)[2]){
    #print(i)
    relative_abundance[,i]=abundance[,i]/colSums(abundance)[i]
    
  }
  
  #relative_abundance=relative_abundance+1e^-6
  #relative_abundance2=relative_abundance+2
  #relative_abundance=relative_abundance*10
  dim(relative_abundance)
  
  
  relative_abundancet=t(relative_abundance)
  relative_abundancet[-1,]
  #relative_abundancet=relative_abundancet*10
  
  
  
  
  relative_abundancet_append=relative_abundancet[-dim(relative_abundancet)[1],]
  #paste(colnames(relative_abundancet_append),'i',sep="_")
  
  colnames(relative_abundancet_append)=paste(colnames(relative_abundancet_append),'i',sep="_")
  
  relative_abundancet_full=cbind(relative_abundancet_append,relative_abundancet[-1,])
  
  relative_abundancet_full=as.data.frame(relative_abundancet_full)
  return(list(relative_abundancet_full,relative_abundance))
  
  #relative_abundancet_full=relative_abundancet_full*10
  ## genearate relative abundance dataframe end.
}



calculate_first_interaction=function(relative_abundance,relative_abundancet_full){
  

    # First around calculate interaction matrix. begin.
    num_species=dim(relative_abundancet_full)[2]/2
    alpa=matrix(NA,dim(relative_abundance)[1],dim(relative_abundance)[1])
    BB=matrix(runif(dim(relative_abundancet_full)[1]*dim(relative_abundancet_full)[2],min=0,max=1), nrow=dim(relative_abundancet_full)[1]) * 1e-8
    CC=matrix(runif(dim(relative_abundancet_full)[1]*dim(relative_abundancet_full)[2],1), nrow=dim(relative_abundancet_full)[1]) * 1e-8
    
    #first run 
    
    maxc=c()
    for(i in 1:num_species){
      #i=2
      #i=1
      maxv=0
      #print(i)
      #print(colnames(relative_abundancet_full)[i])
      #print(colnames(relative_abundancet_full)[i+num_species])
      formular=paste(colnames(relative_abundancet_full)[i+num_species],"~",colnames(relative_abundancet_full)[i],"+",colnames(relative_abundancet_full)[i],"*(1 - ",colnames(relative_abundancet_full)[i],sep=" ")
      #startc=c()
      startl=list()
      startc=c()
      AA=matrix(runif(num_species^2,min=0,max=1), nrow=num_species) *2 * 0.8
      AA=AA-diag(diag(AA))+diag(num_species)
      #
      
      for(j in 1:num_species)
      {
        if(i != j){
          #print(j)
          var=paste('a',j,sep='')
          #print(var)
          #assign(var,1)
          startc=append(startc,var)
          startl=append(startl,AA[i,j])
          #print(c(i,j))
          #print(AA[i,j])
          
          formular=paste(formular,"-",var,"*",colnames(relative_abundancet_full)[j])
        }    
      }
      formular=paste(formular,") * 0.5 + 1e-6", sep=' ')
      names(startl)=startc
      #print(formular)
      
      sig=tryCatch(nls(formular,data=relative_abundancet_full,start=startl,nls.control(warnOnly = T,nDcentral = T)),error=function(e) return(FALSE))
      #
      #
      if(length(sig)==1){
        #### function to test
        iter=1
        while(length(sig)==1){
          
          if(length(sig)==1){
            #print(iter)
            #BB=matrix(runif(num_species*dim(relative_abundancet_full)[1],min=0,max=1), nrow=dim(relative_abundancet_full)[1]) *2 * 1e-6
            relative_abundancet_full2=relative_abundancet_full+BB*2^iter
            sig=tryCatch(nls(formular,data=relative_abundancet_full2,start=startl,nls.control(warnOnly = T)),error=function(e) return(FALSE))
            iter=iter+1
          }
          
          #print(iter)
          if(iter > 2000){
            break
          }
          
        }
        ####
        #print(paste("change the origin data add 2", iter,sep="^"))
        maxv=iter
      }
      
      #length(sig)
      ps=summary(sig)
      
      inter=ps$coefficients[,1]
      inter=append(inter,1,i-1)
      #print(inter)
      alpa[i,]=inter
      maxc=append(maxc,maxv)
    }
    return(list(alpa,maxc,BB))
    # First around calculate interaction matrix. end.
    
}   



calculate_second_interaction=function(relative_abundance,relative_abundancet_full,alpa,maxc,BB){
  
    num_species=dim(relative_abundancet_full)[2]/2
    ########## begin.
    #second run and use optimal initail paramters.
    maxiter=max(maxc)-1
    
    
    alpa2=matrix(NA,dim(relative_abundance)[1],dim(relative_abundance)[1])
    for(i in 1:num_species){
      #i=2
      #print(i)
      #print(colnames(relative_abundancet_full)[i])
      #print(colnames(relative_abundancet_full)[i+num_species])
      formular=paste(colnames(relative_abundancet_full)[i+num_species],"~",colnames(relative_abundancet_full)[i],"+",colnames(relative_abundancet_full)[i],"*(1 - ",colnames(relative_abundancet_full)[i],sep=" ")
      #startc=c()
      startl=list()
      startc=c()
      #
      
      for(j in 1:num_species)
      {
        if(i != j){
          #print(j)
          var=paste('a',j,sep='')
          #print(var)
          #assign(var,1)
          startc=append(startc,var)
          startl=append(startl,alpa[i,j])
          #print(c(i,j))
          #print(AA[i,j])
          
          formular=paste(formular,"-",var,"*",colnames(relative_abundancet_full)[j])
        }    
      }
      formular=paste(formular,") * 0.5 + 1e-6", sep=' ')
      names(startl)=startc
      #print(formular)
      
      #sig=tryCatch(nls(formular,data=relative_abundancet_full,start=startl,nls.control(warnOnly = T)),error=function(e) return(FALSE))
      relative_abundancet_full2=relative_abundancet_full+BB*2^maxiter
      sig=tryCatch(nls(formular,data=relative_abundancet_full2,start=startl,nls.control(warnOnly = T)),error=function(e) return(FALSE))
      #
      #
      if(length(sig)==1){
        #### function to test
        iter=1
        while(length(sig)==1){
          
          if(length(sig)==1){
            #print(iter)
            #BB=matrix(runif(num_species*dim(relative_abundancet_full)[1],min=0,max=1), nrow=dim(relative_abundancet_full)[1]) *2 * 1e-6
            relative_abundancet_full2=relative_abundancet_full+BB*2^iter
            sig=tryCatch(nls(formular,data=relative_abundancet_full2,start=startl,nls.control(warnOnly = T)),error=function(e) return(FALSE))
            iter=iter+1
          }
          
          #print(iter)
          if(iter > 2000){
            break
          }
          
        }
        ####
        print(paste("second change the origin data add 2", iter,sep="^"))
      }
      
      #length(sig)
      ps=summary(sig)
      
      inter=ps$coefficients[,1]
      inter=append(inter,1,i-1)
      #print(inter)
      alpa2[i,]=inter
    }
    #############end
    return(alpa2)

}
   
    print("the result process here.")
    
    #返回的结果是相关互作系数
    #The returned result is the correlation interaction coefficient
    predict_abundance=function(alpa,relative_abundance,specisNum){
      
      alpam=alpa
      #set day 0 the abundance of each bacteria
      #species is 3
      S=specisNum
      #geneartion or days is 8 days
      g=dim(relative_abundance)[2]
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
        
        #print(paste("i",i,sep=":"))
        #i=2
        for (j in 1:S){
          #print(paste("j",j,sep=":"))
          # 注意下式
          k=r[j] * NN[i-1,j] * (1-AA[j,] %*% NN[i-1,]) * step
          NN[i,j]=NN[i-1,j]+k
          # 如果物种丰度为负值，我们就将其赋为0 
          #print(c(i,j))
          if(NN[i,j] < 0 | is.na(NN[i,j])){
            NN[i,j] = 0
          }
        }
        NN[i,]=NN[i,]+Migra
      }
      
      sum=0
      for(ii in 1:dim(NN)[1]){
        #print('relative')
        #print(relative_abundance[,ii])
        #print('NN')
        #print( NN[ii,])
        for(jj in 1:dim(NN)[2]){
          sum=sum+(relative_abundance[jj,ii]-NN[ii,jj])^2
        }
      }
      #print(sum)
      return(sum)
    }
    
    
    
    
    
    
    
    differencev=c()
    interactionlist=list()
    
    step=1
    
    
    
    generate_abundence_list=generate_relatvie_dataframe(abundance)
    relative_abundancet_full=generate_abundence_list[[1]]
    relative_abundance=generate_abundence_list[[2]]
    
    
    
    interaction_1_list=calculate_first_interaction(relative_abundance,relative_abundancet_full)
    alpa=interaction_1_list[[1]]
    maxc=interaction_1_list[[2]]
    BB =interaction_1_list[[3]]
    alpa2=calculate_second_interaction(relative_abundance,relative_abundancet_full,alpa,maxc,BB)
    
    
    #print(alpa)
    predictv=predict_abundance(alpa,relative_abundance,dim(alpa)[1])
    differencev=append(differencev,predictv)
    interactionlist[[step*2-1]]=alpa
    
    #print(alpa2)
    predictv=predict_abundance(alpa2,relative_abundance,dim(alpa)[1])
    differencev=append(differencev,predictv)
    interactionlist[[step*2]]=alpa2
    #########body end
    step=step+1

    loop_step=function(loopnum){
      generate_abundence_list=generate_relatvie_dataframe(abundance)
      relative_abundancet_full=generate_abundence_list[[1]]
      relative_abundance=generate_abundence_list[[2]]
      #
      differencev=c()
      interactionlist=list()
      step=1
      for (i in 1:loopnum){
        print(i)
        interaction_1_list=calculate_first_interaction(relative_abundance,relative_abundancet_full)
        alpa=interaction_1_list[[1]]
        maxc=interaction_1_list[[2]]
        BB =interaction_1_list[[3]]
        alpa2=calculate_second_interaction(relative_abundance,relative_abundancet_full,alpa,maxc,BB)
        
        
        #print(alpa)
        predictv=predict_abundance(alpa,relative_abundance,dim(alpa)[1])
        differencev=append(differencev,predictv)
        interactionlist[[step*2-1]]=alpa
        
        #print(alpa2)
        predictv=predict_abundance(alpa2,relative_abundance,dim(alpa)[1])
        differencev=append(differencev,predictv)
        interactionlist[[step*2]]=alpa2
        
        #
        step=step+1
        
      } 
      
      
      minindex=which.min(differencev)
      interactionlist[[2]][minindex]
      #return(list(interactionlist[[minindex]],min(differencev)))
      return(interactionlist[[minindex]])
  }
#
    alpa=loop_step(100)
    
#alpa=alpa_list[[1]]

row.names(relative_abundance)
row.names(alpa)=row.names(relative_abundance)
colnames(alpa)=row.names(relative_abundance)
alpa
print(alpa)
##############
predict_abundance2=function(alpa,relative_abundance,specisNum){
  
  alpam=alpa
  #set day 0 the abundance of each bacteria
  #species is 3
  S=specisNum
  #geneartion or days is 8 days
  g=dim(relative_abundance)[2]
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
  
  matplot(seq_len(dim(NN)[1]), NN, type='l',lwd=2,
          xlab = "Time", ylab = expression(Predict~Species~Abundance~N[i]),main="Predict")
  legend("topright", colnames( t(relative_abundance) ),col=seq_len( dim( NN )[2] ), lty=seq_len(dim( NN )[2]), lwd=2  )
  matplot(seq_len(dim( t(relative_abundance))[1]),t(relative_abundance),type="l",ylim=c(0,1),lwd=2,
          xlab = "Time", ylab = expression(Real~Species~Abundance~N[i]),main="Real")
  legend("topright", colnames( t(relative_abundance) ),col=seq_len( dim( t(relative_abundance) )[2] ), lty=seq_len(dim( t(relative_abundance) )[2]), lwd=2  )
  sum=0
  for(ii in 1:dim(NN)[1]){
    #print('relative')
    #print(relative_abundance[,ii])
    #print('NN')
    #print( NN[ii,])
    for(jj in 1:dim(NN)[2]){
      sum=sum+(relative_abundance[jj,ii]-NN[ii,jj])^2
    }
  }
  print(sum)
  return(NN)
}

#
##############

predict_abundance2(alpa,relative_abundance,dim(alpa)[1])



    