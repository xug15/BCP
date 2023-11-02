# AA为物种间相互作用强度矩阵
alpamt=matrix(c(1,0.4972781,-7.564028,
                1.649809,1,15.372411,
                2.982275, 3.259204,1
),3,3)
alpam=t(alpamt)
alpam
#set day 0 the abundance of each bacteria
d0abundance=c(0.688642126,0.248462592,0.062895281)
#species is 3
S=3
#geneartion or days is 8 days
g=8
  
  # 根据本文假设，所有物种的增长率相同，我们设置为1，即所有每一代在没有其他因素的影响下增长一倍
  r=matrix(1,S, 1)
  # AA为物种间相互作用强度矩阵，原文写的其取值范围为U[0,2A],所以我们在此分布上进行取样
  AA=alpam

  # 生存矩阵NN，记录结果
  NN=matrix(0, g, S)

  NN[1,]=d0abundance
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

t(NN)
matplot(seq_len(dim(NN)[1]), NN, type='l', ylim = c(0,1),
        xlab = "Time", ylab = expression(Species~Abundance~N[i]))









