generate_p.sample.mat<-function(ndose.A,ndose.B, NN=50000, target, epi)
{
    #### This is a function for generating samples of the toxicity matrix p from the joint prior distribution.
  ll<-0.00
  uu<-0.80
  p.sample.mat<-array(NA,dim=c(NN,ndose.A,ndose.B,ndose.A,ndose.B)) 
  for(j in 1:ndose.A)
  {
    for(k in 1:ndose.B)
    {
      # Under model M_{j,k}.
      
      # Sample p_{j,k}\sim \mathrm{Unif}(\phi-\epsilon, \phi+\epsilon).  
      p.sample.mat[1:NN,j,k,j,k]<-runif(NN,target-epi,target+epi)
      
      #### For k'=k. ####
      if (j>1)
      {
        for (l_j in (j-1):1)
        {
         	p.sample.mat[1:NN,l_j,k,j,k]<-sapply(1:NN,function(nn) 
            runif(1,ll,min(p.sample.mat[nn,l_j+1,k,j,k],target-epi))
          )
        }
      }
      
      if (j<ndose.A)
      {
        for (l_j in (j+1):ndose.A)
        {
         p.sample.mat[1:NN,l_j,k,j,k]<-sapply(1:NN,function(nn)
            runif(1,max(p.sample.mat[nn,l_j-1,k,j,k],target+epi),uu)
          )
        }
      }
      
      #### j'=j ####
      if (k>1)
      {
        for (l_k in (k-1):1)
        {
        	p.sample.mat[1:NN,j,l_k,j,k]<-sapply(1:NN,function(nn) runif(1,ll,min(p.sample.mat[nn,j,l_k+1,j,k],target-epi))
          )
        }
      }
      
      if (k<ndose.B)
      {
        for (l_k in (k+1):ndose.B)
        {
         	p.sample.mat[1:NN,j,l_k,j,k]<-sapply(1:NN,function(nn) runif(1,max(p.sample.mat[nn,j,l_k-1,j,k],target+epi),uu)
          )
        }
      }
      
      # For elements not on j-th row and not on k-th column. 
      if (j>1)
      {
        for (l_j in (j-1):1)
        {
          if (k>1)
          {
            for (l_k in (k-1):1)
            {
              ub <- sapply(1:NN, function(nn) min(p.sample.mat[nn,l_j+1, l_k, j, k], p.sample.mat[nn,l_j, l_k+1, j, k]))              
              p.sample.mat[1:NN,l_j,l_k,j,k]<-sapply(1:NN,function(nn) runif(1,ll,ub[nn]))
            }
          }
          
          if (k<ndose.B)
          {
            for (l_k in (k+1):ndose.B)
            {
              lb<-sapply(1:NN, function(nn) p.sample.mat[nn,l_j,l_k-1,j, k])
              ub<-sapply(1:NN, function(nn) p.sample.mat[nn,l_j+1,l_k,j, k])

              p.sample.mat[1:NN,l_j,l_k,j,k] <- runif(NN,lb,ub)
              vector<-sapply(1:NN, function(nn) return((lb[nn]<target-epi)&&(ub[nn]>target+epi)))
              index<-which(vector==1)
              
              if (length(index)>0)
              {
              	for (i in 1:length(index))
              	{
               		prob<-c(target-epi-lb[index[i]],ub[index[i]]-target-epi)
                	indicator<-sample(x=c(-1,1),size=1,replace=T,prob=prob/sum(prob))
                	if(indicator==-1)
                	{p.sample.mat[index[i],l_j,l_k,j,k]<-runif(1,lb[index[i]],target-epi)}
                	else
                	{p.sample.mat[index[i],l_j,l_k,j,k]<-runif(1,target+epi,ub[index[i]])}
              	}	
              }	
            }
          }
        }
      }
      
      if (j<ndose.A)
      {
        for (l_j in (j+1):ndose.A)
        {
          if (k>1)
          {
            for (l_k in (k-1):1)
            {
              lb<-sapply(1:NN,function(nn) p.sample.mat[nn,l_j-1,l_k,j,k])
              ub<-sapply(1:NN,function(nn) p.sample.mat[nn,l_j,l_k+1,j,k])
              
              p.sample.mat[1:NN,l_j,l_k,j,k] <- runif(NN,lb,ub)
              vector<-(sapply(1:NN, function(nn) return((lb[nn]<target-epi)&&(ub[nn]>target+epi))))
              
              index<-which(vector==1)
              
              if (length(index)>0)
              {
              	for (i in 1:length(index))
              	{
               		prob<-c(target-epi-lb[index[i]],ub[index[i]]-target-epi)
                	indicator<-sample(x=c(-1,1),size=1,replace=T, prob=prob/sum(prob))
                	if(indicator==-1)
                	{p.sample.mat[index[i],l_j,l_k,j,k]<-runif(1,lb[index[i]],target-epi)}
                	else
                	{p.sample.mat[index[i],l_j,l_k,j,k]<-runif(1,target+epi,ub[index[i]])}
              	}
              }
            }
          }
          
          if (k<ndose.B)
          {
            for (l_k in (k+1):ndose.B)
            {
              lb<-sapply(1:NN, function(nn) max(p.sample.mat[nn,l_j-1,l_k, j, k],p.sample.mat[nn,l_j,l_k-1,j,k]))
              # p.sample.mat[1:NN,l_j,l_k,j,k]<-sapply(1:NN,function(nn) runif(1,lb[nn],0.8))	
              p.sample.mat[1:NN,l_j,l_k,j,k]<-sapply(1:NN,function(nn) runif(1,lb[nn],uu))							
            }
          }
        }
      }
    }
  }
  return(p.sample.mat=p.sample.mat)
}
