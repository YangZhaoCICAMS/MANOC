posteriorH<-function(y,n,target,p.sample.mat)
{
  ndose.A <- nrow(y)
	ndose.B <- ncol(y)
	lik<-array(NaN, dim=dim(y))
	for (j in 1:ndose.A)
	{
		for (k in 1:ndose.B)
		{
			likeli<-1
			for (j0 in 1:ndose.A)
			{
				for (k0 in 1:ndose.B)
				{
					likeli<-likeli*p.sample.mat[,j0,k0,j,k]^y[j0,k0]*(1-p.sample.mat[,j0,k0,j,k])^(n[j0,k0]-y[j0,k0])
				}
			}
			lik[j,k] <- mean(likeli)
		}
	}

	lik <- lik/sum(lik)
	return(list(lik=lik))
}
