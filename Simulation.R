# The function for conducting simulation study. # 

simulation <- function(simid, Tox_Prob_Mat, p.sample.mat, samplesize=60, cohortsize=3, target=0.3, alpha=0.35, delta=0.03, eta=0.60)
{
	set.seed(simid)
	cohort <- samplesize/cohortsize
	dlt <- rep(0, cohort)
	dose.level <-array(0,dim=c(cohort+1,2))
	dose.level[1, ] <- c(1,1)

	ndose.A <- nrow(Tox_Prob_Mat)
	ndose.B <- ncol(Tox_Prob_Mat)
	y.mat <- array(0, dim=c(ndose.A, ndose.B))
	n.mat <- array(0, dim=c(ndose.A, ndose.B))

	for (i in 1:cohort)
	{
    	index <- dose.level[i,]
		dlt[i] <- rbinom(1, cohortsize, Tox_Prob_Mat[index[1], index[2]]) 
				
		y.mat[index[1],index[2]] <- y.mat[index[1],index[2]] + dlt[i]
		n.mat[index[1],index[2]] <- n.mat[index[1],index[2]] + cohortsize 
		pos.model<-posteriorH(y=y.mat,n=n.mat,target=target,p.sample.mat=p.sample.mat)$lik
		
		pos.model.modified<-pos.model
		pos.model.modified[which(n.mat==0,arr.ind=TRUE)]<-pos.model.modified[which(n.mat==0,arr.ind=TRUE)]+delta
		
		if (sum(dlt)==0)
		{
			dose.level[i+1,]<-c(min(ndose.A,index[1]+1),min(ndose.B,index[2]+1))
		} 
		else	
		{
			dose.level[i+1,] <- get.next.manoc(pos.model=pos.model.modified,j_curr=index[1],k_curr=index[2],alpha=alpha,eta=eta)$next.dose
		}
	}
	
	### MTD Selection. ####
	pos.model.original<-pos.model
	pos.model<-pos.model*(n.mat>0)
	MTD.comb<-which(pos.model==max(pos.model),arr.ind=TRUE)
				
return(list(MTD.sel=Tox_Prob_Mat[MTD.comb],MTD.comb=MTD.comb,y.mat=y.mat,n.mat=n.mat,dose.level.explored=dose.level[1:cohort,],pos.model=pos.model,pos.model.original=pos.model.original,Tox_Prob_Mat=Tox_Prob_Mat))
}