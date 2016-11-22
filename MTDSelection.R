source("PosteriorProbability.R")

MTDSelection<-function(y=y,n=n,target=target,p.sample.mat=p.sample.mat)
{
    y<-y[nrow(y):1,]
    n<-n[nrow(n):1,]
	pos.model<-posteriorH(y=y,n=n,target=target,p.sample.mat=p.sample.mat)$lik
	### MTD Selection. ####
    pos.model.original<-pos.model[nrow(pos.model):1,]
    pos.model<-pos.model*(n>0)
	MTD.sel<-which(pos.model==max(pos.model),arr.ind=TRUE)

	rownames(pos.model.original)<-paste0("Dose", c(nrow(n):1))
	colnames(pos.model.original)<-paste0("Dose", c(1:ncol(n)))
	
	return(list(pos.model=pos.model.original,MTD.sel=MTD.sel))
}
