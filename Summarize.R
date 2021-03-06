# A function for summarizing the simulation results. 
summarize<-function(sim_Results,nsim,target)
{
    for (simid in 1:nsim)
    {
        sim_Results[[simid]]$y.mat<-sim_Results[[simid]]$y.mat[nrow(sim_Results[[simid]]$y.mat):1,]
        sim_Results[[simid]]$n.mat<-sim_Results[[simid]]$n.mat[nrow(sim_Results[[simid]]$n.mat):1,]
        sim_Results[[simid]]$pos.model<-sim_Results[[simid]]$pos.model[nrow(sim_Results[[simid]]$pos.model):1,]
        sim_Results[[simid]]$pos.model.original<-sim_Results[[simid]]$pos.model.original[nrow(sim_Results[[simid]]$pos.model.original):1,]
        sim_Results[[simid]]$Tox_Prob_Mat<-sim_Results[[simid]]$Tox_Prob_Mat[nrow(sim_Results[[simid]]$Tox_Prob_Mat):1,]
    }
    
	summary.MTD<-array(0,dim=dim(sim_Results[[1]]$Tox_Prob_Mat))
	summary.dlt<-array(0,dim=dim(sim_Results[[1]]$Tox_Prob_Mat))
	summary.patients<-array(0,dim=dim(sim_Results[[1]]$Tox_Prob_Mat))
	
	nPatient<-sum(sim_Results[[1]]$n.mat)
	ind.MTD<-which(sim_Results[[1]]$Tox_Prob_Mat==target,arr.ind=T)
	ind.overtox<-which(sim_Results[[1]]$Tox_Prob_Mat>target,arr.ind=T)
	
	for (simid in 1:nsim)
	{
		summary.MTD[sim_Results[[simid]]$MTD.comb]<-summary.MTD[sim_Results[[simid]]$MTD.comb]+1
		summary.dlt<-summary.dlt+sim_Results[[simid]]$y.mat
		summary.patients<-summary.patients+sim_Results[[simid]]$n.mat
	}
    
	summary.MTD.pctg<-summary.MTD/nsim*100
	summary.patients.pctg<-summary.patients/nsim/nPatient*100
    summary.dlt.pctg<-summary.dlt/nsim/nPatient*100
	
	select.pctg.correct<-sum(summary.MTD.pctg[ind.MTD]) # Correct selection %.
	allocation.pctg.correct<-sum(summary.patients.pctg[ind.MTD]) # Correct allocation %. 
	select.pctg.overtox<-sum(summary.MTD.pctg[ind.overtox]) # Overdose selection %.
	allocation.pctg.overtox<-sum(summary.patients.pctg[ind.overtox]) # Overdose allocation %. 
	AI<-(1-sum(select.pctg.correct/100*abs(sim_Results[[1]]$Tox_Prob_Mat-target))/sum(abs(sim_Results[[1]]$Tox_Prob_Mat-target))/100*ncol(sim_Results[[1]]$Tox_Prob_Mat)*nrow(sim_Results[[1]]$Tox_Prob_Mat))*100
	DLT<-sum(summary.dlt.pctg)

    summary.MTD.pctg<-summary.MTD.pctg[nrow(summary.MTD.pctg):1,]
    summary.patients.pctg<-summary.patients.pctg[nrow(summary.patients.pctg):1,]
    summary.dlt.pctg<-summary.dlt.pctg[nrow(summary.dlt.pctg):1,]
    Tox_Prob_Mat<-sim_Results[[1]]$Tox_Prob_Mat[nrow(sim_Results[[1]]$Tox_Prob_Mat):1,]
    
    summary_stat<-list(MTD_Sel=select.pctg.correct,MTD_Pts=allocation.pctg.correct,AI=AI,Over_Sel=select.pctg.overtox,Over_Pts=allocation.pctg.overtox,DLT=DLT,summary.MTD.pctg=summary.MTD.pctg,summary.patients.pctg=summary.patients.pctg,summary.dlt.pctg=summary.dlt.pctg,Tox_Prob_Mat=Tox_Prob_Mat)
	return(summary_stat)
}
