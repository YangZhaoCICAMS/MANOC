# A function for summarizing the simulation results. 
summarize<-function(sim_Results,nsim)
{
	nPatient<-sum(sim_Results[[1]]$n.mat)
	print(nPatient)
	ind.MTD<-which(sim_Results[[1]]$Tox_Prob_Mat==target,arr.ind=T)
	ind.overtox<-which(sim_Results[[1]]$Tox_Prob_Mat>target,arr.ind=T)
	summary.MTD<-array(0,dim=dim(sim_Results[[1]]$Tox_Prob_Mat))
	summary.dlt<-array(0,dim=dim(sim_Results[[1]]$Tox_Prob_Mat))
	summary.patients<-array(0,dim=dim(sim_Results[[1]]$Tox_Prob_Mat))
	summary.risk.hightox<-rep(NaN,nsim)
	for (simid in 1:nsim)
	{
		summary.MTD[sim_Results[[simid]]$MTD.comb]<-summary.MTD[sim_Results[[simid]]$MTD.comb]+1
		summary.dlt<-summary.dlt+sim_Results[[simid]]$y.mat
		summary.patients<-summary.patients+sim_Results[[simid]]$n.mat
		summary.risk.hightox[simid]<-1*((sum(sim_Results[[simid]]$y.mat))/sum(sim_Results[[simid]]$n.mat)>0.3)
	}
	print("summary.risk.hightox")
	print( summary.risk.hightox )
	print("DLT")
	print(sapply(1:nsim,function(simid) sum(sim_Results[[simid]]$y.mat)/nPatient*100))
	summary.MTD.pctg<-summary.MTD/nsim*100 # Matrix. 
	summary.patients.pctg<-summary.patients/nsim/nPatient*100 # Matrix.
	summary.dlt.pctg<-summary.dlt/nsim/nPatient*100 # Matrix
	summary.risk.hightox.pctg<-mean(summary.risk.hightox)*100
	
	select.pctg.correct<-sum(summary.MTD.pctg[ind.MTD]) # Correct selection %.
	allocation.pctg.correct<-sum(summary.patients.pctg[ind.MTD]) # Correct allocation %. 
	select.pctg.overtox<-sum(summary.MTD.pctg[ind.overtox]) # Overdose selection %.
	allocation.pctg.overtox<-sum(summary.patients.pctg[ind.overtox]) # Overdose allocation %. 
	
	DLT<-sum(summary.dlt.pctg)
	summary_stat<-list(Cor_Sel=select.pctg.correct,Cor_All=allocation.pctg.correct,Over_Sel=select.pctg.overtox,Over_All=allocation.pctg.overtox,Risk_highTox=summary.risk.hightox.pctg,DLT=DLT,summary.MTD.pctg=summary.MTD.pctg,summary.patients.pctg=summary.patients.pctg,summary.dlt.pctg=summary.dlt.pctg,Tox_Prob_Mat=sim_Results[[1]]$Tox_Prob_Mat)
	return(summary_stat)
}