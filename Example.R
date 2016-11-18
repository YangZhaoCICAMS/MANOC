rm(list=ls())

setwd("/Users/lamchikin/Dropbox/NOC_DrugComb/MANOC_master/")
source("NextDoseComb.R")
source("PosteriorProbability.R")
source("Simulation.R")
source("ToxProb_Generate.R")
source("Summarize.R")

require(parallel)

# A toxicity probability matrix. 

Tox_Prob_Mat <- 
matrix(c(
0.08,0.10,0.15,0.30,
0.14,0.20,0.30,0.42,
0.19,0.30,0.45,0.60,
0.30,0.47,0.60,0.70
),
nrow=4, ncol=4
)

samplesize <- 60
cohortsize <- 3
target <- 0.30
epi <- 0.025
NN <- 1000

## Generate the matrices of p according to the prior distribution of each model. ## 
p.sample.mat <- generate_p.sample.mat(ndose.A=nrow(Tox_Prob_Mat),ndose.B=ncol(Tox_Prob_Mat), NN=NN, target=target, epi=epi) 

alpha<-0.35
delta<-0.03
eta<-0.60

nsim <- 100

sim_Results<-mclapply(1:nsim, function(simid) simulation(simid=simid,Tox_Prob_Mat=Tox_Prob_Mat, p.sample.mat=p.sample.mat,samplesize=samplesize, cohortsize=cohortsize, target=target, alpha=alpha, delta=delta, eta=eta), mc.cores=1)

summarize(sim_Results=sim_Results,nsim=nsim)