# MANOC (Multi-agent Nonparametric Overdose Control)
R codes to implement the multi-agent nonparametric overdose control design for dose finding in phase I drug-combination trials.

# Inputs 
- samplesize: The maximum number of patients to be enrolled.  
- cohortsize: The number of patients in each cohort. 
- target: The target toxicity rate. 
- epi: A small positive number epsilon in the model specification.  
- alpha: The prespecified feasible bound.    
- delta: A small increment on the posterior probabilities for the untried dose combinations.
- eta: The dose-switching cutoff.
- NN: The number of samples of **p** generated from its prior distribution.
- nsim: The number of trials simulated under each scenario. 
- Tox_Prob_Mat: The prespecified toxicity probability under each scenario. 

# Description
# Functions
- NextDoseComb.R: Containing a function `get.next.manoc()` for determining the next dose combination given the current dose combination, the posterior model probabilities, alpha and eta. 
- PosteriorProbability.R: Containing a function `posteriorH()` for calculating the posterior model probability for each dose combination.
- Simulation.R: Containing a function `simulation()` for conducting simulation studies. 
- Summarize.R: Containing a function `summarize()` for summarizing the outputs produced by the function `simulation()`.
- ToxProb_Generate.R: Containing a function `generate_p.sample.mat()` for generating samples of the toxicity matrix **p** from its prior distribution. Details can be found in the Appendix of the paper. 

# Examples
We apply the MANOC design to the phase Ib trial with a combination of buparlisib and trametinib.

## MTD Selection
- Suppose at the end of trial, the number of patients treated at each dose combination *n* and the corresponding number of toxicities *y* are 
```rscript
> n
     [,1] [,2] [,3] [,4] [,5]
[1,]    3    0    0    0    3
[2,]    0    3    0    9   39
[3,]    0    0    3    3    0
[4,]    0    0    0    3    0
> 
> y
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    0    0    0    0
[2,]    0    0    0    1   12
[3,]    0    0    0    2    0
[4,]    0    0    0    2    0
```
We can use the following code to select the MTD combination. 
```rscript
rm(list=ls())

setwd("/MANOC_master/")
source("ToxProb_Generate.R")
source("PosteriorProbability.R")
source("MTDSelection.R")
target <- 0.33
epi <- 0.025
NN <- 50000

n<-matrix(c(3,0,0,0, 0,3,0,0, 0,0,3,0, 0,9,3,3, 3,39,0,0),4,5)
y<-matrix(c(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,2,2, 0,12,0,0),4,5)

p.sample.mat <- generate_p.sample.mat(ndose.A=nrow(n),ndose.B=ncol(n), NN=NN, target=target, epi=epi) 
MTDSelection(y=y,n=n,target=target,p.sample.mat=p.sample.mat)
```
The output including the posterior probability of each dose combination and the selected MTD pair are respectively given by
```rscript
> MTDSelection(y=y,n=n,target=target,p.sample.mat=p.sample.mat)
$pos.model
     [,1] [,2] [,3] [,4] [,5]
[1,] 0.00 0.00 0.00 0.00 0.02
[2,] 0.00 0.00 0.00 0.02 0.42
[3,] 0.00 0.00 0.04 0.19 0.05
[4,] 0.01 0.05 0.15 0.06 0.00

$MTD.sel
     row col
[1,]   2   5
```

## Next Dose Level
- If ten cohorts of patients have been enrolled and the corresponding *n* and *y* are
```rscript
> n
     [,1] [,2] [,3] [,4] [,5]
[1,]    3    0    0    0    3
[2,]    0    3    0    9    3
[3,]    0    0    3    3    0
[4,]    0    0    0    3    0
> y
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    0    0    0    0
[2,]    0    0    0    1    2
[3,]    0    0    0    2    0
[4,]    0    0    0    2    0
```
To decide the dose combination at which the eleventh cohort of patients treated, we can use the function `get.next.manoc()`.
```rscript
> rm(list=ls())
> setwd("/MANOC_master/")
> source("ToxProb_Generate.R")
> source("NextDoseComb.R")
> 
> target <- 0.33
> epi <- 0.025
> delta <- 0.05
> alpha<-0.35
> eta<-0.55
> NN <- 50000
> 
> j_curr<-1
> k_curr<-5
>
> n<-matrix(c(3,0,0,0, 0,3,0,0, 0,0,3,0, 0,9,3,3, 3,3,0,0),4,5)
> y<-matrix(c(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,2,2, 0,2,0,0),4,5)
> 
> p.sample.mat <- generate_p.sample.mat(ndose.A=nrow(n),ndose.B=ncol(n), NN=NN, target=target, epi=epi) 
> 
> get.next.manoc(y=y,n=n,target=target,delta=delta,p.sample.mat=p.sample.mat,j_curr=j_curr,k_curr=k_curr,alpha=alpha,eta=eta)
$next.dose
[1] 2 5
```
That is, we assign the eleventh cohort of patients to dose combination (2,5).

## Simulation Studies
```rscript
> rm(list=ls())
> 
> setwd("/MANOC_master/")
> source("NextDoseComb.R")
> source("PosteriorProbability.R")
> source("Simulation.R")
> source("ToxProb_Generate.R")
> source("Summarize.R")
> 
> require(parallel)
> 
> # A toxicity probability matrix. 
> Tox_Prob_Mat <-
+ matrix(c(
+ 0.01,0.01,0.03,0.06,
+ 0.03,0.07,0.15,0.33,
+ 0.06,0.12,0.24,0.42,
+ 0.10,0.21,0.38,0.58,
+ 0.18,0.33,0.53,0.72
+ ),
+ nrow=4, ncol=5
+ )
> 
> samplesize <- 66
> cohortsize <- 3
> target <- 0.33
> epi <- 0.025
> NN <- 50000
> 
> alpha<-0.35
> delta<-0.05
> eta<-0.55
> 
> nsim <- 1000
> 
> ## Generate the matrices of p according to the prior distribution of each model. ## 
> p.sample.mat <- generate_p.sample.mat(ndose.A=nrow(Tox_Prob_Mat),ndose.B=ncol(Tox_Prob_Mat), NN=NN, target=target, epi=epi) 
> 
> sim_Results<-mclapply(1:nsim, function(simid) simulation(simid=simid,Tox_Prob_Mat=Tox_Prob_Mat, p.sample.mat=p.sample.mat,samplesize=samplesize, cohortsize=cohortsize, target=target, alpha=alpha, delta=delta, eta=eta), mc.cores=1)
```

We can use the function `summarize()` to operating characteristics of the MANOC design as a list, including 
```rscript
> summarize(sim_Results=sim_Results,nsim=nsim,target=target)
$Cor_Sel
[1] 39.8

$Cor_All
[1] 22.2

$AI
[1] 70.6

$Over_Sel
[1] 27.6

$Over_All
[1] 22.0

$DLT
[1] 25.1

$summary.MTD.pctg
     [,1] [,2] [,3] [,4] [,5]
[1,]  0.0  0.0  0.0  0.4  4.6
[2,]  0.0  0.0  0.3 13.9 19.1
[3,]  0.0  0.7 12.5 18.0  0.3
[4,]  0.2 20.7  8.9  0.4  0.0

$summary.patients.pctg
     [,1] [,2] [,3] [,4] [,5]
[1,] 4.73 0.09 0.17 1.11 3.92
[2,] 0.20 6.70 3.85 13.0 9.98
[3,] 0.29 3.64 16.6 13.1 0.86
[4,] 1.55 12.2 4.95 2.91 0.14

$summary.dlt.pctg
     [,1] [,2] [,3] [,4] [,5]
[1,] 0.06 0.00 0.01 0.11 0.76
[2,] 0.00 0.49 0.50 2.80 2.98
[3,] 0.01 0.52 3.86 4.93 0.47
[4,] 0.08 3.59 2.07 1.69 0.11

$Tox_Prob_Mat
     [,1] [,2] [,3] [,4] [,5]
[1,] 0.01 0.03 0.06 0.10 0.18
[2,] 0.01 0.07 0.12 0.21 0.30
[3,] 0.03 0.15 0.24 0.38 0.53
[4,] 0.06 0.30 0.42 0.58 0.72
```
Here, the values returned by `summarize()` are
- Cor_Sel: the total correct selection percentage of the MTD combination(s)
- Cor_All: the total percentage of patients treated at the MTD combination(s)
- AI: the accuracy index
- Over_Sel: the total selection percentage of the overtoxic dose combination(s)
- Over_All: the total percentage of patients treated at the overtoxic dose combination(s)
- summary.MTD.pctg: selection percentage at each dose combination
- summary.patients.pctg: the number of patients treated at each dose combination
- summary.dlt.pctg: the number of toxicities observed at each dose combination

## Authors
Chi Kin Lam, Ruitao Lin and Guosheng Yin 

## Reference
Lam, C.K., Lin R. and Yin, G. (2016) “Nonparametric overdose control for dose finding in drug-combination trials”.
