# MANOC (Multi-agent Nonparametric Overdose Control)
R codes to implement the multi-agent nonparametric overdose control design for dose finding in phase I drug-combination trials.

## Inputs 
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

## Functions
- NextDoseComb.R: Containing a function `get.next.manoc()` for determining the next dose combination given the current dose combination, the posterior model probabilities, alpha and eta. 

- PosteriorProbability.R: Containing a function `posteriorH()` for calculating the posterior model probability for each dose combination.

- Simulation.R: Containing a function `simulation()` for conducting simulation studies. 

- Summarize.R: Containing a function `summarize()` for summarizing the outputs produced by the function `simulation()`.

- ToxProb_Generate.R: Containing a function `generate_p.sample.mat()` for generating samples of the toxicity matrix **p** from its prior distribution. Details can be found in the Appendix of the paper. 

## Examples
### Posterior Probability
```rscript
> setwd("/MANOC_master/")
> source("ToxProb_Generate.R")
> source("PosteriorProbability.R")
> 
> target <- 0.30
> epi <- 0.025
> NN <- 10000
> 
> ## Generate the matrices of p. ## 
> p.sample.mat <- generate_p.sample.mat(ndose.A=4,ndose.B=5, NN=NN, target=target, epi=epi) 
> 
> n<-matrix(c(3,0,0,0,0,3,0,0,0,0,3,0,0,9,3,3,3,39,0,0),4,5)
> y<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,2,0,12,0,0),4,5)
> 
> PostProb<-posteriorH(y=y,n=n,target=0.3,p.sample.mat=p.sample.mat)$lik
```

```rscript
> round(PostProb,digits=2)
     [,1] [,2] [,3] [,4] [,5]
[1,] 0.00 0.00 0.00 0.00 0.03
[2,] 0.00 0.00 0.00 0.03 0.45
[3,] 0.00 0.00 0.05 0.17 0.03
[4,] 0.01 0.06 0.14 0.04 0.00
```
### Next Dose Level
```rscript
> rm(list=ls())
> setwd("/MANOC_master/")
> source("ToxProb_Generate.R")
> source("NextDoseComb.R")
> 
> target <- 0.33
> epi <- 0.025
> delta <- 0.05
> NN <- 50000
> 
> j_curr<-1
> k_curr<-5
> 
> alpha<-0.35
> eta<-0.55
>
> n<-matrix(c(3,0,0,0, 0,3,0,0, 0,0,3,0, 0,9,3,3, 3,3,0,0),4,5)
> y<-matrix(c(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,2,2, 0,2,0,0),4,5)
> ## Generate the matrices of p. ## 
> p.sample.mat <- generate_p.sample.mat(ndose.A=nrow(n),ndose.B=ncol(n), NN=NN, target=target, epi=epi) 
> 
> get.next.manoc(y=y,n=n,target=target,delta=delta,p.sample.mat=p.sample.mat,j_curr=j_curr,k_curr=k_curr,alpha=alpha,eta=eta)
$next.dose
[1] 2 5
```


### Simulation Studies
## Authors and Reference
Chi Kin Lam, Ruitao Lin and Guosheng Yin 

## Reference
Lam, C.K., Lin R. and Yin, G. (2016) “Nonparametric overdose control for dose finding in drug-combination trials”.
