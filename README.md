# Multi-agent Nonparametric Overdose Contraol (MANOC)
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
- NextDoseComb.R: Containing a function `get.next.manoc` for determining the next dose combination given the current dose combination, the posterior model probabilities, alpha and eta. 

- PosteriorProbability.R: Containing a function `posteriorH` for calculating the posterior model probability for each dose combination.

- Simulation.R: Containing a function `simulation(simid,Tox_Prob_Mat,p.sample.mat,samplesize,cohortsize,target,alpha,delta,eta)` for conducting simulation studies. 

- Summarize.R: Containing a function `summarize` for summarizing the outputs produced by the function `simulation()`.

- ToxProb_Generate.R: Containing a function `generate_p.sample.mat(\delta)` for generating samples of the toxicity matrix **p** from its prior distribution. Details can be found in the Appendix of the paper. 

## Examples
### Posterior Probability
```
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

```
> round(PostProb,digits=2)
     [,1] [,2] [,3] [,4] [,5]
[1,] 0.00 0.00 0.00 0.00 0.03
[2,] 0.00 0.00 0.00 0.03 0.45
[3,] 0.00 0.00 0.05 0.17 0.03
[4,] 0.01 0.06 0.14 0.04 0.00
```

## Authors and Reference
Chi Kin Lam, Ruitao Lin and Guosheng Yin 

## Reference
Lam, C.K., Lin R. and Yin, G. (2016) “Nonparametric overdose control for dose finding in drug-combination trials”.
