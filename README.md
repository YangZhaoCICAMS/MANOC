# Nonparametric Overdose Control for Dose Finding in Drug-Combination Trials
## Inputs

## Functions
- NextDoseComb.R: Containing a function `get.next.manoc(pos.model,j_curr,k_curr,alpha,eta)` for determining the next dose combination given the current dose combination, the posterior model probabilities, alpha and eta. 

- PosteriorProbability.R: Containing a function `posteriorH(y,n,target,p.sample.mat)` for calculating the posterior model probability for each dose combination.
- Simulation.R: Containing a function `simulation(simid,Tox_Prob_Mat,p.sample.mat,samplesize,cohortsize,target,alpha,delta,eta)' for conducting simulation studies. 

- Summarize.R: Containing a function `summarize` for summarizing the outputs produced by the function `simulation()`

- generate_p.sample.mat

