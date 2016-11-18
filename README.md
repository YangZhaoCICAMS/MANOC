# Multi-agent Nonparametric Overdose Contraol (MANOC)
## Inputs 
- samplesize: The maximum number of patients to be enrolled.  
- cohortsize: The number of patients in each cohort. 
- target: The target toxicity rate. 
- epi: A small positive number epsilon in the model specification.  
- alpha: 
- delta: A small increment for the untried dose combinations.
- eta: dose-switching cutoff.
- nsim:
- NN: The number of samples of **p** generated from its prior distribution.
- Tox_Prob_Mat:

## Functions
- NextDoseComb.R: Containing a function `get.next.manoc(pos.model,j_curr,k_curr,alpha,eta)` for determining the next dose combination given the current dose combination, the posterior model probabilities, alpha and eta. 

- PosteriorProbability.R: Containing a function `posteriorH(y,n,target,p.sample.mat)` for calculating the posterior model probability for each dose combination.

- Simulation.R: Containing a function `simulation(simid,Tox_Prob_Mat,p.sample.mat,samplesize,cohortsize,target,alpha,delta,eta)` for conducting simulation studies. 

- Summarize.R: Containing a function `summarize` for summarizing the outputs produced by the function `simulation()`.

- ToxProb_Generate.R: Containing a function `generate_p.sample.mat(\delta)` for generating samples of the toxicity matrix **p** from its prior distribution. Details can be found in the Appendix of the paper. 

## Examples

## Authors and Reference
Chi Kin Lam, Ruitao Lin and Guosheng Yin 

## Reference
Lam, C.K., Lin R. and Yin, G. (2016) “Nonparametric overdose control for dose finding in drug-combination trials”.
