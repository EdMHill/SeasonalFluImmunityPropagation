# SeasonalFluImmunityPropagation
This repository houses code developed for the analysis presented in the scientific publication "Seasonal influenza: Modelling approaches to capture immunity propagation" by Edward M. Hill, Stavros Petrou, Simon de Lusignan, Ivelina Yonova, and Matt J. Keeling.  

Publication details:
Hill et al. (2019) Seasonal influenza: Modelling approaches to capture immunity propagation. PLOS Computational Biology. 15(10): e1007096. doi: 10.1371/journal.pcbi.1007096.  

Please find below an explainer of the directory structure within this repository.  

# Data

## MortalityData
Population mortality data (from ONS)

## ILIData
Using RCGP data, compute cumulative ILI cases (sum of weekly ILI rate per 100,000) for 2009/10 - 2017/2018 seasons.  
Non-age & age-structured model variants.  
Values are scaled using weekly influenza positivity data.

## RCGPSamplePositivity
Construct arrays storing weekly influenza positivity values.     
Monitored through the RCGP sentinel swabbing scheme in England.

## VaccEfficacy
Population-level vaccine efficacy estimates

## VaccUptake
Daily vacc. uptake proportions for at-risk, low risk, entire population groups  
#### NonAgeStrucModel_DailyVaccUptakeBySeasonCalYr_EMH.xlsx
	Seasonal influenza vaccine uptake by season, for 2008/09 onwards    
#### NonAgeStrucModel_DailyVaccUptakeCalYr_PandemicFluVacc_EMH.xlsx
	Pandemic influenza vaccine uptake for the 2008/09 influenza season  

## WHOFluNet
Data on circulating influenza virus composition from World Health Oragnisation FluNet database

# Results
## EvaluatingModelFit
Scripts to produce bar plots, violin plots and scatter plots.  
ModelSimnData directory contains a file per fit performed, each containing outputs from 1,000 model simulation replicates.

## ForwardProjection
Data files for the two simulated vaccine efficacy scenarios and Fig. 7 plot code.

## MatlabFileExchangeScripts
Community-sourced files that are used to generate the figures.

## ParamInference
Contains scripts to produce end-of-generation threshold value plots and posterior histograms.  
Outputs from the ABC inference scheme reside in the following directories, comprising retained parameter sets and end-of-generation threshold values.  
#### FiveSeasonFitData  
	Output files generated fitting to the empirical data covering 2012/13-2016/17  
#### FourSeasonFitData  
	Output files generated fitting to the empirical data covering 2012/13-2015/16  
#### SixSeasonFitData  
	Output files generated fitting to the empirical data covering 2012/13-2017/18  
#### SynthFitData  
	Output files generated fitting to the synthetic data  

# src

## ParamInference
Perform parameter inference on data using Approximate Bayesian Computation with an Adpative Population Monte Carlo scheme.

#### ABC_APMC
	Directory housing files containing the ABC adaptive population Monte Carlo algorithm

## PosteriorPredictiveModelSimns
Run simulations of the model using parameter sets representing samples from the posterior distribution.

## ModelFns
Functions to numerically solve ODEs and update exposure history array
