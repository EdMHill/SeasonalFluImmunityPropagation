#Purpose:
#Script to run Adaptive Population Monte Carlo Approximate Bayesian Computation
#Run serially & take arguments from command line as function inputs

#Fit multi-strain, non-age structured influenza transmission model
#(with immunity propagation) to empirical data

#Code Author: Ed Hill
#-------------------------------------------------------------------------------

#--------------------------------------------------------------------------
# ADD FILES TO SEARCH PATH FOR ODES/MODEL RUN FUNCTION
#--------------------------------------------------------------------------
include("../ModelFns/RunSeasonalFluModelODEs.jl")
include("../ModelFns/ExpHistUpdate.jl")
include("ABC_APMC/APMC_ParLoop.jl")
include("SeasonalFluModel_SerialInferenceFns.jl")

#-------------------------------------------------------------------------------
# LOAD REQUIRED PACKAGES
#-------------------------------------------------------------------------------
using Distributions
using StatsBase
using Combinatorics
using Iterators
using DifferentialEquations
using XLSX

function SetUpInferenceAPMC(ARGS)

#Take command line arguments, ARGS, assign to variable names
#ARGS input listing
# ARG 1: RunID
# ARG 2: Observed data filename
# ARG 3: SeasonsToSimulate
# ARG 4: N_alpha
# ARG 5: alpha
# ARG 6: MinAcceptRate
# ARG 7: MaxGen
# ARG 8: PerturbVarFn
# ARG 9: PriorFn
# ARG 10: SummStatFn
# ARG 11: SampleFromPriorFn
# ARG 12: ModelSimnFn
# ARG 13: FirstGenFromFileFlag

#To convert strings to numbers, use parse

#--------------------------------------------------------------------------
# Set RunID idx
#--------------------------------------------------------------------------
RunID = parse(Int64, ARGS[1])

#--------------------------------------------------------------------------
# Specify empirical data being fitted to
#--------------------------------------------------------------------------
ObvsDataFileName = ARGS[2]
ObvsDataTemp = readdlm(ObvsDataFileName,',')

#Specify number of influenza seasonsto be simulated
#From 2009/2010 onwards: - > 7 will be up to and inc. 2015/16
 #                        - > 8 up to and inc. 2016/17 etc
SeasonsToSimulate = parse(Int64, ARGS[3])

#Get maximum number of seasons that could be used in inference procedure
#In case all seasons are not used, get number of seasons omitted
MaxSeasonNumToSimulate = size(ObvsDataTemp,1)
SeasonsUnused =  MaxSeasonNumToSimulate - SeasonsToSimulate

#Error check
IntCheck = rem(SeasonsToSimulate,1)==0
if IntCheck == 0 || SeasonsToSimulate <= 0 || SeasonsToSimulate > MaxSeasonNumToSimulate
    error("Incompatible SeasonsToFit value")
end

#Pick out row 4 onwards of ObvsDataTemp.
#Row 4 corresponds to 2012/13 influenza season, row 5 to 2013/14 influenza season etc
ObvsData = ObvsDataTemp[4:end-SeasonsUnused,:]

#--------------------------------------------------------------------------
# SET UP APMC SCHEME RELATED PARAMETERS
#--------------------------------------------------------------------------

#Specify target number of samples and fraction kept at each iteration
#(alpha)
N_alpha = parse(Int64, ARGS[4])
alpha = parse(Float64, ARGS[5])

#Set minimal acceptance rate
MinAcceptRate = parse(Float64, ARGS[6])
MaxGen = parse(Int64, ARGS[7])

#Specify APMC related functions
s_1 = Symbol(ARGS[8]) #Convert string to Symbol
s_2 = Symbol(ARGS[9]) #Convert string to Symbol
s_3 = Symbol(ARGS[10]) #Convert string to Symbol
s_4 = Symbol(ARGS[11]) #Convert string to Symbol
s_5 = Symbol(ARGS[12]) #Convert string to Symbol

#Make Symbols callable functions
PerturbVarFn = getfield(Main, s_1) #Fn specifying perturbation distribution for newly proposed samples
PriorFn = getfield(Main, s_2) #Fn to check whether perturbed samples are within prior bounds
SummStatFn = getfield(Main, s_3) #Error measure, compare observed to simulated data
SampleFromPriorFn = getfield(Main, s_4) #Sampling from prior distribution
ModelSimnFn = getfield(Main, s_5) #Model simulation

#Indicator variable
# 0 - Sample first particle sets (N in total) from prior distribution
# 1 - 1 - Read particle sets from file, the retained particles at end
#       of a previous inference run (N_alpha in total)
FirstGenFromFileFlag = parse(Int64, ARGS[13])

#Run function
RunInferenceAPMC(RunID,ObvsData,SeasonsToSimulate,
    N_alpha,alpha,MinAcceptRate,MaxGen,PerturbVarFn,PriorFn,SummStatFn,
    SampleFromPriorFn,ModelSimnFn,FirstGenFromFileFlag)
end

#--------------------------------------------------------------------------
# PASS TO FUNCTION
#--------------------------------------------------------------------------
SetUpInferenceAPMC(ARGS)
