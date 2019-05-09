#Purpose:
#Multiple simulations with parameters read from file.
#Run in serial

#Model specifics:
#SEIR influenza deterministic transmission dynamic model ODEs
#Mult-strain, non-age structured model with immunity propagation

#Author: Ed Hill
#--------------------------------------------------------------------------

function RunSeasonalFluModelMultiSimns(ParticleSets,SeedRunID,TotalRunNum,SeasonsToSimulate,
							 		   SimnRunType,VaccEffTestFlag,MultiSeasonImmFlag,OutputFName)
#Inputs:
#   ParticleSets - Parameter sets to be simualted
#   SeedRunID - (Integer) Used to give a unique save filename
#   TotalRunNum - (Integer) Number of simulations, with unique parameter sets, to be
#                   performed
#	SeasonsToSimulate - (Integer) Total number of influenza seasons to be simulated (start season of 2009/10)
#	SimnRunType - (Integer)  Influenza vaccine uptake & efficacy, plus use of  ODEBurnIn, ODEH1N1OnlyTime, ODEAlleStrainTime
#							->	1 - exploratory; 2 - historical; 3 - forward simulations
#	VaccEffTestFlag - (Flag variable) Declare if vaccine efficacy varies per simulation run, or fixed for all
#   					-> 0: Fixed for all
#						-> 1: In forward simulated seasons, randomly sampled each run;
#						-> 2: In forward simulated seasons, take max vacc efficacies from 2010-2018 period
#	MultiSeasonImmFlag -  (Binary indicator) Declare if immunity may be propoagetd for more than a single season
#   					-> 0: No, 1: Yes
#	OutputFName - (string) Save location


#Select exposure history susceptibility modification form for vaccine-related states
#0 - Fixed values every season
#1 - Relate to previous season vaccine efficacy
ExpHistVaccType = 1

if ExpHistVaccType != 0 && ExpHistVaccType !=1
   error("ExpHistVaccType must take value 0 or 1")
end

#--------------------------------------------------------------------------
### DISEASE DYNAMICS SIMN/FLAG PARAMETERS
#--------------------------------------------------------------------------
SimnStartDate=9 #Month of year to start simulation on

#Store population-level FOI flag option
#0 - inactive, 1 - active
StoreFlag_PopnFOI = 0

#Run time for ODE model. Take values post burn in
if SimnRunType == 1
    ODEBurnIn = 0*365 #No entity information recorded
    ODEH1N1OnlyTime = 1*365 #Time spect recording entity values, but only H1N1 infection allowed
    ODEAllStrainTime = 8*365 #Infection by all strains now possible
elseif SimnRunType == 2
    ODEBurnIn = 0*365 #No entity information recorded
    ODEH1N1OnlyTime = 1*365 #Time spect recording entity values, but only H1N1 infection allowed
	ODEAllStrainTime = (SeasonsToSimulate-1)*365 #Infection by all strains now possible
elseif SimnRunType == 3
    ODEBurnIn = 0*365 #No entity information recorded
    ODEH1N1OnlyTime = 1*365 #Time spect recording entity values, but only H1N1 infection allowed
    ODEAllStrainTime = (SeasonsToSimulate-1)*365 #Infection by all strains now possible
else
    error("Incorrect RunType entered")
end

ODESampleTime = ODEH1N1OnlyTime + ODEAllStrainTime #Timeframe over which data outputs are stored
MaxTime = ODEBurnIn + ODESampleTime #Total simlation run time

#Time over which only H1N1 type infection allowed
TotalTime_H1N1Only = ODEBurnIn + ODEH1N1OnlyTime

#Get number of seasons
NumOfSeasons = convert(Int64,ODESampleTime/365)

#Specify timestep to use in ode45 scheme
timestep = 1

#Concatenate simulation variables
SimnParam=[SimnStartDate,ODEBurnIn,ODESampleTime,MaxTime,timestep,TotalTime_H1N1Only,StoreFlag_PopnFOI]

#------------------------------------------------------------------------------
###  TRANSMISSION RELATED PARAMETERS
#------------------------------------------------------------------------------
NumOfStrains = 4 #Specify number of strains in the system (A/H1N1, A/H3N2, two B lineages)
beta = [0.4390,0.5026,0.4263,0.4255] #Placeholder beta values. Will be replaced by laoded parameter set values
gamma = 1/3.8*ones(NumOfStrains) #recovery rate
sigma = [1/1.4,1/1.4,1/0.6,1/0.6]  #latent rate
InfectionParam = [beta sigma gamma]

#Set proportion of population that are susceptible at end of season within
#a natural infection exposure history class to remain in that class
#(rather than move to the Naive exposure history class)
MultiSeasonImmPropn = 0

#------------------------------------------------------------------------------
###  DEMOGRAHPY PARAMETERS
#------------------------------------------------------------------------------
#Number of exposure history classes
#One for naive, one for vacc with no natural infection, one per strain
#(unvacc), one per strain with vacc
ExpHistNum = (NumOfStrains*2) + 2

#Set per capita birth&death rates
life_exp = 81.0
b=1/(life_exp*365)
d=1/(life_exp*365)
BirthParam=zeros(1,ExpHistNum)
BirthParam[1]=b
DeathParam=d

#------------------------------------------------------------------------------
###  VACCINATION UPTAKE
#--------------------------------------------------------------------------
if SimnRunType == 1
    #SYNTHETIC DATA

    #Set vaccine uptake rate
    #Entry (1,kk)- day kk
    VaccUptakeBySeason = zeros(365) #Proportion given the vaccine during each day
	#TargetCovPerSeason = 0.15
	#NumVaccDays = 92
    #VaccUptakeBySeason[274:365] = TargetCovPerSeason/NumVaccDays
elseif SimnRunType == 2
    #BASED ON HISTORICAL DATA

	## Import the data - Pandemic flu vacc (2009/2010 season)
	PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeCalYr_PandemicFluVacc_EMH.xlsx","2009_2010","C5:NC5")

	## Import the data - Sesonal flu vacc
	HistoricalSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeBySeasonCalYr_EMH.xlsx","All popn.","C5:NC12")

    #Combine uptake vectors into a single array
    VaccUptakeBySeason_DataFrame = [PandemicFluVaccUptake; HistoricalSeasonalFluVaccUptake]

    # Collate into Array
    VaccUptakeBySeason = Array{Float64, 2}(VaccUptakeBySeason_DataFrame)


elseif SimnRunType == 3
	#RUN ALTERNATIVE VACC. SCHEMES

	## Import the data - Pandemic flu vacc (2009/2010 season)
	PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeCalYr_PandemicFluVacc_EMH.xlsx","2009_2010","C5:NC5")

	## Import the data - Sesonal flu vacc
	HistoricalSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeBySeasonCalYr_EMH.xlsx","All popn.","C5:NC12")

    #Combine uptake vectors into a single array
    VaccUptakeBySeasonDataFrame = [PandemicFluVaccUptake; HistoricalSeasonalFluVaccUptake]

    # Collate into Array
    VaccUptakeBySeasonUpTo2018 = Array{Float64, 2}(VaccUptakeBySeasonDataFrame)

	#Extend to cover up to and including 2049/2050 influenza season
	VaccUptakeBySeason = zeros(SeasonsToSimulate,365)
	VaccUptakeBySeason[1:9,:] = VaccUptakeBySeasonUpTo2018
	VaccUptakeBySeason[10:end,:] = repmat(VaccUptakeBySeasonUpTo2018[end,:]',SeasonsToSimulate-9,1)

else
    error("Incorrect RunType entered")
end

#--------------------------------------------------------------------------
### VACCINATION - LEAKY TRANSMISSION SETTINGS
#--------------------------------------------------------------------------
#Set flag variable for "leaky" transmission being unactive or active
#0 - Infectiousness of infected vacc. group unmodified.
#1 - Infected vacc. group has reduced infectiousness
LeakyTransFlag = 0

#Validity check on LeakyTransFlag value
if LeakyTransFlag !=0 && LeakyTransFlag !=1
    error("Invalid value of $LeakyTransFlag for LeakyTransFlag, should be 0 or 1")
end

#--------------------------------------------------------------------------
### VACCINATION - EFFICACY
#--------------------------------------------------------------------------
#Set leaky vaccine efficacy parameters
#Row i - Strain i (A/H1N1, A/H3N2, two B lineages)
#Cell 1 - Susceptibility reduction (as propn)
#Cell 2 - infectiousness reduction (as propn)
if SimnRunType == 1
    #SYNTHETIC DATA
    if LeakyTransFlag == 0
        LeakyVaccVarBySeason = zeros(NumOfStrains)
    elseif LeakyTransFlag == 1
        alpha = zeros(NumOfStrains)
        delta = zeros(NumOfStrains)
        LeakyVaccVarBySeason = [alpha delta]
    end
elseif SimnRunType == 2
    #BASED ON HISTORICAL DATA
	VaccEfficacyDataFrame = XLSX.readdata("../../Data/VaccEfficacy/VaccEfficacy_AllPopn.xlsx","MidPoint","C3:F11")

    # Collate into Array
    VaccEfficacy = Array{Float64, 2}(VaccEfficacyDataFrame)

    #Assign to LeakyVaccVar variable
    if LeakyTransFlag == 0
        LeakyVaccVarBySeason = VaccEfficacy'
    elseif LeakyTransFlag == 1
        alpha = VaccEfficacy'
        delta = VaccEfficacy'

        LeakyVaccVarBySeason = [alpha,delta]
    end

elseif SimnRunType == 3
	#RUN ALTERNATIVE VACC. SCHEMES
	VaccEfficacyDataFrame = XLSX.readdata("../../Data/VaccEfficacy/VaccEfficacy_AllPopn.xlsx","MidPoint","C3:F11")

    # Collate into Array
    VaccEfficacyUpTo2018 = Array{Float64, 2}(VaccEfficacyDataFrame)

	#Extend to cover up to and including 2029/2030 influenza season
	VaccEfficacy = zeros(SeasonsToSimulate,NumOfStrains)
	VaccEfficacy[1:9,:] = VaccEfficacyUpTo2018

    #Assign to LeakyVaccVar variable
    if LeakyTransFlag == 0
        LeakyVaccVarBySeason = VaccEfficacy'
    elseif LeakyTransFlag == 1
        alpha = VaccEfficacy'
        delta = VaccEfficacy'

        LeakyVaccVarBySeason = [alpha,delta]
    end
else
    error("Incorrect RunType entered")
end

#--------------------------------------------------------------------------
### INITIAL INF. PROPORTION
#--------------------------------------------------------------------------
#Column i for strain i
#Row 1 for when only H1N1 infection allowed
#Row 2 when infection by any strain allowed
InfPropn_StartOfSeason = [1e-5 0 0 0;
                            2.5e-6 2.5e-6 2.5e-6 2.5e-6]

#--------------------------------------------------------------------------
### GIVE DETAILS ON LOADING INITIAL SUSCEPTIBLE CLASS CONDITIONS FROM FILE IF NEEDED
#--------------------------------------------------------------------------
ICFromFile = [[0],""]

#--------------------------------------------------------------------------
### ASSIGN TO VARIABLE IF THESE SET OF SIMULATIONS ARE FORWARD SIMULATIONS
#--------------------------------------------------------------------------
if SimnRunType == 3
	ForwardSimnFlag = 1
else #Other values of SimnRunType used. Not forward projecting based on historical data.
	ForwardSimnFlag = 0
end

#--------------------------------------------------------------------------
### NUMBER OF INFLUENZA SEASONS WORTH OF DATA BEING CONSIDERED
#--------------------------------------------------------------------------
if ForwardSimnFlag == 0 #check fit to historical data.
	#Discount first three seasons, 2009/10 to 2011/12
	#These seasons not fit to during inference procedure
	RetainedSeasonNum = SeasonsToSimulate - 3
else
	#Running forward simulation. All seasons will be saved to output
	RetainedSeasonNum = SeasonsToSimulate
end

#--------------------------------------------------------------------------
# AGGREGATE FIXED PARAMETERS
#--------------------------------------------------------------------------
FixedModelParams = [SimnRunType,ExpHistVaccType,StoreFlag_PopnFOI,
                    SimnParam,NumOfStrains,ExpHistNum,InfectionParam,MultiSeasonImmPropn,BirthParam,
                    DeathParam,VaccUptakeBySeason,LeakyTransFlag,LeakyVaccVarBySeason,
                    InfPropn_StartOfSeason,ICFromFile,RetainedSeasonNum]

#--------------------------------------------------------------------------
### INITIALISE CELL TO STORE OUTPUTS OF INTEREST
#--------------------------------------------------------------------------
OutputSimnCell = Array{Array{Float64},2}(TotalRunNum,2)

#--------------------------------------------------------------------------
### PERFORM TOTALRUNNUM AMOUNT OF MODEL RUNS
#--------------------------------------------------------------------------
for ii = 1:TotalRunNum

    #----------------------------------------------------------------------
    ### Read in parameter values from file
    #----------------------------------------------------------------------
	if ForwardSimnFlag == 0
		ParticleSetVal = ParticleSets[ii,:]
	else
		ParticleSetVal = ParticleSets[1,:]
		if VaccEffTestFlag == 1

			#Sample efficacies (per strain) uniformly from 2010/2011-2017/2018 empirical data
			VaccEfficacy2010To2018 = VaccEfficacyUpTo2018[2:end,:]

			#----------------------------------------------------------------------
		    ### Sample vaccine efficacyfor 2018/19 to 2029/30 influenza seasons on each model simulation
		    #----------------------------------------------------------------------
			for StrainType = 1:4
				VaccEfficacy[10:end,StrainType] = sample(VaccEfficacy2010To2018[:,StrainType],SeasonsToSimulate-9; replace=true, ordered=false)
			end


		else #VaccEffTestFlag =2
			#Get maximum efficacy value per strain, for seasons 2010/2011-2017/2018
			MaxEffVals = findmax(VaccEfficacyUpTo2018[2:end,:],1)[1]
			VaccEfficacy[10:end,:] = repmat(MaxEffVals,SeasonsToSimulate-9,1)
		end

		#Pass updated vaccine efficacy array into FixedModelParams
		FixedModelParams[13] = VaccEfficacy'
	end

    #----------------------------------------------------------------------
    ### CALL FUNCTION TO RUN MODEL. Outputs two-tuple. First entry
    #----------------------------------------------------------------------

	#Two-tuple output assinged to row ii of OutputSimnCell
	#If ForwardSimnFlag == 0,
	#First entry: Influenza positive GP consultations per 100,000 population by strain
	#Second entry: Entry per timestep giving proportion of population infected

	#If ForwardSimnFlag == 1,
	#First entry: Seasonal influenza incidence rate per 100,000 population (by strain and season)
	#Second entry: Entry per timestep giving proportion of population infected

	if MultiSeasonImmFlag == 0
		OutputSimnCell[ii,:] = RunModelSimn(FixedModelParams,ParticleSetVal)
	elseif MultiSeasonImmFlag == 1
    	OutputSimnCell[ii,:] = RunModelSimnWithMultiSeasonImm(FixedModelParams,ParticleSetVal)
	end
end

#----------------------------------------------------------------------
### SAVE TO OUTPUT VARIABLE
#----------------------------------------------------------------------
OutputFile = matopen(OutputFName,"w")
write(OutputFile, "SimnData", OutputSimnCell)
close(OutputFile)

end

function ModelSimnCallFn(ARGS)

#Take command line arguments, ARGS, assign to variable names
#To convert strings to numbers, use parse

#--------------------------------------------------------------------------
# Read in parameter values from file
#--------------------------------------------------------------------------
ParamInputFile = ARGS[1]
ParticleSets = readdlm(ParamInputFile,'\t')

#----------------------------------------------------------------------
### Specify output file save location
#----------------------------------------------------------------------
SeedRunID = parse(Int64, ARGS[2])
OutputFName = "../../Results/ModelSimnOutputs/ModelSimnRun_#$(SeedRunID).mat"

#--------------------------------------------------------------------------
#Read in vairables for number of simulations and length of each run (in number of influenza seasons)
#--------------------------------------------------------------------------
TotalRunNum = parse(Int64, ARGS[3])
SeasonsToSimulate = parse(Int64, ARGS[4])

#--------------------------------------------------------------------------
#Read in variables for simulation type
#--------------------------------------------------------------------------
SimnRunType = parse(Int64, ARGS[5])
VaccEffTestFlag = parse(Int64, ARGS[6])
MultiSeasonImmFlag = parse(Int64, ARGS[7])

if SimnRunType != 3
	VaccEffTestFlag == 0 #Automatically set VaccEffTestFlag to be 0 irrespective of command line input
end

#--------------------------------------------------------------------------
#Variable value error checks
#--------------------------------------------------------------------------
if SimnRunType!=1 && SimnRunType!=2 && SimnRunType!=3
	error("Value of SimnRunType is $SimnRunType. SimnRunType must take value 1,2 or 3.")
end

if VaccEffTestFlag!=0 && VaccEffTestFlag!=1 && VaccEffTestFlag!=2
	error("Value of VaccEffTestFlag is $VaccEffTestFlag. VaccEffTestFlag must take value 0, 1, or 2.")
end

if MultiSeasonImmFlag!=0 && MultiSeasonImmFlag!=1
	error("Value of MultiSeasonImmFlag is $MultiSeasonImmFlag. MultiSeasonImmFlag must take value 0, 1, or 2.")
end

#--------------------------------------------------------------------------
#Call function to run multiple simulations
#--------------------------------------------------------------------------
RunSeasonalFluModelMultiSimns(ParticleSets,SeedRunID,TotalRunNum,SeasonsToSimulate,
							 SimnRunType,VaccEffTestFlag,MultiSeasonImmFlag,
							 OutputFName)
end


#--------------------------------------------------------------------------
### ADD FILES TO SEARCH PATH FOR ODES/MODEL RUN FUNCTION
include("../ModelFns/RunSeasonalFluModelODEs.jl")
include("../ModelFns/ExpHistUpdate.jl")
include("../ParamInference/SeasonalFluModel_SerialInferenceFns.jl")

# Load libraries
using StatsBase
using Combinatorics
using Iterators
using DifferentialEquations
using DataFrames
using MAT
using XLSX

#--------------------------------------------------------------------------
# PASS TO FUNCTION
#--------------------------------------------------------------------------
ModelSimnCallFn(ARGS)
