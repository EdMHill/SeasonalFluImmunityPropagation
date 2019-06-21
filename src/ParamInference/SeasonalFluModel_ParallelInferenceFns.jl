#Purpose:
#Contains functions to perform Adaptive Population Monte Carlo Approximate Bayesian Computation
#Run in parallel

#Fit multi-strain, non-age structured influenza transmission model
#(with immunity propagation) to empirical data

#Code Author: Ed Hill
#-------------------------------------------------------------------------------

### GROUPING OF FUNCTIONS TO BE PASSED INTO THE APMC SCHEME ###
# (i) FUNCTION TO PERTURB SAMPLES WHEN GENERATING PROPOSED PARAMETER SETS
# (ii) FUNCTIONS TO DRAW SAMPLES FROM PRIOR
# (iii) SPECIFY PRIOR DENSITY
# (iv) SUMMARY STATISTIC CALCULATION FUNCTIONS
# (v) FUNCTION TO RUN MODEL SIMULATION & PRODUCE DESIRED OUTPUTS TO FEED INTO SUMMARY STATISTIC FUNCTION
# (vi) FUNCTIONS USED IN MODEL SIMULATION

#-------------------------------------------------------------------------------
# (i) FUNCTION TO PERTURB SAMPLES WHEN GENERATING PROPOSED PARAMETER SETS
#-------------------------------------------------------------------------------

@everywhere function  OLCMPerturbVarFn(RetainedParams,RetainedWeights,SurvivorParticleIdx,GenT)
# Optimal local covariance matrix (OLCM)
# Uses a multivariate normal distribution with a covariance matrix based on a subset of the particles from the
# previous iteration, whose distances are smaller than the threshold of the current iteration

#Inputs:
#   RetainedParams,RetainedWeights - Current set of samples with associated weights
# 	SurvivorParticleIdx - Boolean true/false states defining the subset of the particles from the previous iteration, whose distances are smaller than the threshold of the current iteratio
#	GenT - Generation number, used to determine if should sample using a global covariance
#Outputs:
#   C - Variance-covariance matrix. Three dimensional, slice per particle

if GenT == 1 #First generation,
	C_SingleSlice = 2.0*cov(RetainedParams,AnalyticWeights(RetainedWeights[:]),corrected=true)::Array{Float64}

	#Check if C is non-singular. Modify if singular.
	while rank(C_SingleSlice) != ParamNum #Check
		C_SingleSlice = C_SingleSlice + 1e-12I
	end

	#Repeat covariance matrix, once per parameter set
	ParticleNum = size(RetainedParams,1) #Number of particles in generation matches number of rows in parameter array
	C = repeat(C_SingleSlice, outer=(1,1,ParticleNum))
else

	#Get parameter sets from the previous iteration, whose distances are smaller than the threshold of the current iteration
	# i.e. "those that survived", row ID of RetainedParams array dennoted by SurvivorParticleIdx
	SPP_thetas = RetainedParams[vec(SurvivorParticleIdx),:] 	#SPP, subset of previous particles

	#Pick out weights for subset of the particles from the
	# previous iteration, whose distances are smaller than the threshold of the current iteration
	# Normalise those weights
	SPP_Weights = RetainedWeights[vec(SurvivorParticleIdx)]
	SPP_NormWeights = SPP_Weights/sum(SPP_Weights)

	#Compute mean of the survivor particle population
    m = sum(SPP_thetas.*SPP_NormWeights,1)

	#Initialise variance-covariance array
	RetainedParamsDims = size(RetainedParams)
	ParamNum = RetainedParamsDims[2]  #Number of parameters matches number of columns in parameter array
	ParticleNum = RetainedParamsDims[1] #Number of particles in generation matches number of rows in parameter array
	C = zeros(ParamNum, ParamNum,ParticleNum)

	#Loop through each particle and compute covariance
	for kk = 1:1:ParticleNum
		Current_C = zeros(ParamNum, ParamNum)
		for jj = 1:1:ParamNum
			for ii = jj:1:ParamNum
				Current_C[ii, jj] = sum(SPP_NormWeights.*(SPP_thetas[:,ii] - m[ii]).*(SPP_thetas[:, jj] - m[jj])) + (m[ii] - RetainedParams[kk,ii])*(m[jj] - RetainedParams[kk,jj])

				#Covariance matrix is symmetric.
				#Assign transposed off-diagonal value
				if ii != jj
					Current_C[jj, ii] = Current_C[ii, jj]
				end
      		end
		end


		#Check if Current_C is non-singular. Blow it up if singular.
		while rank(Current_C) != ParamNum #Check
			Current_C = Current_C + 1e-12I
		end

		#Assign revised covariance to C (collection of covariance arrays)
		C[:,:,kk] = Current_C
	end

	return C
end
end

#-------------------------------------------------------------------------------
# (ii) FUNCTION TO DRAW SAMPLES FROM PRIOR
#-------------------------------------------------------------------------------

@everywhere function  SampleFirstGenFn_FromFile(N::Int64)
#Inputs:
#   N - Number of samples from LHC parameter set to take forward as first generation
    ParticlesFromPrior = readdlm("FILENAME_RetainedParams.txt",'\t')

    #Assign weights to particle
    Particle_Weights = readdlm("FILENAME_RetainedWeights.txt",'\t')

    #Get SurvivorParticleIdx
    SurvivorParticleIdx = readdlm("FILENAME_SurvivorParticleIdx.txt",'\t')

	  #Get summary statistic measure for each paramter set
    Particle_SummStat = readdlm("FILENAME_RetainedSummStat.txt",'\t')

    return ParticlesFromPrior, Particle_Weights, SurvivorParticleIdx, Particle_SummStat
end

#FITTING TO 2012/2013 - 2015/2016 INFLUENZA SEASONS
#Draw samples from specified prior ranges, with ascertainment prob. per
#season
@everywhere function  SampleFirstGenFn_FourSeasonFit(N)
#Inputs:
#   N - Number of samples from prior to take forward as first generation

	#Distributions to be sampled from for each variable
	d = [Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632), #Transmissibility
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05)] #Strain modifier for susceptibility

    ParamToFitNum = length(d) #Total number of parameters to be fitted

    #Sample and assign to array
    ParticlesFromPrior = zeros(N,ParamToFitNum)
    for ii = 1:ParamToFitNum
		ParticlesFromPrior[:,ii] = rand(d[ii],N)
    end

    #Assign weights to particle
    Particle_Weights = ones(N)

    return ParticlesFromPrior, Particle_Weights
end

#FITTING TO 2012/2013 - 2016/2017 INFLUENZA SEASONS
#Draw samples from specified prior ranges, with ascertainment prob. per
#season
@everywhere function  SampleFirstGenFn_FiveSeasonFit(N)
#Inputs:
#   N - Number of samples from prior to take forward as first generation

	#Distributions to be sampled from for each variable
	d = [Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632), #Transmissibility
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05)] #Strain modifier for susceptibility

    ParamToFitNum = length(d) #Total number of parameters to be fitted

    #Sample and assign to array
    ParticlesFromPrior = zeros(N,ParamToFitNum)
    for ii = 1:ParamToFitNum
		ParticlesFromPrior[:,ii] = rand(d[ii],N)
    end

    #Assign weights to particle
    Particle_Weights = ones(N)

    return ParticlesFromPrior, Particle_Weights
end

#FITTING TO 2012/2013 - 2017/2018 INFLUENZA SEASONS
#Draw samples from specified prior ranges, with ascertainment prob. per
#season
@everywhere function  SampleFirstGenFn_SixSeasonFit(N)
#Inputs:
#   N - Number of samples from prior to take forward as first generation

	#Distributions to be sampled from for each variable
	d = [Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632), #Transmissibility
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05)] #Strain modifier for susceptibility

    ParamToFitNum = length(d) #Total number of parameters to be fitted

    #Sample and assign to array
    ParticlesFromPrior = zeros(N,ParamToFitNum)
    for ii = 1:ParamToFitNum
		ParticlesFromPrior[:,ii] = rand(d[ii],N)
    end

    #Assign weights to particle
    Particle_Weights = ones(N)

    return ParticlesFromPrior, Particle_Weights
end


#FITTING TO 2012/2013 - 2017/2018 INFLUENZA SEASONS
#Draw samples from specified prior ranges, with ascertainment prob. per
#season
#Includes parameter for carry-over immunity resulting from natural infection
@everywhere function  SampleFirstGenFn_SixSeasonFitPlusMultiSeasonImm(N)
#Inputs:
#   N - Number of samples from prior to take forward as first generation

	#Distributions to be sampled from for each variable
	d = [Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632), #Transmissibility
		Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05)] #Strain modifier for susceptibility

    ParamToFitNum = length(d) #Total number of parameters to be fitted

    #Sample and assign to array
    ParticlesFromPrior = zeros(N,ParamToFitNum)
    for ii = 1:ParamToFitNum
		ParticlesFromPrior[:,ii] = rand(d[ii],N)
    end

    #Assign weights to particle
    Particle_Weights = ones(N)

    return ParticlesFromPrior, Particle_Weights
end

#----------------------------------------------------------------------
# (iii) SPECIFY PRIOR DENSITY
#----------------------------------------------------------------------

#Per season ascertainment probability, fit to 2012/2013 to 2015/2016
#seasons inclusive
@everywhere function Prior_FourSeasonFit(x)

#Inputs:
#   x - proposed particle value to be tested

#Outputs:
#   PriorProb - Prior likelihood value of current particle (parameter set)
#   InBoundFlag - Denotes whether proposed particle values are within prior bounds

#Get number of dimensions of x
Particle_nDims = ndims(x)

#--------------------------------------------------------------------------
#COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
#use)

#Distributions each varaible were sampled from
d = [Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632), #Transmissibility
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05)] #Strain modifier for susceptibility

if Particle_nDims == 1
	ParticleNum = 1
	ParamNum = length(x)

	PriorProb = pdf.(d[1], x[1])

    for ii = 2:ParamNum
        PriorProb = PriorProb.*pdf.(d[ii], x[ii])
    end

	#Alter flag value for implausible parameter set to 0
	if PriorProb==0
		InBoundFlag = 0::Int64
	else
		InBoundFlag = 1::Int64
	end
else
	ParticleNum = size(x,1)
	ParamNum = size(x,2)

	PriorProb = pdf.(d[1], x[:,1])
    for ii = 2:ParamNum
		PriorProb = PriorProb.*pdf.(d[ii], x[:,ii]);
    end

	#Alter flag value for implausible parameter set to 0
	InBoundFlag = ones(Int64,ParticleNum)
	InBoundFlag[PriorProb.==0] = 0
end

return PriorProb,InBoundFlag

end

#Per season ascertainment probability, fit to 2012/2013 to 2016/2017
#seasons inclusive
@everywhere function Prior_FiveSeasonFit(x)

#Inputs:
#   x - proposed particle value to be tested

#Outputs:
#   PriorProb - Prior likelihood value of current particle (parameter set)
#   InBoundFlag - Denotes whether proposed particle values are within prior bounds

#Get number of dimensions of x
Particle_nDims = ndims(x)

#--------------------------------------------------------------------------
#COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
#use)
#Distributions each varaible were sampled from
d = [Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632), #Transmissibility
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05)] #Strain modifier for susceptibility


if Particle_nDims == 1
	ParticleNum = 1
	ParamNum = length(x)

	PriorProb = pdf.(d[1], x[1])

    for ii = 2:ParamNum
        PriorProb = PriorProb.*pdf.(d[ii], x[ii])
    end

	#Alter flag value for implausible parameter set to 0
	if PriorProb==0
		InBoundFlag = 0::Int64
	else
		InBoundFlag = 1::Int64
	end
else
	ParticleNum = size(x,1)
	ParamNum = size(x,2)

	PriorProb = pdf.(d[1], x[:,1])
    for ii = 2:ParamNum
		PriorProb = PriorProb.*pdf.(d[ii], x[:,ii]);
    end

	#Alter flag value for implausible parameter set to 0
	InBoundFlag = ones(Int64,ParticleNum)
	InBoundFlag[PriorProb.==0] = 0
end

return PriorProb,InBoundFlag

end

#Per season ascertainment probability, fit to 2012/2013 to 2017/2018
#seasons inclusive
@everywhere function Prior_SixSeasonFit(x)

#Inputs:
#   x - proposed particle value to be tested

#Outputs:
#   PriorProb - Prior likelihood value of current particle (parameter set)
#   InBoundFlag - Denotes whether proposed particle values are within prior bounds

#Get number of dimensions of x
Particle_nDims = ndims(x)

#--------------------------------------------------------------------------
#COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
#use)
#Distributions each varaible were sampled from
d = [Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632), #Transmissibility
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05)] #Strain modifier for susceptibility

if Particle_nDims == 1
	ParticleNum = 1
	ParamNum = length(x)

	PriorProb = pdf.(d[1], x[1])

    for ii = 2:ParamNum
        PriorProb = PriorProb.*pdf.(d[ii], x[ii])
    end

	#Alter flag value for implausible parameter set to 0
	if PriorProb==0
		InBoundFlag = 0::Int64
	else
		InBoundFlag = 1::Int64
	end
else
	ParticleNum = size(x,1)
	ParamNum = size(x,2)

	PriorProb = pdf.(d[1], x[:,1])
    for ii = 2:ParamNum
		PriorProb = PriorProb.*pdf.(d[ii], x[:,ii]);
    end

	#Alter flag value for implausible parameter set to 0
	InBoundFlag = ones(Int64,ParticleNum)
	InBoundFlag[PriorProb.==0] = 0
end

return PriorProb,InBoundFlag

end

#Per season ascertainment probability, fit to 2012/2013 to 2017/2018
#seasons inclusive
#Includes parameter for carry-over immunity resulting from natural infection
@everywhere function Prior_SixSeasonFitPlusMultiSeasonImm(x)

#Inputs:
#   x - proposed particle value to be tested

#Outputs:
#   PriorProb - Prior likelihood value of current particle (parameter set)
#   InBoundFlag - Denotes whether proposed particle values are within prior bounds

#Get number of dimensions of x
Particle_nDims = ndims(x)

#--------------------------------------------------------------------------
#COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
#use)
#Distributions each varaible were sampled from
d = [Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632),Uniform(0.2632,3*0.2632), #Transmissibility
		Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05),Uniform(0,0.05)] #Strain modifier for susceptibility

if Particle_nDims == 1
	ParticleNum = 1
	ParamNum = length(x)

	PriorProb = pdf.(d[1], x[1])

    for ii = 2:ParamNum
        PriorProb = PriorProb.*pdf.(d[ii], x[ii])
    end

	#Alter flag value for implausible parameter set to 0
	if PriorProb==0
		InBoundFlag = 0::Int64
	else
		InBoundFlag = 1::Int64
	end
else
	ParticleNum = size(x,1)
	ParamNum = size(x,2)

	PriorProb = pdf.(d[1], x[:,1])
    for ii = 2:ParamNum
		PriorProb = PriorProb.*pdf.(d[ii], x[:,ii]);
    end

	#Alter flag value for implausible parameter set to 0
	InBoundFlag = ones(Int64,ParticleNum)
	InBoundFlag[PriorProb.==0] = 0
end

return PriorProb,InBoundFlag

end

#-------------------------------------------------------------------------------
# (iv) SPECIFY SUMMARY STATISTIC CALCULATION
#-------------------------------------------------------------------------------
function  SummStatFun(ObservedData,SimnData)

    #Disaggregate simulation data
    x = SimnData[1]
    Total_I = SimnData[2]

	#Compute poisson deviance
	SummStatIterNum = length(ObservedData)
	TempSum = 0.0
    for ii = 1:SummStatIterNum
        if ObservedData[ii].!=0
            TempSum = TempSum + (ObservedData[ii]*log(ObservedData[ii]/x[ii])) - (ObservedData[ii] - x[ii])
        end
    end
    SummStatVal_Overall = 2*(TempSum + sum(x[ObservedData.==0]))

    #Perform temporal check
    NumOfSeasons = convert(Int64,size(Total_I,1)/366) #%Number of seasons obtained by dividing number of daily records by days in yr (+1 to account for day 0 recording1)
    Total_I_Array = zeros(NumOfSeasons,366)
    for jj = 1:NumOfSeasons
       StartIdx = ((jj-1)*366) + 1
       EndIdx = jj*366

       Total_I_Array[jj,:] = Total_I[StartIdx:EndIdx]
    end

	#Find day of seasonal year in which infection is at peak
	MaxInfValBySeason,inds = findmax(Total_I_Array[4:end,:],2)
	MaxInfIdxBySeason = map(x->ind2sub(Total_I_Array[4:end,:], x)[2], inds)

    #If peak outside Sep-Feb, set TemporalFlag to 0 and put amended
    #SummStatVal as infinity
    if sum(MaxInfIdxBySeason.>182) == 0 #181 days September-February. Plus account for initial value
        AmendedSummStatVal = SummStatVal_Overall
    else
        AmendedSummStatVal = Inf
    end

    return AmendedSummStatVal

end

#-------------------------------------------------------------------------------
# (v) FUNCTION TO RUN MODEL SIMULATION & PRODUCE DESIRED OUTPUTS TO FEED INTO
# SUMMARY STATISTIC FUNCTION
#-------------------------------------------------------------------------------
@everywhere function  RunModelSimn(FixedModelParams,x)
#Inputs:
#   FixedModelParams - Parameters with consistent values across all runs
#   x - Values of parameters that are being inferred

#Outputs:
#   SimnData - Model output to be fed into summary statistic function

    #Disaggregate FixedModelParams inputs
    SimnRunType = FixedModelParams[1]
    ExpHistVaccType = FixedModelParams[2]
    StoreFlag_PopnFOI = FixedModelParams[3]
    SimnParam = FixedModelParams[4]
    NumOfStrains = FixedModelParams[5]
    ExpHistNum = FixedModelParams[6]
    InfectionParam = FixedModelParams[7]
    MultiSeasonImmPropn = FixedModelParams[8]
    BirthParam = FixedModelParams[9]
    DeathParam = FixedModelParams[10]
    VaccUptakeBySeason = FixedModelParams[11]
    LeakyTransFlag = FixedModelParams[12]
    LeakyVaccVarBySeason = FixedModelParams[13]
    InfPropn_StartOfSeason = FixedModelParams[14]
    ICFromFile = FixedModelParams[15]
    RetainedSeasonNum = FixedModelParams[16]

    #Update beta!
    InfectionParam[:,1] = x[1:4]

    #Update ascertainment prob

    #Update exposure history
    #Build exposure history array. Assign to variable
    ExpHistArrayParams = x[5:7]
    ExpHistArray = BuildExpHistArray(NumOfStrains,ExpHistArrayParams)
	ExpHistArrayFnInputs = [ExpHistArray,ExpHistArrayParams]

	#Run the model!
	StoreArrays = RunSeasonalFluModel(SimnParam,
   	    InfectionParam,BirthParam,DeathParam,VaccUptakeBySeason,LeakyVaccVarBySeason,LeakyTransFlag,
   	    ExpHistArrayFnInputs,ExpHistNum,ExpHistVaccType,MultiSeasonImmPropn,
           InfPropn_StartOfSeason,ICFromFile,SimnRunType)

	#Disaggregate StoreArrays
	T = StoreArrays[1]::Array{Float64,1}
	C = StoreArrays[2]::Array{Float64,2}
	E_NotV = StoreArrays[5]::Array{Float64,2}
	E_V = StoreArrays[6]::Array{Float64,2}
	I_NotV = StoreArrays[7]::Array{Float64,2}
	I_V = StoreArrays[8]::Array{Float64,2}

	#Visits in week (or month) t, denoted c_{t}, difference in the cumulative proportion
	#of consultations over the previous timestep
	#i.e. c_{t} = p(C(t)-C(t-1)),
	CumulCaseCount = [C[1:366:end,:];C[end,:]'] #Get cumul. case count at specified intervals
	StrainSeasonRateTotalSum = CumulCaseCount[2:end,:]-CumulCaseCount[1:end-1,:]
	StrainSeasonRateTotalSum_FitCheckSeasons = StrainSeasonRateTotalSum[4:end,:] #Get model simulated values for 2012/2013 onward

	#----------------------------------------------------------------------
	### Get ascertainable cases, based on ascertainment prob type
	#----------------------------------------------------------------------
    AscertainProb = x[8:end] #Parameter values from LHC file
    StrainSeasonRatePerStrain = zeros(RetainedSeasonNum,4) #Initialise storage array

	#Update ascertainable amount per strain
	FluToGPconversion = 100000*AscertainProb #Scale to rate per 100,000 popualation
	for jj = 1:length(AscertainProb)
		StrainSeasonRatePerStrain[jj,:] = StrainSeasonRateTotalSum_FitCheckSeasons[jj,:]*FluToGPconversion[jj]
	end

    #Compute infected temporal profile
    Total_I = sum(I_NotV,2) + sum(I_V,2) + sum(E_NotV,2) + sum(E_V,2)

    #Assign outputs to be used in Summary Statisitc function to array
    SimnData = [StrainSeasonRatePerStrain,Total_I]

	return SimnData
end

#Includes fitting parameter for carry-over immunity resulting from natural infection
@everywhere function  RunModelSimnWithMultiSeasonImm(FixedModelParams,x)
#Inputs:
#   FixedModelParams - Parameters with consistent values across all runs
#   x - Values of parameters that are being inferred

#Outputs:
#   SimnData - Model output to be fed into summary statistic function

    #Disaggregate FixedModelParams inputs
    SimnRunType = FixedModelParams[1]
    ExpHistVaccType = FixedModelParams[2]
    StoreFlag_PopnFOI = FixedModelParams[3]
    SimnParam = FixedModelParams[4]
    NumOfStrains = FixedModelParams[5]
    ExpHistNum = FixedModelParams[6]
    InfectionParam = FixedModelParams[7]
    MultiSeasonImmPropn = FixedModelParams[8]
    BirthParam = FixedModelParams[9]
    DeathParam = FixedModelParams[10]
    VaccUptakeBySeason = FixedModelParams[11]
    LeakyTransFlag = FixedModelParams[12]
    LeakyVaccVarBySeason = FixedModelParams[13]
    InfPropn_StartOfSeason = FixedModelParams[14]
    ICFromFile = FixedModelParams[15]
    RetainedSeasonNum = FixedModelParams[16]

    #Update beta!
    InfectionParam[:,1] = x[1:4]

    #Update ascertainment prob

    #Update exposure history
    #Build exposure history array. Assign to variable
    ExpHistArrayParams = x[5:7]
    ExpHistArray = BuildExpHistArray(NumOfStrains,ExpHistArrayParams)
	ExpHistArrayFnInputs = [ExpHistArray,ExpHistArrayParams]

	MultiSeasonImmPropn = x[8]::Float64

	#Run the model!
	StoreArrays = RunSeasonalFluModel(SimnParam,
   	    InfectionParam,BirthParam,DeathParam,VaccUptakeBySeason,LeakyVaccVarBySeason,LeakyTransFlag,
   	    ExpHistArrayFnInputs,ExpHistNum,ExpHistVaccType,MultiSeasonImmPropn,
           InfPropn_StartOfSeason,ICFromFile,SimnRunType)


	#Disaggregate StoreArrays
	T = StoreArrays[1]::Array{Float64,1}
	C = StoreArrays[2]::Array{Float64,2}
	E_NotV = StoreArrays[5]::Array{Float64,2}
	E_V = StoreArrays[6]::Array{Float64,2}
	I_NotV = StoreArrays[7]::Array{Float64,2}
	I_V = StoreArrays[8]::Array{Float64,2}


	#Visits in week (or month) t, denoted c_{t}, difference in the cumulative proportion
	#of consultations over the previous timestep
	#i.e. c_{t} = p(C(t)-C(t-1)),
	CumulCaseCount = [C[1:366:end,:];C[end,:]'] #Get cumul. case count at specified intervals
	StrainSeasonRateTotalSum = CumulCaseCount[2:end,:]-CumulCaseCount[1:end-1,:]
	StrainSeasonRateTotalSum_FitCheckSeasons = StrainSeasonRateTotalSum[4:end,:] #Get model simulated values for 2012/2013 onward


	#----------------------------------------------------------------------
	### Get ascertainable cases, based on ascertainment prob type
	#----------------------------------------------------------------------
    AscertainProb = x[9:end] #Parameter values from LHC file
    StrainSeasonRatePerStrain = zeros(RetainedSeasonNum,4) #Initialise storage array

	#Forward simulation modifications
	#Save counts unmodified by ascertainment probabilities
    if SimnRunType == 3
        StrainSeasonRatePerStrain = StrainSeasonRateTotalSum_FitCheckSeasons
	else
		#Update ascertainable amount per strain
		FluToGPconversion = 100000*AscertainProb #Scale to rate per 100,000 popualation
		for jj = 1:length(AscertainProb)
			StrainSeasonRatePerStrain[jj,:] = StrainSeasonRateTotalSum_FitCheckSeasons[jj,:]*FluToGPconversion[jj]
		end
    end

    #Compute infected temporal profile
    Total_I = sum(I_NotV,2) + sum(I_V,2) + sum(E_NotV,2) + sum(E_V,2)

    #Assign outputs to be used in Summary Statisitc function to array
    SimnData = [StrainSeasonRatePerStrain,Total_I]

	return SimnData
end

#----------------------------------------------------------------------
# (vi) FUNCTIONS USED IN SUMMARY STATISTIC CALCULATION
#----------------------------------------------------------------------

@everywhere function BuildExpHistArray(NumOfStrains,ExpHistArrayParams)
#PURPOSE: Construct exposure history array based on current parameter set
#Inputs:
#   NumOfStrains - As titled
#   ExpHistArrayParams - File containing parameter sets to be run
#Outputs:
#   ExpHistArray - interaction array between exposure history and susceptibility to
#   the current season strain variant

	#Number of exposure history classes
	#One for naive, one for vacc with no natural infection, one per strain
	#(unvacc), one per strain with vacc
	ExpHistNum = (NumOfStrains*2) + 2

	#Column per exposure history. Row per strain.
	ExpHistArray = ones(NumOfStrains,ExpHistNum)
	HalfExpHistNum = convert(Int,ExpHistNum/2)

	#Column 1 - no exposure in previous season
	#Columns 2-5 - natural infection by one of the strains
	#Column 6  - vacc. in previous season, no natural infection
	#Columns 7-10 - Natural infection by ones of the sratins AND vaccinated

	#Right hand columns, vaccinated previous season
	ExpHistArray[:,HalfExpHistNum+1:end] = ExpHistArrayParams[3]

	for ii = 2:HalfExpHistNum
	    #Update entries: natural infection
	    ExpHistArray[ii-1,ii] = ExpHistArrayParams[1]
	    ExpHistArray[ii-1,ii+HalfExpHistNum] = ExpHistArrayParams[1]
	end

	#Update entries: influenza B corss-reactivity
	ExpHistArray[NumOfStrains-1,[HalfExpHistNum end]] = ExpHistArrayParams[2]
	ExpHistArray[NumOfStrains,[HalfExpHistNum-1 end-1]] = ExpHistArrayParams[2]

	return ExpHistArray
end


function RunInferenceAPMC(RunID,ObvsData,SeasonsToSimulate,
                        N_alpha,alpha,MinAcceptRate,MaxGen,PerturbVarScale,
						PriorFn,SummStatFn,SampleFromPriorFn,ModelSimnFn,
						FirstGenFromFileFlag)

#-------------------------------------------------------------------------------
# SET UP APMC SCHEME RELATED PARAMETERS
#-------------------------------------------------------------------------------

#Calculate number of samples before retention phase
ScalingFactor = 1/alpha
N = convert(Int64,round(N_alpha*ScalingFactor))

#-------------------------------------------------------------------------------
# DEFINE AND GROUP MODEL SIMULATION FIXED PARAMETERS
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SPECIFY TYPE OF RUN THROUGH FLAG VARIABLE
#-------------------------------------------------------------------------------

#TYPE OF RUN
# (INFLUENCES VACCINE UPTAKE/EFFICACY, & USE OF ODEBurnIn, ODEH1N1OnlyTime, ODEAlleStrainTime)
#1 - exploratory; 2 - historical; 3 - alternative vacc. scheme
SimnRunType = 2 #Fixed to 2, as fitting to historical data

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
ODEBurnIn = 0*365 #No entity information recorded
ODEH1N1OnlyTime = 1*365 #Time spect recording entity values, but only H1N1 infection allowed
ODEAllStrainTime = (SeasonsToSimulate-1)*365 #Infection by all strains now possible

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
beta = [0.4390,0.5026,0.4263,0.4255] #Placeholder beta values. Will be replaced during inference process.
gamma = 1/3.8*ones(NumOfStrains) #recovery rate
sigma = [1/1.4,1/1.4,1/0.6,1/0.6]  #latent rate
InfectionParam = [beta sigma gamma]

#Set proportion of population that are susceptible at end of season within
#a natural infection exposure history class to retain pre-existing immunity
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

#BASED ON HISTORICAL DATA

## Import the data - Pandemic flu vacc (2009/2010 season)
PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeCalYr_PandemicFluVacc_EMH.xlsx","2009_2010","C5:NC5")

## Import the data - Sesonal flu vacc
HistoricalSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeBySeasonCalYr_EMH.xlsx","All popn.","C5:NC12")

#Combine uptake vectors into a single array
VaccUptakeBySeason_DataFrame = [PandemicFluVaccUptake; HistoricalSeasonalFluVaccUptake]

# Collate into Array
VaccUptakeBySeason = Array{Float64, 2}(VaccUptakeBySeason_DataFrame)

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
### NUMBER OF INFLUENZA SEASONS WORTH OF DATA BEING CONSIDERED
#--------------------------------------------------------------------------
RetainedSeasonNum = size(ObvsData,1)

#--------------------------------------------------------------------------
# AGGREGATE FIXED PARAMETERS
#--------------------------------------------------------------------------
FixedModelParams = [SimnRunType,ExpHistVaccType,StoreFlag_PopnFOI,
                    SimnParam,NumOfStrains,ExpHistNum,InfectionParam,MultiSeasonImmPropn,BirthParam,
                    DeathParam,VaccUptakeBySeason,LeakyTransFlag,LeakyVaccVarBySeason,
                    InfPropn_StartOfSeason,ICFromFile,RetainedSeasonNum]

#--------------------------------------------------------------------------
# SET END-OF-GENERATION OUTPUT TEXT FILE INFO
#--------------------------------------------------------------------------
OutputFileName = "../../Results/ParamInference/APMCOutputFiles/EndOfGenAPMC_"

#--------------------------------------------------------------------------
# RUN APMC SCHEME
#--------------------------------------------------------------------------
RetainedParams,RetainedWeights,RetainedSummStat,SurvivorParticleIdx,GenT = APMC_ParLoop(ObvsData,FirstGenFromFileFlag,
	N,alpha,N_alpha,MinAcceptRate,MaxGen,PerturbVarScale,
    PriorFn,SummStatFn,SampleFromPriorFn,FixedModelParams,ModelSimnFn,RunID,OutputFileName)

#-------------------------------------------------------------------------------
# SAVE TO FILE
#-------------------------------------------------------------------------------

#Specify filename to write sample values
FName_GenT = "../../Results/ParamInference/APMCOutputFiles/APMCsamples_JulRun#$(RunID)_GenT.txt"
FName_RetainedParams = "../../Results/ParamInference/APMCOutputFiles/APMCsamples_JulRun#$(RunID)_RetainedParams.txt"
FName_RetainedSummStat = "../../Results/ParamInference/APMCOutputFiles/APMCsamples_JulRun#$(RunID)_RetainedSummStat.txt"
FName_RetainedWeights = "../../Results/ParamInference/APMCOutputFiles/APMCsamples_JulRun#$(RunID)_RetainedWeights.txt"
FName_SurvivorParticleIdx = "../../Results/ParamInference/APMCOutputFiles/APMCsamples_JulRun#$(RunID)_SurvivorParticleIdx.txt"

#--------------------------------------------------------------------------
### SAVE SUMMARY STATISTICS TO FILE
#--------------------------------------------------------------------------
writedlm(FName_GenT,GenT)
writedlm(FName_RetainedParams,RetainedParams)
writedlm(FName_RetainedSummStat,RetainedSummStat)
writedlm(FName_RetainedWeights,RetainedWeights)
writedlm(FName_SurvivorParticleIdx,SurvivorParticleIdx)

end
