function APMC_ParLoop(ObservedData,FirstGenFromFileFlag::Int64,N::Int64,alpha::Float64,N_alpha::Int64,
                MinAcceptRate::Float64,MaxGen::Int64,PerturbVarFn::Function,
                PriorFn::Function,SummStatFn::Function,SampleFirstGenFn::Function,
                FixedModelParams,ModelSimnFn::Function,RunID,OutputFileName)
#Inputs:
#   ObservedData - (float array) The data!
#   FirstGenFromFileFlag - (binary) Indicator for method of obtaining initial sample (1 from file with particles from previous APMC run, 0 to draw from prior))
#   N - (Int scalar) Number of particles in generation pre retention phase
#   alpha - (Float scalar) Proportion of samples kept in retention phase
#   N_alpha -  (Int scalar) Number of retained samples
#   MinAcceptRate - (Float scalar) Stopping criterion: Threshold for propn of newly genertated samples that
#                   must attain superior summary statistic to threshold summary statistic
#                   value epsilon for simulation to continue
#   PerturbVarScale - (Float) Factor to scale the weighted empirical variance, used as covariance matrix in Gaussian distribution when perturbing samples
#   MaxGen - (Int scalar) Forcibly stop inference procedure if this amount of geneations
#                   is reached
#   PerturbVarFn - (fn handle) Used as covariance matrix in Gaussian distribution when perturbing samples
#   PriorFun - (fn handle) Prior calculation & check if perturbed sample is in prior bounds!
#   SummStatFun - (fn handle) summary statistic calculation
#   SampleFirstGenFn - (fn handle) specifying generation of first generation
#                           of particles & associated weights
#                      -> generated from prior distribution, equal weights
#                      -> import samples/weights from file
#   FixedModelParams - (Tuple) Batch of quantities to be passed into model simulation function
#   ModelSimnFn - (fn handle) Perform model simulation
#   RunID - (Int scalar) Identifier for APMC run
#   OutputFileName - (string) Directory location for end-of-generation data to be
#                     saved to

#Outputs:
# RetainedParams,RetainedWeights,RetainedSummStat - Result: N_alpha
#           particles with assoicated weights and summary statistic values
# SurvivorParticleIdx - Logical (true/false) vector, gives those N_alpha particles
#                       NOT newly generated in final generation
# GenT - Number of generations until stopping criterion satisfied.

#-------------------------------------------------------------------------------
#RUN INITIAL LOOP
#-------------------------------------------------------------------------------

if FirstGenFromFileFlag==1 #Get first generation particles from previous inference run
    #Initialisation
  RetainedParams,Particle_Weights_Temp,SurvivorParticleIdxTemp,SummStat_Temp = SampleFirstGenFn(N_alpha)

  #Convert SurvivorParticleIdx to BitArray
  SurvivorParticleIdx = convert(BitArray,SurvivorParticleIdxTemp)

  #Get number of parameters being inferred
  ParamInferredNum = size(RetainedParams,2)
  FullParamSet = zeros(N,ParamInferredNum) #Initialise full parameter set array


  #Assign weights to particle
  Particle_Weights = SharedArray{Float64}(N)
  Particle_Weights[1:end] = ones(N)
  Particle_Weights[1:N_alpha] = Particle_Weights_Temp

  #Initialise summary statistic storage vector
  SummStat = SharedArray{Float64}(N)
  SummStat[1:end] = zeros(N)
  SummStat[1:N_alpha] = SummStat_Temp

  #Return the first alpha-quartile of SummStat
  epsilon = SummStat[N_alpha]

  #Obtain first-generation of accepted particles
  RetainedWeights = Particle_Weights[1:N_alpha]
  RetainedSummStat = SummStat[1:N_alpha]

  #Set covariance matrix from accepted samples (based on information from file)
  GenNum = 1
  C = PerturbVarFn(RetainedParams,RetainedWeights,SurvivorParticleIdx,GenNum)
else #Get first generation particles from prior
    #Initialisation
  FullParamSet,Particle_Weights_Temp = SampleFirstGenFn(N)

  #Assign weights to particle
  Particle_Weights = SharedArray{Float64}(N)
  Particle_Weights[1:end] = Particle_Weights_Temp

  #Initialise summary statistic storage vector
  SummStat = SharedArray{Float64}(N)

  # First tolerance, sample from prior
  @sync @parallel for ii = 1:N

      #Run model with designated parameter set
      SimnData = ModelSimnFn(FixedModelParams,FullParamSet[ii,:])

      #Generate summary statistic for current prior set
      SoS = SummStatFn(ObservedData,SimnData)

      if isnan(SoS)     #Convert implausible (NaN) outputs to inf
          SummStat[ii] = Inf
      else
          SummStat[ii] = SoS
      end
  end

  #Find particle idxs that are to be retained (below epsilon)
  EpsilonIdx = sortperm(SummStat)
  ParamRetained = EpsilonIdx[1:N_alpha]

  #Obtain first-generation of accepted particles
  RetainedParams = FullParamSet[ParamRetained,:]
  RetainedWeights = Particle_Weights[ParamRetained]
  RetainedSummStat = SummStat[ParamRetained]

  #Return the first alpha-quartile of SummStat
  epsilon = RetainedSummStat[end]
  #epsilon = percentile(SummStat,alpha*100)

  #Set covariance matrix from accepted samples (in first generation set)
  GenNum = 0
  SurvivorParticleIdx = 0 #Placeholder value. As first generation SurvivorParticleIdx not used within PerturbVarFn
  C = PerturbVarFn(RetainedParams,RetainedWeights,SurvivorParticleIdx,GenNum)
end

#Reorder FullParamSet, Weights, SummStat. Place retained values before unretained
FullParamSet[1:N_alpha,:] = RetainedParams
Particle_Weights[1:N_alpha] = RetainedWeights
SummStat[1:N_alpha] = RetainedSummStat

#Set acceptance rate to 1
AcceptRate = 1.

#Set up generation index
GenT = 1

#-------------------------------------------------------------------------------
# Write param values, weights, etc to output files
#-------------------------------------------------------------------------------
ParamSetFileName = "$(OutputFileName)_ParamSets_Run#$(RunID).txt"
WeightsFileName = "$(OutputFileName)_Weights_Run#$(RunID).txt"
SummStatFileName = "$(OutputFileName)_SummStat_Run#$(RunID).txt"
ThresholdValFileName = "$(OutputFileName)_ThresholdVal_Run#$(RunID).txt"
AccRateFileName = "$(OutputFileName)_AccRate_Run#$(RunID).txt"
SurvivorParticleIdxFileName = "$(OutputFileName)_SurvivorParticleIdx_Run#$(RunID).txt"

writedlm(ParamSetFileName,RetainedParams)
writedlm(WeightsFileName,RetainedWeights)
writedlm(SummStatFileName,RetainedSummStat)
writedlm(ThresholdValFileName,epsilon)
writedlm(AccRateFileName,AcceptRate)
writedlm(SurvivorParticleIdxFileName,SurvivorParticleIdx)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# RUN APMC LOOP
#-------------------------------------------------------------------------------
ParamNum = size(RetainedParams,2) #Number of parameters

while (AcceptRate > MinAcceptRate) && (GenT < MaxGen)
    #Update generation index
    GenT = GenT + 1

    #----------------------------------------------------------------------
    # Sample N-N_alpha particles from the previous weighted sample
    #----------------------------------------------------------------------

    #Generate random numbers to be used for drawing from previous parameter set
    #generation
    PickParticleRand = rand(N - N_alpha)

    #Construct vector of particle weights.
    CumSumWeights = cumsum(RetainedWeights)/sum(RetainedWeights)

    #Draw from previous parameter set generation
    ParamSetsDrawn = Array{Float64,2}(N-N_alpha,ParamNum)
    CovarSetsDrawn = Array{Array{Float64,2}}(N-N_alpha)
    for ii = 1:length(PickParticleRand)
        #Return the index of the first value in CumSumWeights greater than or equal to PickParticleRand[ii]
        ParamDrawIdx = searchsortedfirst(CumSumWeights, PickParticleRand[ii])
        ParamSetsDrawn[ii,:] = RetainedParams[ParamDrawIdx,:]

        #If using local covariance, store selected covariance (associated with sample chosen to be perturbed)
        if ndims(C) == 3
          CovarSetsDrawn[ii] = C[:,:,ParamDrawIdx]
        end

    end
    #----------------------------------------------------------------------


    #----------------------------------------------------------------------
    # Generate new particles with the proposal distribution
    #----------------------------------------------------------------------
    PerturbedParamSets = SharedArray{Float64,2}(N-N_alpha,ParamNum)
    @sync @parallel for jj = 1:N-N_alpha #Resample out of range param. sets

        if ndims(C) == 1 #C is vector of component-wise variances
            d2=MvNormal(ParamSetsDrawn[jj,:],sqrt.(C))
        elseif ndims(C) == 2 #C is a global covariance
            d2=MvNormal(ParamSetsDrawn[jj,:],C)
        elseif ndims(C) == 3 #C is a local covariance
            d2=MvNormal(ParamSetsDrawn[jj,:],CovarSetsDrawn[jj])
        else
            error("Covariance array has a dimension above 3. Not compatible. Program ended.")
        end

        PerturbedParamSets[jj,:] = rand(d2)
    end

    #Check if perturbed parameter set is in prior bounds.
    PriorFnOutputTuple = PriorFn(PerturbedParamSets)
    InBoundFlag_NewParamSet = PriorFnOutputTuple[2]::Array{Int64,1}

    #Indicate if particle is contained within prior
    InBoundFlag_AllParamSet = [ones(Int64,N_alpha);InBoundFlag_NewParamSet]
    #----------------------------------------------------------------------

    #Assign perturbed parameter sets to concatenated paramter set array
    #PerturbedParamSets
    FullParamSet[(N_alpha + 1):N,:] = PerturbedParamSets
    @sync @parallel for ii = (N_alpha + 1):N
        #tic()
        #Assign current particle/param set to variable
        CurrentParamSet = FullParamSet[ii,:]::Array{Float64,1}

        #Assign validity of current particle to variable (does it lie in prior
        #dist?)
        InBoundFlag_CurrentParamSet = InBoundFlag_AllParamSet[ii]::Int64

        #If particle not in prior bounds, set weight to zero and SummStat value to
        #Inf. Otherwise, run model.
        if InBoundFlag_CurrentParamSet == 0
            SummStat[ii] = Inf
            Particle_Weights[ii] = 0
        else

            #Run model with designated parameter set
            SimnData = ModelSimnFn(FixedModelParams,CurrentParamSet)

            #Generate summary statistic for current prior set
            SummStatVal_NoPrior = SummStatFn(ObservedData,SimnData)

            if isnan(SummStatVal_NoPrior)     #Convert implausible (NaN) outputs to inf
                SummStat[ii] = Inf
            else
                SummStat[ii] = SummStatVal_NoPrior
            end

            #-------------------------------------------------------------------
            #Assign weights to particle
            #-------------------------------------------------------------------

            #Assign weight to particle - initialise count variables
            WeightsDenom = 0.
            RetainedWeightsSum = sum(RetainedWeights)
            NormWeights = RetainedWeights./RetainedWeightsSum

            #Assign weights to particle - variance variables
            if ndims(C) == 1 #C is vector of component-wise variances
                sigma = sqrt.(C)

                #Assign weights to particle - evaluation loop
                for jj = 1:N_alpha
                    d1 = Normal(RetainedParams[jj,1], sigma[1]) #Declare distribution
                    KernelVal = pdf(d1,CurrentParamSet[1]) #Get pdf for current parameter set
                    for kk=2:ParamNum
                        d1 = Normal(RetainedParams[jj,kk], sigma[kk])
                        KernelVal = KernelVal*pdf(d1,CurrentParamSet[kk])
                    end

                    WeightsDenom = WeightsDenom + (NormWeights[jj]*KernelVal)

                end
            elseif ndims(C) == 2 # C is a global variance-covariance array

                #Assign weights to particle - evaluation loop
                for jj = 1:N_alpha
                    d1 = MvNormal(RetainedParams[jj,:], C) #Declare distribution
                    KernelVal = pdf(d1,CurrentParamSet) #Get pdf for current parameter set

                    WeightsDenom = WeightsDenom + (NormWeights[jj]*KernelVal)
                end
            elseif ndims(C) == 3  # C is a local variance-covariance array
                #Assign weights to particle - evaluation loop
                for jj = 1:N_alpha
                    d1 = MvNormal(RetainedParams[jj,:], C[:,:,jj]) #Declare distribution
                    KernelVal = pdf(d1,CurrentParamSet) #Get pdf for current parameter set

                    WeightsDenom = WeightsDenom + (NormWeights[jj]*KernelVal)
                end
            else
                error("Covariance array has a dimension above 3. Not compatible. Program ended.")
            end

            #For particles lying outside prior range, amend weights to be zero
            PriorVal::Float64, InBoundFlagCheck::Int64 = PriorFn(CurrentParamSet)

            Particle_Weights[ii] = PriorVal/WeightsDenom
            #-------------------------------------------------------------------
        end

    end


    #Update acceptance rate
    AcceptRate = ((sum(SummStat[N_alpha+1:end].<epsilon))/(N-N_alpha))::Float64

    #Find particle idxs that are to be retained
    EpsilonIdx = sortperm(SummStat)

    ParamRetained = EpsilonIdx[1:N_alpha]

    #Obtain those particles that were present at beginnng of generation,
    #and have "survived" (been retained) at end of generation
    SurvivorParticleIdx = ParamRetained.<=N_alpha

    #Obtain next-generation of accepted particles
    RetainedParams = FullParamSet[ParamRetained,:]
    RetainedWeights = Particle_Weights[ParamRetained]
    RetainedSummStat = SummStat[ParamRetained]

    #Return the first alpha-quartile of SummStat
    epsilon = RetainedSummStat[end]

    #Set covariance matrix from accepted samples
    C = PerturbVarFn(RetainedParams,RetainedWeights,SurvivorParticleIdx,GenT)

    #Reorder FullParamSet, Weights, SummStat. Place retained values before unretained
    FullParamSet[1:N_alpha,:] = RetainedParams
    Particle_Weights[1:N_alpha] = RetainedWeights
    SummStat[1:N_alpha,:] = RetainedSummStat

    #-------------------------------------------------------------------------------
    # Write param values, weights, etc to output files
    #-------------------------------------------------------------------------------
    io1 = open(ParamSetFileName, "a")
    io2 = open(WeightsFileName, "a")
    io3 = open(SummStatFileName, "a")
    io4 = open(ThresholdValFileName, "a")
    io5 = open(AccRateFileName, "a")
    io6 = open(SurvivorParticleIdxFileName, "a")

    writedlm(io1,RetainedParams)
    writedlm(io2,RetainedWeights)
    writedlm(io3,RetainedSummStat)
    writedlm(io4,epsilon)
    writedlm(io5,AcceptRate)
    writedlm(io6,SurvivorParticleIdx)

    close(io1); close(io2); close(io3); close(io4); close(io5); close(io6);

    #-------------------------------------------------------------------------------

    println("GenT number $GenT complete")
end

return [RetainedParams,RetainedWeights,RetainedSummStat,SurvivorParticleIdx,GenT]

end
