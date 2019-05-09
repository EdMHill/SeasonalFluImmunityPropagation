# Purpose:
# Functions for solving system of ODES and running seasonal influenza transmission model


#Model specifics:
# Deterministic
# Multi-strain, non-age structured influenza transmission dynamic model
# Multi-strain model. Four strains in total (A/H1N1, A/H3N2, two B lineages)
# Vaccination status included
# Exposure history included

# Author: Ed Hill
#-------------------------------------------------------------------------------

# Define ODE equations
function SeasonalFluModelODEs_LeakyVacc(t,pop,dPop,InfectionParam::Array{Float64,2},BirthParam,
    DeathParam,DayOfYear,ExpHistArray,
    vacc_per_day,LeakyVaccVar::Array{Float64},LeakyTransFlag,NumOfStrains)

    #Inputs
    # infection_param - Array. Row per strain. Columns as follows:
    #    --> beta (Col 1) - Transmission rate, sigma (Col 2) - rate of loss of latency, gamma (Col 3) - rate of loss of infectiousness
    # BirthParam/DeathParam - Birth&death parameter
    # DayofYear - Actual day of year, with 0-1 corresponding to Jan 1st, 1-2 Jan 2nd etc
    # ExpHistArray - Two dimensional array. Column per exposure history. Row per strain.
    #                Interaction between exposure history and susceptibility to
    #                  the current season strain variant.
    # vacc_per_day - rate of immunisation
    # LeakyVaccVar - Two dimensional array. Row per strain.
    #                  Two columns per row:
    #                   --> Column 1 - Vaccination efficacy (on susceptibility),
    # LeakyTransFlag - Flag variable to specify if infectiousness of infected
    #                  vacc. groups is reduced (0 - Not reduced, 1 - Reduced)
    # NumOfStrains - Number of strains

    #Outputs
    # dPop_vec - Return values of derivatives at time t for values pop
        #Disaggregate infection param
        beta = InfectionParam[:,1] #Transmission rate
        sigma = InfectionParam[:,2] #rate of loss of latency
        gamma = InfectionParam[:,3] #rate of loss of infectiousness

        #Disaggregate DemogParam
        B = BirthParam #Birth rate
        D = DeathParam #Death rate

        #Get current day of year, may be used for seasonality, vaccine rate etc
        currentday = Int(ceil(mod(DayOfYear,365))) #Cast as integer, otherwise will be float & throw error
        if currentday == 0
            currentday = 365 #Deal with case where at precise end of calendar year
        end
        v = vacc_per_day[currentday] #Get vaccine rate for current day

        #Get number of exposure histories to be accounted for
        ExpHistNum = size(ExpHistArray,2)::Int64

        #Get susceptible compartments
        S_NotV=pop[1:ExpHistNum]
        S_V=pop[ExpHistNum+1:ExpHistNum*2]

        #Initialise numbers in each class and ODE variables
        SuscepClassEndIdx=ExpHistNum*2
        E_NotV=pop[SuscepClassEndIdx+1:SuscepClassEndIdx+NumOfStrains]
        I_NotV=pop[SuscepClassEndIdx+NumOfStrains+1:SuscepClassEndIdx+(2*NumOfStrains)]
        R_NotV=pop[SuscepClassEndIdx+(2*NumOfStrains)+1:SuscepClassEndIdx+(3*NumOfStrains)]
        E_V=pop[SuscepClassEndIdx+(3*NumOfStrains)+1:SuscepClassEndIdx+(4*NumOfStrains)]
        I_V=pop[SuscepClassEndIdx+(4*NumOfStrains)+1:SuscepClassEndIdx+(5*NumOfStrains)]
        R_V=pop[SuscepClassEndIdx+(5*NumOfStrains)+1:SuscepClassEndIdx+(6*NumOfStrains)]


        #Initialise ODE indexing variables
        SusODENum = 2*ExpHistNum

        #Get proportion unvaccinated
        UnvaccPropn = sum(S_NotV) + sum(E_NotV) + sum(I_NotV) + sum(R_NotV)

        #Rate of immunisation
        mu = v[1]/UnvaccPropn

        #--------------------------------------------------------------------------
        ### LEAKY VACCINE VARIABLES
        #Disaggragte LeakyVaccVar, scale force of infection if appropriate
        if LeakyTransFlag == 0 #Unmodfiied infectiousness
            alpha = LeakyVaccVar[:,1]
            I = I_NotV + I_V #Get total number of infecteds. Vector, entry per strain
        elseif LeakyTransFlag == 1 #Reduced infectiousness active
            alpha = LeakyVaccVar[:,1]
            delta = LeakyVaccVar[:,2]
            I = I_NotV + (1.0-delta).*I_V #Get scaled force of infection. Vector, entry per strain
        else
            error("Invalid value of $LeakyTransFlag for LeakyTransFlag, should be 0 or 1")
        end

        #---------------------------------------------------------------------------
        ### THE ODES
        ### Iterate through each set of exposure history and
        ### associated susceptible compartments
        for i=1:ExpHistNum

            #Note: In below, sum(beta.*I.*ExpHistArray[:,i]) corresponds to force of infection
            #from all strains against susceptibles

            #Non-vaccinated population, S compartment ODE
            dPop[i] = B[i] - sum(beta.*I.*ExpHistArray[:,i])*S_NotV[i] - D*S_NotV[i] - mu*S_NotV[i]

            #Vaccinated population, S compartment ODE
            dPop[ExpHistNum+i] = - sum((1.0-alpha).*beta.*I.*ExpHistArray[:,i])*S_V[i] - D*S_V[i] + mu*S_NotV[i]
        end

        ### Iterate through each set of strains
        for i = 1:NumOfStrains
            #---------------------------------------------------------------------------
            ### ODE for each non-susceptible, non-vaccinated compartment (E_NotV, I_NotV, R_NotV)
            dPop[SusODENum+i]= beta[i]*sum(ExpHistArray[i,:].*S_NotV)*I[i] - sigma[i]*E_NotV[i] - D*E_NotV[i] - mu*E_NotV[i]
            dPop[SusODENum+NumOfStrains+i]= sigma[i]*E_NotV[i] - gamma[i]*I_NotV[i] - D*I_NotV[i] - mu*I_NotV[i]
            dPop[SusODENum+(2*NumOfStrains)+i]= gamma[i]*I_NotV[i] - D*R_NotV[i] - mu*R_NotV[i]

            #---------------------------------------------------------------------------
            ### ODE for each non-susceptible, vaccinated compartment (E_V, I_V, R_V)
            dPop[SusODENum+(3*NumOfStrains)+i]= beta[i]*sum(ExpHistArray[i,:].*S_V)*(1.0-alpha[i])*I[i] - sigma[i]*E_V[i] - D*E_V[i] + mu*E_NotV[i]
            dPop[SusODENum+(4*NumOfStrains)+i]= sigma[i]*E_V[i] - gamma[i]*I_V[i] - D*I_V[i] + mu*I_NotV[i]
            dPop[SusODENum+(5*NumOfStrains)+i]= gamma[i]*I_V[i] - D*R_V[i] + mu*R_NotV[i]

            #---------------------------------------------------------------------------
            ### ODE for cumulative number of cases
            dPop[SusODENum+(6*NumOfStrains)+i]= sigma[i]*(E_NotV[i] + E_V[i])
        end
end


function RunSeasonalFluModel(SimnParam,
    InfectionParam,BirthParam,DeathParam,VaccUptakePerSeason,LeakyVaccVarBySeason,LeakyTransFlag,
    ExpHistArrayFnInputs,ExpHistNum,ExpHistVaccType,MultiSeasonImmPropn,
    SeasonStartPropnInf,ICFromFile,SimnRunType)
    #Inputs
    # SimnParam - Specify month of year simulation will begin (1-Jan, 2-Feb,..., 12-Dec),
    #             simn length and ode45 timestep value
    # infection_param - Array. Row per strain. Columns as follows:
    #    --> beta (Col 1) - Transmission rate, sigma (Col 2) - rate of loss of latency, gamma (Col 3) - rate of loss of infectiousness
    # BirthParam/DeathParam - Birth&death parameter
    # VaccUptakePerSeason - rate of immunisation, daily uptake rates
    # LeakyVaccVarBySeason - For exploratory runs (SimnRunType == 1)Two dimensional array. Row per strain.
    #                  Two columns per row:
    #                   --> Column 1 - Vaccination efficacy (on susceptibility),
    #                   --> Column 2 - Infectiousness reduction for those vaccinated
    #                    - For historical/forward simulation runs (SimnRunType == 2/3) Two entry cell.
    #                   --> Cell 1 - Vaccination efficacy (on susceptibility),
    #                   --> Cell 2 - Infectiousness reduction for those vaccinated
    #                   --> Col per strain. --> Season per row
    # LeakyTransFlag - Flag variable to specify if infectiousness of infected
    #                  vacc. groups is reduced (0 - Not reduced, 1 - Reduced)
    # ExpHistArrayFnInputs - Two entry cell.
    #                   --> Cell 1 - 2D array. Column per exposure history. Row per strain.
    #                Interaction between exposure history and susceptibility to
    #                  the current season strain variant.
    #                  --> Cell 2 - Exposure history params used to populate
    #                  ExpHistArray
    # ExpHistNum - Number of exposure history classes in use
    # ExpHistVaccType - Exposure history susceptibility modification form for vaccine-related states
    # MultiSeasonImmPropn - Constant: Proportion of population that are susceptible at end of season within
    #                       a natural infection exposure history class that remain in that expsoure
    #                       history category (rather than move to the Naive exposure history class)
    # SeasonStartPropnInf - Vector: % of population who are initialised as infected with each strain at
    #                       beginning of each flu season
    # ICFromFile - Indicate whether initial conditions will be obtained from
    #               external file and, if so, the file name
    # SimnRunType - Type of run: 1 - exploratory; 2 - historical; 3 - alternative vacc. schemes

    #Outputs
    # Compartment counts at each timestep

    #---------------------------------------------------------------------------
    ### SIMULATION RELATED VARIABLES

    #Disaggregate simulation paramter input vector
    SimnStartDate=convert(Int64,SimnParam[1])
    ODEBurnIn=SimnParam[2]
    ODESampleTime=SimnParam[3]
    MaxTime=SimnParam[4]
    timestep=SimnParam[5]
    TotalTime_H1N1Only = SimnParam[6]
    StoreFlag_PopnFOI = SimnParam[7]

    #Values to pass to tspan in ODE fn
    DaysPerMonth=[31 28 31 30 31 30 31 31 30 31 30 31]
    CumulDaysCount=zeros(length(DaysPerMonth))
    CumulDaysCount[2:end]=cumsum(DaysPerMonth[1:end-1])
    #------------------------------------------------------------------------------
    ###  SET UP INITIAL CONDITIONS

    #Get number of strains to be accounted for
    NumOfStrains = size(SeasonStartPropnInf,2)

    #Set up initial infected, pick out specified row from SeasonStartPropnInf
    #Row 1 for when only H1N1 infection allowed
    #Row 2 when infection by any strain allowed
    if TotalTime_H1N1Only > 0 #Initially, only H1N1 strain infection allowed
        I0_NotV = SeasonStartPropnInf[1,:]
    else #Can immediately allow infection by all strain types
        I0_NotV = SeasonStartPropnInf[2,:]
    end

    #Set initial conditions for ODEs
    #Initialise time and storage arrays for each set of compartments
    ICFromFileFlag = ICFromFile[1]::Array{Int64,1} #Specify if initial conditions should be obtained from file
    if ICFromFileFlag == 1
        S0_NotV  = ICFromFile[2] #IF POSSIBLE, IMPORT WITH INITIAL INF PROPN ALREADY ACCOUNTED FOR!
    else #All suscep. are put in naive, no exposure in previous year class
        S0_NotV = zeros(ExpHistNum)
        S0_NotV[1] = 1 - sum(I0_NotV) #All susceptibles assigned to naive group
    end

    #All other compartments begin empty
    S0_V = zeros(ExpHistNum)
    E0_NotV = zeros(NumOfStrains)
    R0_NotV = zeros(NumOfStrains)

    E0_V = zeros(NumOfStrains)
    I0_V = zeros(NumOfStrains)
    R0_V = zeros(NumOfStrains)
    C0 = zeros(NumOfStrains)

    #Initialise time and storage arrays for each set of compartments
    T0=0.0
    NumYrsRecorded = ceil(ODESampleTime/365)
    CompartmentTraceSize = (ODESampleTime/timestep) + NumYrsRecorded #Add NumYrsRecorded to include initial conditions at start of each flu season
    CompartmentTraceSize = convert(Int64,CompartmentTraceSize)

    NumYrsSimTotal = ceil(MaxTime/365)
    FOI_CompartmentTraceSize = (MaxTime/timestep) + NumYrsSimTotal #Add NumYrsSimTotal to include initial condition
    FOI_CompartmentTraceSize = convert(Int64,FOI_CompartmentTraceSize)

    Store_S_NotV = zeros(CompartmentTraceSize,ExpHistNum)
    Store_S_V = zeros(CompartmentTraceSize,ExpHistNum)
    Store_E_NotV = zeros(CompartmentTraceSize,NumOfStrains) #Row per timestep, column per strain!
    Store_I_NotV = zeros(CompartmentTraceSize,NumOfStrains)
    Store_R_NotV = zeros(CompartmentTraceSize,NumOfStrains)
    Store_E_V = zeros(CompartmentTraceSize,NumOfStrains)
    Store_I_V = zeros(CompartmentTraceSize,NumOfStrains)
    Store_R_V = zeros(CompartmentTraceSize,NumOfStrains)
    Store_C = zeros(CompartmentTraceSize,NumOfStrains)  #Cumulative number of cases
    Store_T = zeros(CompartmentTraceSize)

    #If required, initialise storage array for population-level FOI
    if StoreFlag_PopnFOI == 1
       Store_PopnFOI = zeros(FOI_CompartmentTraceSize,NumOfStrains)
    end

    #---------------------------------------------------------------------------
    ### INITIALISE VACCINE UPTAKE AND EFFICACY VALUES (BASED ON SIMN. RUN TYPE)
    if SimnRunType == 1
        #SYNTHETIC DATA, no amendments to uptake/efficacy data needed
        vacc_per_day = VaccUptakePerSeason
        LeakyVaccVar = LeakyVaccVarBySeason
    elseif SimnRunType == 2
        #BASED ON HISTORICAL DATA
        if TotalTime_H1N1Only > 0
            #Initially use 2009/2010 data (first row of arrays)
            RecordedSeasonIdx = 1
        else
            #If no H1N1 only period, use 2010/2011 data (second row of vacc. arrays)
            RecordedSeasonIdx = 2
        end

        vacc_per_day = VaccUptakePerSeason[RecordedSeasonIdx,:]

        if LeakyTransFlag == 0
            LeakyVaccVar = LeakyVaccVarBySeason[:,RecordedSeasonIdx]
        elseif LeakyTransFlag == 1
            alpha = LeakyVaccVarBySeason[1][:,RecordedSeasonIdx]
            delta = LeakyVaccVarBySeason[2][:,RecordedSeasonIdx]
            LeakyVaccVar = [alpha delta]
        end
    elseif SimnRunType == 3
        #RUN ALTERNATIVE VACC. SCHEMES
        if TotalTime_H1N1Only > 0
            #Initially use 2009/2010 data (first row of arrays)
            RecordedSeasonIdx = 1
        else
            RecordedSeasonIdx = 2
        end

        vacc_per_day = VaccUptakePerSeason[RecordedSeasonIdx,:]

        if LeakyTransFlag == 0
            LeakyVaccVar = LeakyVaccVarBySeason[:,RecordedSeasonIdx]
        elseif LeakyTransFlag == 1
            alpha = LeakyVaccVarBySeason[1][:,RecordedSeasonIdx]
            delta = LeakyVaccVarBySeason[2][:,RecordedSeasonIdx]
            LeakyVaccVar = [alpha delta]
        end
    else
        error("Unknown SimnRunType value, %d",SimnRunType)
    end

    #Get ExpHistArray
    ExpHistArray = ExpHistArrayFnInputs[1]
    ExpHistArrayParams = ExpHistArrayFnInputs[2]

    ExpHistArray = ExpHistUpdate(ExpHistVaccType,ExpHistArray,ExpHistArrayParams,LeakyVaccVar[:,1],
                                                ExpHistNum,NumOfStrains)

    #------------------------------------------------------------------------------
    ### MAIN MODEL LOOP

    #Constrain S,E,I,R components to be nonnegative (if using ode45!)
    ODE_EqnNum_Sus = ExpHistNum*2
    ODE_EqnNum_NonSus = 7*NumOfStrains
    ODE_TotalEqnNum = ODE_EqnNum_NonSus+ODE_EqnNum_Sus #Get total number of ODEs in system

    #Initialise initial condition vector
    IC = zeros(ODE_TotalEqnNum)
    IC[1:ExpHistNum] =  S0_NotV
    IC[ExpHistNum+1:ODE_EqnNum_Sus] =  S0_V
    IC[ODE_EqnNum_Sus+1:end] =  [E0_NotV I0_NotV R0_NotV E0_V I0_V R0_V C0]

    #Specify current month
    MonthIdx = SimnStartDate

    InitialStoreFlag = 0 #Flag to specify if any values have been stored yet
    StoreStartIdx = 1 #Initialise row indexing counter for storage arrays

    InitialStoreFlag_PopnFOI = 0 #Flag to specify if any popnFOI values have been stored yet
    StorePopnFOI_StartIdx = 1 #Initialise row indexing counter for popnFOI storage arrays

    # Initialise model
    while T0<MaxTime
        #tspan = (0.0,365.0)
        tspanForSol = (T0,T0+DaysPerMonth[MonthIdx]) #Specify time bounds to used in ODE solver
        tspan = T0:timestep:T0+DaysPerMonth[MonthIdx] #Timestep outputs

        #---------------------------------------------------------------------------
        ### INITIAL CONDITIONS

        #Error check
        if any(isnan,IC)
            error("NaN found in initial conditions")
        end
        #Normalise IC so it sums to one!
        IC[1:end-NumOfStrains] = IC[1:end-NumOfStrains]/sum(IC[1:end-NumOfStrains])

        #---------------------------------------------------------------------------
        ### RUN ODE SOLVER
          prob = ODEProblem((t,pop,dpop) ->SeasonalFluModelODEs_LeakyVacc(t,pop,dpop,
          InfectionParam,BirthParam,DeathParam,t+CumulDaysCount[SimnStartDate],
          ExpHistArray,vacc_per_day,LeakyVaccVar,LeakyTransFlag,NumOfStrains),IC,tspanForSol)

          if MonthIdx == SimnStartDate #If start of flu season, store initial value, save_start = true
              sol = solve(prob, Tsit5(), saveat=timestep, abstol = 1e-8, reltol = 1e-10,
              isoutofdomain=(t,u) -> any(x -> x < 0, u))
          else
              sol = solve(prob, Tsit5(), saveat=timestep, save_start = false, abstol = 1e-8, reltol = 1e-10,
              isoutofdomain=(t,u) -> any(x -> x < 0, u))
          end


         #---------------------------------------------------------------------------
         ### PERFORM STORAGE OF COMPARTMENT VALUES
         SusClassTotal=ExpHistNum #number of equations for S^N or S^V
         NonSusClassTotal = NumOfStrains #number of equations per E,I,R,C state

         # Amalgamate output in arrays for current time window
         T = tspan'


         # Collate into Array
        OutputNum=length(sol)::Int64
        pop = zeros(OutputNum,ODE_TotalEqnNum)
        for k = 1:OutputNum
            pop[k,:] = sol[k]
        end
         #Separate columns into susceptible and non-susceptible compartments
         SusNotVPop = pop[:,1:SusClassTotal]::Array{Float64,2}
         SusVPop = pop[:,SusClassTotal+1:(SusClassTotal*2)]::Array{Float64,2}
         ENotV_Pop = pop[:,SusClassTotal*2+1:(SusClassTotal*2) + NonSusClassTotal]::Array{Float64,2}
         INotV_Pop = pop[:,(SusClassTotal*2)+NonSusClassTotal+1:(SusClassTotal*2) + 2*NonSusClassTotal]::Array{Float64,2}
         RNotV_Pop = pop[:,(SusClassTotal*2)+(2*NonSusClassTotal)+1:(SusClassTotal*2) + 3*NonSusClassTotal]::Array{Float64,2}
         EV_Pop = pop[:,(SusClassTotal*2)+(3*NonSusClassTotal)+1:(SusClassTotal*2) + 4*NonSusClassTotal]::Array{Float64,2}
         IV_Pop = pop[:,(SusClassTotal*2)+(4*NonSusClassTotal)+1:(SusClassTotal*2) + 5*NonSusClassTotal]::Array{Float64,2}
         RV_Pop = pop[:,(SusClassTotal*2)+(5*NonSusClassTotal)+1:(SusClassTotal*2) + 6*NonSusClassTotal]::Array{Float64,2}
         C_Pop = pop[:,(SusClassTotal*2)+(6*NonSusClassTotal)+1:(SusClassTotal*2) + 7*NonSusClassTotal]::Array{Float64,2}

         #If beyond burn in, store non-susceptible compartment values
         if T[end]>ODEBurnIn
             if InitialStoreFlag == 0
                 StoreEndIdx = length(tspan)
                 Store_T[StoreStartIdx:StoreEndIdx] = tspan'
                 InitialStoreFlag = 1 #Amend flag value so will not pass through this loop again!
             elseif MonthIdx == SimnStartDate #Store initial condition if at start of flu season
                 StoreEndIdx = StoreStartIdx + length(tspan) - 1
                 #Need to realign indexing. Adding on length(tspan) makes index 1 too large

                 Store_T[StoreStartIdx:StoreEndIdx] = tspan'
             else
                 StoreEndIdx = StoreStartIdx + length(tspan) - 2
                 #Subtract 2 as will not store initial condition,
                 Store_T[StoreStartIdx:StoreEndIdx] = tspan[2:end]'
             end

             #Store susceptible compartment values
            Store_S_NotV[StoreStartIdx:StoreEndIdx,:] = SusNotVPop
            Store_S_V[StoreStartIdx:StoreEndIdx,:] = SusVPop

             #Store non-susceptible compartment values
             Store_E_NotV[StoreStartIdx:StoreEndIdx,:] = ENotV_Pop
             Store_I_NotV[StoreStartIdx:StoreEndIdx,:] = INotV_Pop
             Store_R_NotV[StoreStartIdx:StoreEndIdx,:] = RNotV_Pop
             Store_E_V[StoreStartIdx:StoreEndIdx,:] = EV_Pop
             Store_I_V[StoreStartIdx:StoreEndIdx,:] = IV_Pop
             Store_R_V[StoreStartIdx:StoreEndIdx,:] = RV_Pop
             Store_C[StoreStartIdx:StoreEndIdx,:] = C_Pop

             StoreStartIdx = StoreEndIdx + 1 #Update StartIdx value
         end

         #----------------------------------------------------------------------
         #Store total number of infected, can then be used to recover risk group
         #information!

         if StoreFlag_PopnFOI == 1
             if InitialStoreFlag_PopnFOI == 0
                 StorePopnFOI_EndIdx = length(tspan)

                 #Scale force of infection if appropriate
                 if LeakyTransFlag == 0 #Unmodfiied infectiousness
                     Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:] =
                         INotV_Pop + IV_Pop #Get total number of infecteds. Vector, entry per strain
                 elseif LeakyTransFlag == 1 #Reduced infectiousness active
                     delta = LeakyVaccVar[:,2]
                     Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:] = INotV_Pop + (1-delta').*IV_Pop
                     #Get scaled force of infection. Vector, entry per strain
                 else
                     error("Invalid value of %d for LeakyTransFlag, should be 0 or 1", LeakyTransFlag)
                 end

                 InitialStoreFlag_PopnFOI = 1 #Amend flag value so will not pass through this loop again!
              elseif MonthIdx == SimnStartDate #Store initial condition if at start of flu season
                 StorePopnFOI_EndIdx = StorePopnFOI_StartIdx + length(tspan) - 1
                 #Need to realign indexing. Adding on length(tspan) makes index 1 too large

                 #Scale force of infection if appropriate
                 if LeakyTransFlag == 0 #Unmodfiied infectiousness
                     Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:] =
                         INotV_Pop + IV_Pop #Get total number of infecteds. Vector, entry per strain
                 elseif LeakyTransFlag == 1 #Reduced infectiousness active
                     delta = LeakyVaccVar[:,2]
                     Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:] = INotV_Pop + (1-delta').*IV_Pop
                     #Get scaled force of infection. Vector, entry per strain
                 else
                     error("Invalid value of %d for LeakyTransFlag, should be 0 or 1", LeakyTransFlag)
                 end

             else
                 StorePopnFOI_EndIdx = StorePopnFOI_StartIdx + length(tspan) - 2 #Subtract 2 as will not store initial condition,

                 #Scale force of infection if appropriate
                 if LeakyTransFlag == 0 #Unmodfiied infectiousness
                     Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:] =
                         INotV_Pop+ IV_Pop #Get total number of infecteds. Vector, entry per strain
                 elseif LeakyTransFlag == 1 #Reduced infectiousness active
                     delta = LeakyVaccVar[:,2]
                     Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:] = INotV_Pop + (1-delta').*IV_Pop
                     #Get scaled force of infection. Vector, entry per strain
                 else
                     error("Invalid value of %d for LeakyTransFlag, should be 0 or 1", LeakyTransFlag)
                 end
             end
             StorePopnFOI_StartIdx = StorePopnFOI_EndIdx + 1 #Update StartIdx value
         end

         #---------------------------------------------------------------------------
         ### PERFORM END OF MONTH/SEASON MOVEMENTS
         #Store susceptible compartment values
         S_NotV = SusNotVPop[end,:]::Vector{Float64}
         S_V = SusVPop[end,:]::Vector{Float64}

         E_NotV = ENotV_Pop[end,:]::Vector{Float64}
         I_NotV= INotV_Pop[end,:]::Vector{Float64}
         R_NotV = RNotV_Pop[end,:]::Vector{Float64}
         E_V = EV_Pop[end,:]::Vector{Float64}
         I_V = IV_Pop[end,:]::Vector{Float64}
         R_V = RV_Pop[end,:]::Vector{Float64}
         C = C_Pop[end,:]::Vector{Float64}

         if MonthIdx == 8  #End of season, reinitialise infecteds
            #First sum over the susceptible compartments
             S_Unvacc_MultiSeasonImmEligable = sum(S_NotV[2:NumOfStrains+1]) + sum(S_NotV[NumOfStrains+3:end])

             #vaccinated individuals, fraction to maintain previous infection
             #history!
             S_Vacc_MultiSeasonImmEligable = sum(S_V[2:NumOfStrains+1]) + sum(S_V[NumOfStrains+3:end])

            #---------------------------------------------------------------------------
            ### MAPPINGS TO EACH EXPOSURE HISTORY CLASS


                #No exposure through infection/succesful vacc. (in effect, sum over S^N classes)
                #Account for any "foldover" residual immunity due to natural infection
                S0_NotV[1] = S_NotV[1] + S_NotV[6] + ((1-MultiSeasonImmPropn)*S_Unvacc_MultiSeasonImmEligable)

                #Natural infection by strain i, no vaccination: E_NotV(i) + I_NotV(i) + R_NotV(i) term
                #Account for "foldover" residual immunity due to natural infection: MultiSeasonImmPropn*S_NotV term
                for ii = 1:NumOfStrains
                    S0_NotV[1+ii] = E_NotV[ii] + I_NotV[ii] + R_NotV[ii] + (MultiSeasonImmPropn*(S_NotV[1+ii] + S_NotV[6+ii]))
                end

                #Vaccination, with no natural infection
                #Account for any "foldover" residual immunity due to
                # previous natural infection: (1-MultiSeasonImmPropn)*S_Vacc_MultiSeasonImmEligable
                S0_NotV[6] = S_V[1] + S_V[6] + ((1-MultiSeasonImmPropn)*S_Vacc_MultiSeasonImmEligable)

                #Natural infection by strain i, and vaccinated
                #Account for any "foldover" residual immunity due to natural
                #infection: MultiSeasonImmPropn*S_V(7:10)
                for ii = 1:NumOfStrains
                    S0_NotV[6+ii] = E_V[ii] + I_V[ii] + R_V[ii] + (MultiSeasonImmPropn*(S_V[1+ii]+S_V[6+ii]))
                end

            #-----------------------------------------------------------------------
            ### Initialise infectious popn, update susceptible compartments

            #Row 1 for when only H1N1 infection allowed
            #Row 2 when infection by any strain allowed
            if T[end]<TotalTime_H1N1Only #Initially, only H1N1 strain infection allowed
                I0_NotV = SeasonStartPropnInf[1,:]
            else #Can immediately allow infection by all strain types
                I0_NotV = SeasonStartPropnInf[2,:]
            end

            TotalSeasonStartPropnInf = sum(I0_NotV) #Proportion of popn initially infected
            RelMassPerExpHist=S0_NotV/sum(S0_NotV[:]) #Density in each exposure history Sus. compartment


            #Modify susceptibles assigned to exposure history j
            #Split initial infectious weighted by exposure history
            S0_NotV = S0_NotV - (RelMassPerExpHist*TotalSeasonStartPropnInf)

            #-----------------------------------------------------------------------
            ### Reset all IC for E and R non-vaccination compartments, and all Vacc. compartments
            E0_NotV = zeros(NumOfStrains)
            R0_NotV = zeros(NumOfStrains)

            S0_V = zeros(ExpHistNum) #Reset S^V class to zeros
            E0_V = zeros(NumOfStrains)
            I0_V = zeros(NumOfStrains)
            R0_V = zeros(NumOfStrains)

            #Keep tracking cumulative number of infected, no reset required
            C0 = C

            #Construct initial conditions vector
            IC[1:ExpHistNum] =  S0_NotV
            IC[ExpHistNum+1:ODE_EqnNum_Sus] =  S0_V
            IC[ODE_EqnNum_Sus+1:end] =  [E0_NotV I0_NotV R0_NotV E0_V I0_V R0_V C0]

            #---------------------------------------------------------------------------
            ### IF BEYOND TotalTime_H1N1Only, UPDATE VACCINE UPTAKE AND EFFICACY VALUES (BASED ON SIMN. RUN TYPE)
            if TotalTime_H1N1Only == 0
                ExpArrayAmendTimeThreshold = ODEBurnIn
            else
                ExpArrayAmendTimeThreshold = TotalTime_H1N1Only
            end

            if T[end]>=ExpArrayAmendTimeThreshold && SimnRunType!=1
                #Note, SYNTHETIC DATA, no amendments to uptake/efficacy data needed

                ExpHistArray = ExpHistUpdate(ExpHistVaccType,ExpHistArray,ExpHistArrayParams,LeakyVaccVar[:,1],
                                                            ExpHistNum,NumOfStrains)

                #Update Idx to pick out relevant season data from array
                if RecordedSeasonIdx < size(VaccUptakePerSeason,1)
                    RecordedSeasonIdx = RecordedSeasonIdx + 1
                end
                if SimnRunType == 2 || SimnRunType == 3
                    vacc_per_day = VaccUptakePerSeason[RecordedSeasonIdx,:]::Array{Float64,1}

                    if LeakyTransFlag == 0
                        LeakyVaccVar = LeakyVaccVarBySeason[:,RecordedSeasonIdx]::Array{Float64,1}
                    elseif LeakyTransFlag == 1
                        alpha = LeakyVaccVarBySeason[1][:,RecordedSeasonIdx]::Array{Float64,1}
                        delta = LeakyVaccVarBySeason[2][:,RecordedSeasonIdx]::Array{Float64,1}
                        LeakyVaccVar = [alpha delta]::Array{Float64,2}
                    end
                else
                    error("Unknown SimnRunType value, %d",SimnRunType)
                end
            end
         else #End of month, but not end of season.
            #Set up initial conditions (final values in each state vector))
            IC = pop[end,:]::Array{Float64,1}
         end

         T0 = T[end]::Float64

         #---------------------------------------------------------------------------
         ### UPDATE MONTH INDEX
         MonthIdx = mod(MonthIdx + 1,12)
         if MonthIdx == 0
             MonthIdx = 12
         end
    end

    if StoreFlag_PopnFOI == 0
        return Store_T, Store_C, Store_S_NotV,Store_S_V,Store_E_NotV,Store_E_V,Store_I_NotV,Store_I_V,
        Store_R_NotV,Store_R_V
    else
        return Store_T, Store_C, Store_S_NotV,Store_S_V,Store_E_NotV,Store_E_V,Store_I_NotV,Store_I_V,
        Store_R_NotV,Store_R_V,Store_PopnFOI
     end
end
