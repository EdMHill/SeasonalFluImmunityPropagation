%Purpose:
%Compute posterior parameter distributions summary statistics,
%obtained using AMPC ABC inference scheme

%Author: Ed Hill
%--------------------------------------------------------------------------

%Add required directories to path
addpath('../MatlabFileExchangeScripts')

%% Get data
clear variables

%--------------------------------------------------------------------------
%%% DECLARE REQUESTED PERCENTILE VALUES
%--------------------------------------------------------------------------
ReqPrctiles = [2.5, 50, 97.5];

%--------------------------------------------------------------------------
%%% LOAD PARTICLE SETS & WEIGHTS: FOUR SEASON FIT
%--------------------------------------------------------------------------
PosteriorParams = dlmread('FourSeasonFitData/ParameterSets_FourSeasonFit.txt');
ParamWeights = dlmread('FourSeasonFitData/ParticleWeights_FourSeasonFit.txt');

OutputPrctileVals_FourSeasonFit = PostDistStats(PosteriorParams,ParamWeights,ReqPrctiles)

%--------------------------------------------------------------------------
%%% LOAD PARTICLE SETS & WEIGHTS: FIVE SEASON FIT
%--------------------------------------------------------------------------
PosteriorParams = dlmread('FiveSeasonFitData/ParameterSets_FiveSeasonFit.txt');
ParamWeights = dlmread('FiveSeasonFitData/ParticleWeights_FiveSeasonFit.txt');

OutputPrctileVals_FiveSeasonFit = PostDistStats(PosteriorParams,ParamWeights,ReqPrctiles)

%--------------------------------------------------------------------------
%%% LOAD PARTICLE SETS & WEIGHTS: SIX SEASON FIT
%--------------------------------------------------------------------------
PosteriorParams = dlmread('SixSeasonFitData/ParameterSets_SixSeasonFit.txt');
ParamWeights = dlmread('SixSeasonFitData/ParticleWeights_SixSeasonFit.txt');

OutputPrctileVals_SixSeasonFit = PostDistStats(PosteriorParams,ParamWeights,ReqPrctiles)

%--------------------------------------------------------------------------
%%% LOAD PARTICLE SETS & WEIGHTS: SIX SEASON FIT 
%%%(INCLUDING IMMUNITY PROPAGATION)
%--------------------------------------------------------------------------
PosteriorParams = dlmread('SixSeasonFitData/ExtendedModel_ParameterSets_SixSeasonFit.txt');
ParamWeights = dlmread('SixSeasonFitData/ExtendedModel_ParticleWeights_SixSeasonFit.txt');

OutputPrctileVals_SixSeasonFitIncImmProp = PostDistStats(PosteriorParams,ParamWeights,ReqPrctiles)

%--------------------------------------------------------------------------
%%% LOAD PARTICLE SETS & WEIGHTS: SYNTHETIC DATA SIX SEASON FIT
%--------------------------------------------------------------------------
PosteriorParams = dlmread('SynthFitData/ParameterSets_SynthDataFit.txt');
ParamWeights = dlmread('SynthFitData/ParticleWeights_SynthDataFit.txt');

OutputPrctileVals_SynthData = PostDistStats(PosteriorParams,ParamWeights,ReqPrctiles)

%--------------------------------------------------------------------------
%%% LOAD PARTICLE SETS & WEIGHTS: SYNTHETIC DATA SIX SEASON FIT 
%%%(INCLUDING IMMUNITY PROPAGATION)
%--------------------------------------------------------------------------
PosteriorParams = dlmread('SynthFitData/ExtendedModel_ParameterSets_SynthDataFit.txt');
ParamWeights = dlmread('SynthFitData/ExtendedModel_ParticleWeights_SynthDataFit.txt');

OutputPrctileVals_SynthDataIncImmProp = PostDistStats(PosteriorParams,ParamWeights,ReqPrctiles)


function OutputPrctileVals = PostDistStats(PosteriorParams,ParamWeights,ReqPrctiles)
%Inputs
%PosteriorParams - (2D Array) Row per simulation run. Column per parameter.
%ParamWeights - (Column vector)
%ReqPrctiles - (1D vector)

%Outputs
%OutputPrctileVals - (2D array) Row per requested prctile value. Column per parameter.

%--------------------------------------------------------------------------
%%% VARIABLE INITIALISATION
%--------------------------------------------------------------------------

%Number of parameters (columns of StatInput)
ParamNum = size(PosteriorParams,2);

%Number of prctile evaluations requested
PrctileNum = numel(ReqPrctiles);

%Convert percentage to proportion
ReqPropn = ReqPrctiles/100;

%Initialise output array
OutputPrctileVals = zeros(PrctileNum,ParamNum);

%--------------------------------------------------------------------------
%%% NORMALISE WEIGHTS
%--------------------------------------------------------------------------
NormParamWeights = ParamWeights/sum(ParamWeights);

%--------------------------------------------------------------------------
%%% WEIGHTED PERCENTILE COMPUTATION
%--------------------------------------------------------------------------

%Sort each column of StatInput into ascending order
[SortedStatInput,SortIdx] = sort(PosteriorParams,1);

%For each parameter, compute cumulative sum of weights
CumSumWeightsByParam = cumsum(NormParamWeights(SortIdx),1);

%Iterate over requested prctile values.
for ii = 1:PrctileNum

    if ReqPropn(ii) < 0.5
        %------------------------------------------------------------------
        %%% PRCTILES BELOW 50 (use infimum. Greatest lower bound)
        %------------------------------------------------------------------
        
        for jj = 1:ParamNum           
              SelectedParamIdx = find(CumSumWeightsByParam(:,jj) <= ReqPropn(ii),1,'last');
              OutputPrctileVals(ii,jj) = SortedStatInput(SelectedParamIdx,jj); 
        end
        
    elseif ReqPropn(ii) == 0.5
        
        %------------------------------------------------------------------
        %%% MEDIAN
        %------------------------------------------------------------------
        for jj = 1:ParamNum        
            OutputPrctileVals(ii,jj) = weightedMedian(PosteriorParams(:,jj),ParamWeights);
        end
    else
        %------------------------------------------------------------------
        %%% PRCTILES ABOVE 50 (use supremum. Least upper bound)
        %------------------------------------------------------------------
        
        for jj = 1:ParamNum           
              SelectedParamIdx = find(CumSumWeightsByParam(:,jj) >= ReqPropn(ii),1,'first');
              OutputPrctileVals(ii,jj) = SortedStatInput(SelectedParamIdx,jj); 
        end
    end
end

end