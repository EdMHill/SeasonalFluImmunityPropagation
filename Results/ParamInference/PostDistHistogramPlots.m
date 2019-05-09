%Purpose:
%Plot posterior parameter distributions, obtained AMPC ABC inference
%scheme

%Author: Ed Hill
%--------------------------------------------------------------------------


%% Get data
clear variables

%--------------------------------------------------------------------------
% DEFINE FUNCTIONS INPUT VALUES
%--------------------------------------------------------------------------
SynthDataFlag = 0; %Indicator variable. 0: empirical data; 1: synthetic data


MultiSeasonImmFlag = 0;  %Model type used. 
                         %0: Original, immunity propagation lasts a single season; 
                         %1:Extended, allowing immunity to be propagated for more than a single season


%--------------------------------------------------------------------------
%%% DECLARE NUMBER OF SEASONS FITTED TO VARIABLE
%%% DATA FITTED FROM 2012/13 INFLUENZA SEASON ONWARDS
%%% TAKES VALUE OF 4 (2012/13 - 2015/16),5 (2012/13 - 2016/17) OR 6 (2012/13 - 2017/18)
%--------------------------------------------------------------------------
SeasonsFitted = 5;
       
%If want to plot outputs fitting to data with the extended model or syntehtic data, 
%only six seaon fit is available
%By deafult, revert SeasonsToPlot to 6
if SynthDataFlag == 1 || MultiSeasonImmFlag == 1
    SeasonsFitted = 6;
end

%--------------------------------------------------------------------------
% ERROR CHECKS
%--------------------------------------------------------------------------
if SynthDataFlag ~=0 && SynthDataFlag ~=1
    error('Misspecified SynthDataFlag value of %f. SynthDataFlag must take value 0 or 1.',SynthDataFlag)
end

if MultiSeasonImmFlag ~=0 && MultiSeasonImmFlag ~=1
    error('Misspecified MultiSeasonImmFlag value of %f. MultiSeasonImmFlag must take value 0 or 1.',MultiSeasonImmFlag)
end

if SeasonsFitted ~=4 && SeasonsFitted ~=5 && SeasonsFitted ~=6
    error('SeasonsToPlot must take a value of 4,5 or 6')
end


%--------------------------------------------------------------------------
% CALL FUNCTION
%--------------------------------------------------------------------------
ConstructPosteriorHistograms(SynthDataFlag,MultiSeasonImmFlag,SeasonsFitted)

function ConstructPosteriorHistograms(SynthDataFlag,MultiSeasonImmFlag,SeasonsFitted)
%   SynthDataFlag - (indicator variable) 0: empirical data; 1: synthetic data
%   MultiSeasonImmFlag - (indicator variable) Model type used. 
%                                             0: Original, immunity propagation lasts a single season; 
%                                             1: Extended, allowing immunity to be propagated for more than a single season
%   SeasonsToPlot - (scalar, integer) Number of influenza seasons worth of data fit to in the inference scheme 



%--------------------------------------------------------------------------
%%% SPECIFY SYNTHETIC DATA PARAMETER VALUES & FILE LOCATIONS OF DATA
%--------------------------------------------------------------------------
if SynthDataFlag == 1 %use outputs from fitting to synthetic data
    
    %Parameter values from which synthetic data were generated
    ParamTrueVal_beta = [0.42,0.43,0.415,0.405];
    ParamTrueVal_AscertainProb = [0.0015,0.0007,0.0023,0.0021,0.0004,0.003];
      
    if MultiSeasonImmFlag == 0
        ParamTrueVal_ExpHist = [0.5837,0.4351,0.7026];
    else
        ParamTrueVal_ExpHist = [0.5837,0.4351,0.7026,0.1];
    end
    
    ParamTrueVal = [ParamTrueVal_beta ParamTrueVal_ExpHist ParamTrueVal_AscertainProb];
    
    SynthDataVar = {SynthDataFlag,ParamTrueVal};
    
    %Files where inferred parameter sets are stored
    %Choose file based on whether original or extended model is
    %specified
    if MultiSeasonImmFlag == 0
        FileName = 'SynthFitData/ParameterSets_SynthDataFit.txt';
    else
        FileName = 'SynthFitData/ExtendedModel_ParameterSets_SynthDataFit.txt';
    end
    SeasonsFitted = 6; %By deafult, revert SeasonsFitted to 6

else %use outputs from fitting to empirical data
    ParamTrueVal = []; %Specify ParamTrueVal an empty variable
    SynthDataVar = {SynthDataFlag,ParamTrueVal}; 
    
    if SeasonsFitted == 4
        FileName = 'FourSeasonFitData/ParameterSets_FourSeasonFit.txt';
    elseif SeasonsFitted == 5
        FileName = 'FiveSeasonFitData/ParameterSets_FiveSeasonFit.txt';
    else
        %Choose file based on whether original or extended model is
        %specified
        if MultiSeasonImmFlag == 0
            FileName = 'SixSeasonFitData/ParameterSets_SixSeasonFit.txt';
        else
            FileName = 'SixSeasonFitData/ExtendedModel_ParameterSets_SixSeasonFit.txt';
        end
    end
end


%--------------------------------------------------------------------------
%%% LOAD PARTICLE SETS
%--------------------------------------------------------------------------
PosteriorParamsTemp = dlmread(FileName);

%Swap B/Yamagata and B/Victoria transmissibility param value columns. 
PosteriorParams = PosteriorParamsTemp;
PosteriorParams(:,3) = PosteriorParamsTemp(:,4); %B/Victoria, was column 4, move to column 3
PosteriorParams(:,4) = PosteriorParamsTemp(:,3); %B/Yamagata, was column 3, move to column 4

%--------------------------------------------------------------------------
%%% DECLARE PARAMETER LABELS
%--------------------------------------------------------------------------
if SeasonsFitted == 4
    x_label = {'\beta_{A(H1N1)pdm09}','\beta_{A(H3N2)}','\beta_{B/Victoria}','\beta_{B/Yamagata}',...
            'a (Nat. Inf. Mod. Sus.)','b (Inf. B cross-reactivity)','\xi',...
            '\epsilon_{2012/13}','\epsilon_{2013/14}','\epsilon_{2014/15}','\epsilon_{2015/16}'};
elseif SeasonsFitted == 5
    x_label = {'\beta_{A(H1N1)pdm09}','\beta_{A(H3N2)}','\beta_{B/Victoria}','\beta_{B/Yamagata}',...
            'a (Nat. Inf. Mod. Sus.)','b (Inf. B cross-reactivity)','\xi',...
            '\epsilon_{2012/13}','\epsilon_{2013/14}','\epsilon_{2014/15}','\epsilon_{2015/16}',...
            '\epsilon_{2016/17}'};
elseif SeasonsFitted == 6
    if MultiSeasonImmFlag == 0
        x_label = {'\beta_{A(H1N1)pdm09}','\beta_{A(H3N2)}','\beta_{B/Victoria}','\beta_{B/Yamagata}',...
            'a (Nat. Inf. Mod. Sus.)','b (Inf. B cross-reactivity)','\xi',...
            '\epsilon_{2012/13}','\epsilon_{2013/14}','\epsilon_{2014/15}','\epsilon_{2015/16}',...
            '\epsilon_{2016/17}','\epsilon_{2017/18}'};
    else
        x_label = {'\beta_{A(H1N1)pdm09}','\beta_{A(H3N2)}','\beta_{B/Victoria}','\beta_{B/Yamagata}',...
            'a (Nat. Inf. Mod. Sus.)','b (Inf. B cross-reactivity)','\xi','\delta',...
            '\epsilon_{2012/13}','\epsilon_{2013/14}','\epsilon_{2014/15}','\epsilon_{2015/16}',...
            '\epsilon_{2016/17}','\epsilon_{2017/18}'};
    end
end
       

%--------------------------------------------------------------------------
%%% SET UP PRIOR DISTRIBUTION PLOTTING VARIABLES
%--------------------------------------------------------------------------
ParamNum = numel(x_label);

%Specify prior distribution bounds
beta_lb = 0.2632*ones(4,1); beta_ub = 3*0.2632*ones(4,1);
AscertainProb_lb = 0*ones(6,1); AscertainProb_ub = 0.05*ones(6,1);
if MultiSeasonImmFlag == 0
    ExpHist_lb = [0;0;0]; ExpHist_ub = [1;1;1];
else
    ExpHist_lb = [0;0;0;0]; ExpHist_ub = [1;1;1;1];
end

%Concatenate bounds for all parameters
AllParam_lb = [beta_lb;ExpHist_lb;AscertainProb_lb];
AllParam_ub = [beta_ub;ExpHist_ub;AscertainProb_ub];

%Sepcify number of prior pdf points to plot
PriorPdfSamplePts = 100; %Number of points to compute value of prior pdf at
PriorPdfSampleVals = zeros(ParamNum,PriorPdfSamplePts); %Storage array for prior val sample points
PriorPdfVals = zeros(ParamNum,PriorPdfSamplePts); %Storage array for prior values

%Evaluate prior pdf at requested number of points
for ii = 1:ParamNum
    PriorPdfSampleVals(ii,:) = linspace(AllParam_lb(ii),AllParam_ub(ii),PriorPdfSamplePts);
    PriorPdfVals(ii,:) = unifpdf(PriorPdfSampleVals(ii,:),AllParam_lb(ii),AllParam_ub(ii));
end


%--------------------------------------------------------------------------
%%% SPECIFY PLOTTING VARIABLES
%--------------------------------------------------------------------------
SubPlotRowNum = 3;

if SeasonsFitted == 4 || SeasonsFitted == 5
    SubPlotColNum = 4;
else %Equivalent to SeasonsFitted == 6
    SubPlotColNum = 5;
end

FigSizeWidthScale = 3.5;
FigSizeHeightScale = 2.5;

%Specify subpanels that will be used (those omitted will remain empty)
SubPanelFitNumber_FourSeasons = 1:12;
SubPanelFitNumber = SubPanelFitNumber_FourSeasons;

if SeasonsFitted == 4
    SubPanelFitNumber = 1:11;
elseif SeasonsFitted == 5
    SubPanelFitNumber = 1:12;
else
    %Choose file based on whether original or extended model is
    %specified
    if MultiSeasonImmFlag == 0
        SubPanelFitNumber = 1:13;
    else
        SubPanelFitNumber = [1 2 3 4 6:15];
    end
end

%--------------------------------------------------------------------------
%%% DECLARE HISTOGRAM NORMALISATION OPTION
%--------------------------------------------------------------------------
% NormOption = 1; %1 - pdf; otherwise, plots relative probability
% PriorDataFlag = 1; %Flag variable: 0 (1) -> Do not (Do) include prior distribution on plots
% PriorDataVar = {PriorDataFlag,PriorPdfSampleVals,PriorPdfVals};

% %Plot histograms with priors
% PostDistHistograms(PosteriorParams,x_label,SynthDataVar,...
%                                 SubPlotRowNum,SubPlotColNum,...
%                                 FigSizeWidthScale,FigSizeHeightScale,...
%                                 NormOption)

%Plot histograms omitting priors
NormOption = 0; %1 - pdf; otherwise, plots relative probability
PriorDataFlag = 0;
PriorDataVar = {PriorDataFlag,PriorPdfSampleVals,PriorPdfVals};
PostDistHistograms(PosteriorParams,SeasonsFitted,MultiSeasonImmFlag,x_label,SynthDataVar,PriorDataVar,...
                                SubPlotRowNum,SubPlotColNum,SubPanelFitNumber,...
                                FigSizeWidthScale,FigSizeHeightScale,...
                                NormOption)
end

function PostDistHistograms(PosteriorParams,SeasonsFitted,MultiSeasonImmFlag,...
                                        x_label,SynthDataVar,PriorDataVar,...
                                        SubPlotRowNum,SubPlotColNum,SubPanelFitNumber,...
                                            FigSizeWidthScale,FigSizeHeightScale,...
                                            NormOption)
%Inputs:
%   PosteriorParams - (array, float) Particle sets generated by inference procedure
%   SeasonsFitted - (scalar, integer) Number of influenza seasons worth of data fit to in the inference scheme 
%   MultiSeasonImmFlag - (indicator variable) Model type used. 
%                                             0: Original, immunity propagation lasts a single season; 
%                                             1: Extended, allowing immunity to be propagated for more than a single season
%   x_label - (string) Parameter labels
%   SynthDataVar - (cell) Two tuple cell.
%       -> First entry: Flag variable. 1 if true parameter values know. 0 otherwise
%       -> Second entry: Vector of true parameter values
%   PriorDataVar - (cell) Two tuple cell.
%       -> First entry: Flag variable. 1 if want to plot prior distribution. 0 otherwise
%       -> Second entry: Vector of prior distribution values (x & f(x))
%   SubPlotNumRow, SubPlotColNum - (Integers) Specify partitioning of subplot grid
%   SubPanelFitNumber - (Integers) Specify panels to be used in figure plot
%   FigSizeWidthScale,FigSizeHeightScale - (floats) Specify scaling of figure size 
%   SaveFileName - (string) Specify location to save plots
%   NormOption - (scalar) histogram normalisation option

%Unpack SynthDataVar
SynthDataFlag = SynthDataVar{1};
ParamTrueVal = SynthDataVar{2};

%Unpack PriorDataVar 
PriorDataFlag = PriorDataVar{1};
PriorPdfSampleVals = PriorDataVar{2};
PriorPdfVals = PriorDataVar{3};

%Set up figure
figure('Color',[1 1 1]);
clf
position = [10, 10, FigSizeWidthScale*550, FigSizeHeightScale*450];
set(0, 'DefaultFigurePosition', position);

%Obtain number of Params from input data
ParamNum = numel(x_label);

for ii = 1:ParamNum
    subplot(SubPlotRowNum,SubPlotColNum,SubPanelFitNumber(ii)) 
    hold on
    
    if NormOption == 1
        %Construct histogram, plot of pdf
        h = histogram(PosteriorParams(:,ii),50,'Normalization','pdf','DisplayName','Posterior');
    else
        %Construct histogram, plot of relative probability
        h = histogram(PosteriorParams(:,ii),20,'Normalization','probability','DisplayName','Posterior');
    end
    
    %Get maximum bar height of histogram
    MaxBarHeight = max(h.Values);
    MaxYVal = MaxBarHeight*1.05;
    
    %Plot modal value
    [~,ModalIdx] = max(h.Values);
    ModalVal = mean(h.BinEdges(ModalIdx:ModalIdx+1));
    %p1 = plot([ModalVal ModalVal],[0 MaxYVal],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth', 1.5,'DisplayName','Mode');   
    
    %Plot median line
    MedianVal = median(PosteriorParams(:,ii));
    p2 = plot([MedianVal MedianVal],[0 MaxYVal],'Color',[0.8 0 0],'LineStyle','-','LineWidth', 1.5,'DisplayName','Median');
    
    if SynthDataFlag == 1
        %Plot "true" value line
        p3 = plot([ParamTrueVal(ii) ParamTrueVal(ii)],[0 MaxYVal],'Color',[0 0 0],'LineWidth',1.5,'DisplayName','True');
    end
    
    if PriorDataFlag == 1 && NormOption == 1
        %Plot prior distribution pdf   
        p4 = plot(PriorPdfSampleVals(ii,:),PriorPdfVals(ii,:),'Color',[0 0.8 0],'LineStyle','-.','LineWidth',1.5,'DisplayName','Prior');
    end
    
    %Set subplot properties
    set(gca,'Fontsize',16,'LineWidth',1);
    
    %Set axes labels
    xlabel(x_label(ii))
    if mod(SubPanelFitNumber(ii),SubPlotColNum) == 1
        if NormOption == 1
            ylabel('pdf');
        else
            ylabel('Proportion');
        end
    end
    
    %Set y-axis limits
    ylim([0 MaxYVal])
    
%     %Set x-axis limits for transmissibility params
%     if ii < 5
%        xlim([0.379 0.441])
%     elseif ii >7
%         xlim([0.9e-3 4.1e-3])
%     end
    
    
    %Set x-axis limits (if want to plot prior distribution)
    if PriorDataFlag == 1 && NormOption == 1
        xlim([min(PriorPdfSampleVals(ii,:)) max(PriorPdfSampleVals(ii,:))]);
    end
    box on
    
end


%Add legend
if SynthDataFlag == 1 && PriorDataFlag == 1 && NormOption == 1
    %Add "true" value line and prior distribution
    %leg1 = legend([h,p1,p2,p3,p4]);
    leg1 = legend([h,p2,p3,p4],'Fontsize',16);
elseif SynthDataFlag == 1 %"True" value, no prior
    %leg1 = legend([h,p1,p2,p3]);
    leg1 = legend([h,p2,p3],'Fontsize',16);
elseif PriorDataFlag == 1 && NormOption == 1 %Add prior distribution
    %leg1 = legend([h,p1,p2,p4]);
    leg1 = legend([h,p2,p4],'Fontsize',16);
else %No "true" value or prior
    %leg1 = legend([h,p1,p2]);
    leg1 = legend([h,p2],'Fontsize',16);
end

%Position legend based on model type and number of seasons fitted to
if SeasonsFitted == 4
    set(leg1,...
        'Position',[0.75 0.205 0.07 0.05],...
        'FontSize',16);
elseif SeasonsFitted == 5
    set(leg1,...
    'Position',[0.853 0.862 0.07 0.05],...
        'FontSize',16);
else %SeasonsFitted == 6 
    if MultiSeasonImmFlag == 0
        set(leg1,...
            'Position',[0.622 0.205 0.07 0.05],...
            'FontSize',16);
    else
        set(leg1,...
            'Position',[0.78 0.80 0.07 0.05],...
            'FontSize',16);
    end
end


end                            
                            