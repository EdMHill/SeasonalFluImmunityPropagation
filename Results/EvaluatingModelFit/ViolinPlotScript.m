%Purpose:
%Violin plots of flu positive GP consultations 
%Includes bootstrap data

%Produced from files containing multiple run outputs
%Per run, GP consultations for influeza per strain (per 100,000 population)

%Author: Ed Hill
%--------------------------------------------------------------------------

clear variables

%--------------------------------------------------------------------------
% ADD PATH DEPENDENCIES
%--------------------------------------------------------------------------
addpath('../MatlabFileExchangeScripts') 
addpath('../MatlabFileExchangeScripts/Violinplot-Matlab') 


%--------------------------------------------------------------------------
% DEFINE FUNCTIONS INPUT VALUES
%--------------------------------------------------------------------------
SynthDataFlag = 0; %Indicator variable. 0: empirical data; 1: synthetic data


MultiSeasonImmFlag = 0;  %Model type used. 
                         %0: Original, immunity propagation lasts a single season; 
                         %1:Extended, allowing immunity to be propagated for more than a single season

                         
%--------------------------------------------------------------------------
%%% DECLARE NUMBER OF SEASONS TO PLOT
%%% DATA FITTED FROM 2012/13 INFLUENZA SEASON ONWARDS
%%% TAKES VALUE OF 4 (2012/13 - 2015/16),5 (2012/13 - 2016/17) OR 6 (2012/13 - 2017/18)
%--------------------------------------------------------------------------
SeasonsToPlot = 6;
       
%If want to plot outputs fitting to data with the extended model or syntehtic data, 
%only six seaon fit is available
%By deafult, revert SeasonsToPlot to 6
if SynthDataFlag == 1|| MultiSeasonImmFlag == 1
    SeasonsToPlot = 6;
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

if SeasonsToPlot ~=4 && SeasonsToPlot ~=5 && SeasonsToPlot ~=6
    error('SeasonsToPlot must take a value of 4,5 or 6')
end

%--------------------------------------------------------------------------
% CALL FUNCTION
%--------------------------------------------------------------------------
ConstructViolinPlot(SynthDataFlag,MultiSeasonImmFlag,SeasonsToPlot)

function ConstructViolinPlot(SynthDataFlag,MultiSeasonImmFlag,SeasonsToPlot)
%   SynthDataFlag - (indicator variable) 0: empirical data; 1: synthetic data
%   MultiSeasonImmFlag - (indicator variable) Model type used. 
%                                             0: Original, immunity propagation lasts a single season; 
%                                             1: Extended, allowing immunity to be propagated for more than a single season
%   SeasonsToPlot - (scalar, integer) Number of influenza seasons worth of data fit to in the inference scheme 


%--------------------------------------------------------------------------
%% GET DATA THAT WAS FITTED TO BY INFERENCE PROCEDURE
%--------------------------------------------------------------------------
if SynthDataFlag == 0 %use outputs from fitting to empirical data
    DataFileName = '../../Data/ILIData/EmpData_InfluenzaPositiveGPConsultRateByStrain_2009to2018.csv';
    ObvsDataTemp = dlmread(DataFileName);
    if SeasonsToPlot == 4
        ObvsData = ObvsDataTemp(4:end-2,:);
        ModelSimnFileName = 'ModelSimnData/ModelSimns_FourSeasonFit.mat';
    elseif SeasonsToPlot == 5
        ObvsData = ObvsDataTemp(4:end-1,:);
        ModelSimnFileName = 'ModelSimnData/ModelSimns_FiveSeasonFit.mat';
    elseif SeasonsToPlot == 6
        ObvsData = ObvsDataTemp(4:end,:);       
        if MultiSeasonImmFlag == 0
            ModelSimnFileName = 'ModelSimnData/ModelSimns_SixSeasonFit.mat';
        else
            ModelSimnFileName = 'ModelSimnData/ExtendedModelSimns_SixSeasonFit.mat';
        end
    else
        error('Misspecified SeasonsToPlot value of %f. SeasonsToPlot must take a value of 4,5 or 6',SeasonsToPlot)
    end
elseif SynthDataFlag == 1 %use outputs from fitting to synthetic data      
    if MultiSeasonImmFlag == 0
        DataFileName = '../../Data/ILIData/SynthData_InfluenzaPositiveGPConsultRateByStrain_2009to2018.csv';
        ModelSimnFileName = 'ModelSimnData/ModelSimns_SynthDataFit.mat';
    else
        DataFileName = '../../Data/ILIData/ExtendedModel_SynthData_InfluenzaPositiveGPConsultRateByStrain_2009to2018.csv';
        ModelSimnFileName = 'ModelSimnData/ExtendedModelSimns_SynthDataFit.mat';
    end
    
    %Load synthetic data
    ObvsDataTemp = dlmread(DataFileName);
    ObvsData = ObvsDataTemp(4:end,:);
end

%--------------------------------------------------------------------------
%%% Import bootstrap replicates
%--------------------------------------------------------------------------
BootStrapFileData = load('../../Data/ILIData/InfluenzaPositiveGPConsultRateByStrain_2009to2018_BootstrapSamples.mat');

BootstrapRepsData = BootStrapFileData.StrainCaseSeasonRate_BootstrapReps;

% Get bootstrap data into a single, 3D array
%Column per strain
%  -> Col 1: A/H1; Col 2: A/H3; Col 3: B/Yam; Col 4: B/Vic
%3rd dimension slice per season
% -> Slice 1: 2012/2013 season; Slice 2: 2013/2014 season etc
BootstrapReps = size(BootstrapRepsData{1},1);
BootstrapDataAmalg = zeros(BootstrapReps,4,SeasonsToPlot);
for jj = 1:BootstrapReps
    for kk = 1:SeasonsToPlot
        % Access BootstrapRepsData from entry 4 onward (corresponding to 2012/13 season and beyond)
        BootstrapDataAmalg(jj,:,kk) = BootstrapRepsData{kk+3}(jj,:);
    end
end


%--------------------------------------------------------------------------
%% GET DATA FROM M3 SIMULATION RUNS
%--------------------------------------------------------------------------
%Load required simulation data
load(ModelSimnFileName,'SimnData') %Load SimnData variable

StrainSeasonRatePerStrainAllRuns = SimnData(:,1);

%% For each subtype/lineage per season, compute desired percentiles of specified output statistic

%--------------------------------------------------------------------------
%%% Get simulated data into a single, 3D array
%%% Acces SimnData cell, which has a cell per simulation run
%--------------------------------------------------------------------------
%Column per strain
%  -> Col 1: A/H1; Col 2: A/H3; Col 3: B/Yam; Col 4: B/Vic
%3rd dimension slice per season
% -> Slice 1: 2012/2013 season; Slice 2: 2013/2014 season etc
RunNum = size(StrainSeasonRatePerStrainAllRuns,1); %Number of replicates equals rows of StrainSeasonRatePerStrainAllRuns
SimnDataAmalg = zeros(RunNum,4,SeasonsToPlot);
for jj = 1:RunNum
    for kk = 1:SeasonsToPlot
        SimnDataAmalgTemp = StrainSeasonRatePerStrainAllRuns{jj}(kk,:);
        
        %Modify data for non-fitted seasons (i.e. those being predicted/forward simulated)
        %Compute proportion per subtype/lineage (in place of GP visit counts)
        SimnDataAmalg(jj,:,kk) = SimnDataAmalgTemp;       
    end
end

%--------------------------------------------------------------------------
%% For each subtype/lineage per season, plot observed and simulated data
%--------------------------------------------------------------------------

%Initialise figure
fig = figure('Color',[1 1 1]);
clf
position = [100, 100, 2*550, 2*450];
set(0, 'DefaultFigurePosition', position);

%Specify colour to be used per strain, set of RGB vectors
ColourByStrain = [1 0.5 0;
    0 0.8 0;
    0 0 0.8;
    0.5 0.5 0.5];

ColourByStrain_Bootstrap = [0.5 0.5 0.5];

%Specify figure properties dependent upon number of seasons plotted
% - Years
% - x-axis properties 
    if SeasonsToPlot == 4
        Year = [2013 2014 2015 2016];
        FluSeasonLabels = {'12/13','13/14','14/15','15/16'};
    elseif SeasonsToPlot == 5
        Year = [2013 2014 2015 2016 2017];
        FluSeasonLabels = {'12/13','13/14','14/15','15/16','16/17'};
    elseif SeasonsToPlot == 6
        Year = [2013 2014 2015 2016 2017 2018];
        FluSeasonLabels = {'12/13','13/14','14/15','15/16','16/17','17/18'};
    else
        error('Y axes limit details not set for SeasonsFitted value of %f', SeasonsFitted);
    end


%Iterate through each strain type, plot data
StrainNum = 4;
TitleVec = {'A(H1N1)pdm09','A(H3N2)','B/Yamagata','B/Victoria'};
for ii = 1:StrainNum
    
    %plot violins, colour based on strain type, 
    %with white symbol denoting median, red symbol denoting observed data

    %Panel 1: A(H1N1)pdm09
    %Panel 2: A(H3N2)
    %Panel 3: B/Victoria
    %Panel 4: B/Yamagata
    if ii < 3 %Top panels are for the two A subtypes
        subplot(2,2,ii)
    elseif ii == 3 %Plot B/Yamgata as bottom right panel
        subplot(2,2,4)
    else %Final strain plotted is B/Victoria
        subplot(2,2,3)
    end
    hold on
    
    %Want column per season, Row per simulation replicate
    StrainSpecificData = squeeze(SimnDataAmalg(:,ii,:)); %Remove singleton dimension
    BootstrapStrainData = squeeze(BootstrapDataAmalg(:,ii,:)); %Remove singleton dimension
    
    %violins_SimnData = violinplot(StrainSpecificData, FluSeasonLabels,'ShowData',false);
    %violins_BootstrapData = violinplot(BootstrapStrainData, FluSeasonLabels,'ShowData',false);

    %Get 95% interval for observed data, from bootstrap replicates
    BootstrapStrainData_CI = prctile(BootstrapStrainData,[2.5 97.5],1);
    
    %Plot my violin plot form, 
    %Third input (XLoc) specifying the position of each violin entity along the x-axis
    %Fourth input (XPosLabels) specifying the position of the axes labels
    XlocSpacing = 1;
    Xloc = 1:XlocSpacing:(XlocSpacing*5)+1;
    XPosLabels = 1:XlocSpacing:(XlocSpacing*5)+1;
    violins_SimnData = violinplotEH(StrainSpecificData, FluSeasonLabels,Xloc,XPosLabels,'ShowData',false);
   
    
    %Plot simulation data
    jj = 1;
    while jj <= SeasonsToPlot
        %Set violin plot region properties
        violins_SimnData(jj).ViolinColor = ColourByStrain(ii,:);
        violins_SimnData(jj).ViolinAlpha = 1.0; %shading transparency
        violins_SimnData(jj).EdgeColor = [0 0 0];
        
        %Set IQR properties (is a filled polygon! To hide, match ViolinColor)
        violins_SimnData(jj).BoxColor = ColourByStrain(ii,:); 

        %Set whisker line properties
        violins_SimnData(jj).WhiskerPlot.Color = [0 0 0];       
        
        %Set median marker properties
        violins_SimnData(jj).MedianPlot.MarkerFaceColor = [1 1 1];
        violins_SimnData(jj).MedianPlot.MarkerEdgeColor = [0 0 0];
        violins_SimnData(jj).MedianPlot.LineWidth = 1;
                    
        jj = jj + 1;
    end
    
    %Plot bootstrap data
    if SynthDataFlag == 0
        jj = 1;
        while jj <= SeasonsToPlot
            plot([jj jj],BootstrapStrainData_CI(:,jj),'r','LineWidth',2.5)
            
            jj = jj + 1;
        end
    end
    
    %Plot empirical data
    p1 = plot(1:1:SeasonsToPlot,ObvsData(:,ii),'x','Color','r','MarkerSize',15,'DisplayName','Observed data');
    
    if ii == 1 || ii == 4
        %Specify y-axis properties
        ylabel({'Influenza positive GP visits'; '(per 100,000)'})
    end

        
    %Set title
    title(TitleVec{ii})
    
    %Specify x-axis properties
    if ii > 2      
        xlabel('Influenza season');
    end
    
    %Set y-axis limits
    ylim([0 140])
    
    %Set x-axis limits
    if SeasonsToPlot == 4
        xlim([0 4.9])
    elseif SeasonsToPlot == 5
        xlim([0 5.9])
    elseif SeasonsToPlot == 6
        xlim([0 7])
    end
    
    %Add legend to first panel
    p2 = plot(nan,nan,'Color','r','DisplayName','95% bootstrapped CI','LineWidth',2);
    p3 = scatter(nan,nan,'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
        'DisplayName','Simulated median');
    p4 = fill(nan,nan,[0.5 0.5 0.5],'DisplayName','Predicted distribution','LineWidth',1);
    if ii == 1
        if SynthDataFlag == 0 %use outputs from fitting to empirical data
            leg1 = legend([p1,p2,p3,p4],'Location','NorthWest');
        elseif SynthDataFlag == 1 %use outputs from fitting to synthetic data
            leg1 = legend([p1,p3,p4],'Location','NorthWest');
        end
    end
    
    %Specify general axis properties
    set(gca,'FontSize',15)
    set(gca,'LineWidth',1.5)
    
    %Add surrounding box
    box on
    
    hold off
end

% %%% Add labels to panel %%
% annotation(fig,'textbox',...
%     [0.08 0.935 0.04 0.0377777772810724],...
%     'String',{'(a)'},...
%     'LineStyle','none','FontSize',15,'FontWeight','bold');
% 
% 
% annotation(fig,'textbox',...
%     [0.52 0.935 0.04 0.0377777772810724],...
%     'String',{'(b)'},...
%     'LineStyle','none','FontSize',15,'FontWeight','bold');
% 
% annotation(fig,'textbox',...
%     [0.08 0.46 0.04 0.0377777772810724],...
%     'String',{'(c)'},...
%     'LineStyle','none','FontSize',15,'FontWeight','bold');
% 
% annotation(fig,'textbox',...
%     [0.52 0.46 0.04 0.0377777772810724],...
%     'String',{'(d)'},...
%     'LineStyle','none','FontSize',15,'FontWeight','bold');
%--------------------------------------------------------------------------
end