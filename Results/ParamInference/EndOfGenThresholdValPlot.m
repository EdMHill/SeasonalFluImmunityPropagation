%Purpose:
%Plot ABC inference scheme end-of-generation threshold profile

%Author: Ed Hill
%--------------------------------------------------------------------------


%% Get data
clear variables

%--------------------------------------------------------------------------
%%% IMPORT END-OF-GENERATION THRESHOLD VALUE DATA
%--------------------------------------------------------------------------
FileName = 'SixSeasonFitData/ABCEndOfGenThresholdVals_SixSeasonFit.txt';
EndOfGenThresholdVals = dlmread(FileName);

%--------------------------------------------------------------------------
%%% DECLARE MAXIMUM GENERATION VALUE TO PLOT UP TO
%--------------------------------------------------------------------------
TotalGens = numel(EndOfGenThresholdVals);

%--------------------------------------------------------------------------
%%% DECLARE PORTION OF GENERATIONS TO BE PLOTTED IN INSET PANEL
%--------------------------------------------------------------------------
InsetPanelPortion = 0.25;

%--------------------------------------------------------------------------
%%% SPECIFY PLOTTING VARIABLES
%--------------------------------------------------------------------------
FigSizeWidthScale = 2;
FigSizeHeightScale = 1.5;


%--------------------------------------------------------------------------
%%% CONSTRUCT MAIN PANEL FIGURE
%--------------------------------------------------------------------------
%Set up figure
figure('Color',[1 1 1]);
clf
position = [10, 50, FigSizeWidthScale*550, FigSizeHeightScale*450];
set(0, 'DefaultFigurePosition', position);
hold on

%Plot theshold value at conclusion of each generation
plot(1:1:TotalGens,EndOfGenThresholdVals,'Linewidth',1.5)

%Specify x-axis limits
xlim([0 TotalGens]);

%State axes labels
xlabel('Generations completed')
ylabel({'End-of-generation';'summary staistic threshold'})

%Set plot properties
set(gca,'Fontsize',16,'LineWidth',1);

%--------------------------------------------------------------------------
%%% CONSTRUCT INSET PANEL FIGURE
%--------------------------------------------------------------------------

% create smaller axes in top right, and plot on it
axes('Position',[.5 .5 .4 .4])
box on

%Get entry corresponding to beginning of last quarter of generations
InsetStartGen = (1 - InsetPanelPortion) * TotalGens;

%Plot theshold value at conclusion of each generation
plot(InsetStartGen:1:TotalGens,EndOfGenThresholdVals(InsetStartGen:end),'Linewidth',1.5)

%Specify x-axis limits
xlim([InsetStartGen TotalGens]);

%Set plot properties
set(gca,'Fontsize',16,'LineWidth',1);