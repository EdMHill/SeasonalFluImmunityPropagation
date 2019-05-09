%Purpose:
%Plot temporal profiles (of influenza infection prevalence)
%from en masse model simulations

%Author: Ed Hill
%--------------------------------------------------------------------------

clear variables

%% Get data from model simulation

%Load required simulation data
InputDataFileName = 'ModelSimnData/ModelSimns_SixSeasonFit.mat';
load(InputDataFileName,'SimnData')

%Load cells containing infected timeseries (second column of SimnData)
InfTimeProfileAllRuns = SimnData(:,2);

%Set SimnRunType (1 -> historical; 2 -> forward simn)
SimnRunType = 1;

%Specify number of seasons to be plotted (starting from 2009/2010 influenza season)
TotalYrs = 9;

%% Amass runs into single array

%Initialise storage array
%  -> Row per run
%  -> column per timestep
RunNum = numel(InfTimeProfileAllRuns);
LengthT = numel(InfTimeProfileAllRuns{1});
InfTimeProfileAllRunsArray = zeros(RunNum,LengthT);

%Iterate through each individual run, add timeseries info to InfTimeProfileAllRunsArray
for ii = 1:RunNum
    InfTimeProfileAllRunsArray(ii,:) = InfTimeProfileAllRuns{ii};
end

%% Plot all timeseries on single plot
figure('Color',[1 1 1]);
clf
position = [100, 100, 2*550, 450];
set(0, 'DefaultFigurePosition', position);
hold on

%Plot temporal profiles for each simulation replicate
plot(1:LengthT,InfTimeProfileAllRunsArray(100:100:end,:),'Color',[0.5,0.5,0.5],'LineWidth',0.5);

if SimnRunType == 1  %Historical seasons only
        
    %Specify x ticks
    xticks(0:365:TotalYrs*365)
    
    %Construct x tick labels
    XVals = cell(TotalYrs+1,1);
    XVals{1} = '0';
    for ii = 1:TotalYrs
        XVals{ii+1} = num2str(ii);
        xticklabels(XVals)
    end
    
    %Add lines denoting mid-season mark
    MidSeasonMark = (0:366:(TotalYrs-1)*366) + (366/2);
    for ii = 1:TotalYrs
        plot([MidSeasonMark(ii) MidSeasonMark(ii)],[0 1],'--','Color',[0.8,0,0],'LineWidth',1)
    end
    
    %Specify x tick labels
    %xticks([0 1*365 2*365 3*365 4*365 5*365 6*365 7*365 8*365])   
    xticks(0:365:(TotalYrs-1)*365)
    xticklabels({'0','1','2','3','4','5','6','7','8'})
    
    %Add lines denoting mid-season mark
    MidSeasonMark = (0:366:(TotalYrs-2)*366) + (366/2);
    for ii = 1:8
        plot([MidSeasonMark(ii) MidSeasonMark(ii)],[0 1],'--','Color',[0.8,0,0],'LineWidth',1)
    end
    
end

%Specify axes labels
xlabel('Elapsed time (years)')
ylabel('Proportion infected')

%Specify ylim
ylim([0 max(InfTimeProfileAllRunsArray(:))+0.01])

%Specify plot properties
set(gca,'FontSize',16)
set(gca,'LineWidth',1)

box on