%Purpose:
%Type vs Type, subtype vs subtype, lineage vs lineage scatter plots
%Depict simulated data and highlight observed data

%Author: Ed Hill
%--------------------------------------------------------------------------

clear variables

%% Specify data that was fit to
SynthDataFlag = 1; %Indicator variable. 0: empirical data; 1 - synthetic data

%--------------------------------------------------------------------------
%%% Read in observed data
%--------------------------------------------------------------------------

%Declare number of seasons to be plotted
SeasonsToPlot = 6; %From 2012/13 influenza season onward
                    %value of 6 covers 2012/13-2017/18

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
        ModelSimnFileName = 'ModelSimnData/ModelSimns_SixSeasonFit.mat';
    else
        error('Misspecified SeasonsToPlot value of %f. SeasonsToPlot must take a value of 4,5 or 6',SeasonsToPlot)
    end
elseif SynthDataFlag == 1 %use outputs from fitting to synthetic data
    DataFileName = '../../Data/ILIData/SynthData_InfluenzaPositiveGPConsultRateByStrain_2009to2018.csv';
    ObvsDataTemp = dlmread(DataFileName);
    ObvsData = ObvsDataTemp(4:end,:);
    
    ModelSimnFileName = 'ModelSimnData/ModelSimns_SynthDataFit.mat';
    
    SeasonsToPlot = 6; %By deafult, revert SeasonsToPlot to 6    
else
    error('Misspecified SynthDataFlag value, set as %f. SynthDataFlag must take value 0 or 1.',SynthDataFlag);
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
Bootstrap_50PercentReps = BootstrapReps*0.5; %Number of replicates used for formulating regions of points around point estimate
Bootstrap_95PercentReps = BootstrapReps*0.95;

BootstrapDataAmalg = zeros(BootstrapReps,4,SeasonsToPlot);
for jj = 1:BootstrapReps
    for kk = 1:SeasonsToPlot
        % Access BootstrapRepsData from entry 4 onward (corresponding to 2012/13 season and beyond)
        BootstrapDataAmalg(jj,:,kk) = BootstrapRepsData{kk+3}(jj,:);
    end
end

%Calculate rate for each flu type, per bootstrap sample
AllInfA_Bootstrap = sum(BootstrapDataAmalg(:,1:2,:),2); %Sum A/H1 and A/H3 columns
AllInfB_Bootstrap = sum(BootstrapDataAmalg(:,3:4,:),2); %Sum B/Yam and B/Vic columns

%Get maximum and minimum for Inf A and Inf B across all seasons
MaxInfA_Bootstrap = max(max(AllInfA_Bootstrap));
MinInfA_Bootstrap = min(min(AllInfA_Bootstrap));

MaxInfB_Bootstrap = max(max(AllInfB_Bootstrap));
MinInfB_Bootstrap = min(min(AllInfB_Bootstrap));

%Get maximum and minimum for each strain type/lineage
MaxByStrain_Bootstrap = max(max(BootstrapDataAmalg),[],3);
MinByStrain_Bootstrap = min(min(BootstrapDataAmalg),[],3);


%--------------------------------------------------------------------------
%%% Read in simulated epidemiological data
%--------------------------------------------------------------------------
load(ModelSimnFileName,'SimnData') %Load SimnData variable
RunNum = size(SimnData,1); %Number of replicates equals rows of SimnData

%--------------------------------------------------------------------------
%%% Get simulated data into a single, 3D array
%%% Acces SimnData cell, which has a cell per simulation run
%--------------------------------------------------------------------------
%Column per strain
%  -> Col 1: A/H1; Col 2: A/H3; Col 3: B/Yam; Col 4: B/Vic
%3rd dimension slice per season
% -> Slice 1: 2012/2013 season; Slice 2: 2013/2014 season etc
SimnDataAmalg = zeros(RunNum,4,SeasonsToPlot);
for jj = 1:RunNum
    for kk = 1:SeasonsToPlot
        SimnDataAmalg(jj,:,kk) = SimnData{jj,1}(kk,:);
    end
end

%Calculate axes limits per panel to be used across all figures
AllInfASimn = sum(SimnDataAmalg(:,1:2,:),2); %Sum A/H1 and A/H3 columns
AllInfBSimn = sum(SimnDataAmalg(:,3:4,:),2); %Sum B/Yam and B/Vic columns

%Get maximum and minimum for Inf A and Inf B across all seasons
MaxInfA_Simn = max(max(AllInfASimn));
MinInfA_Simn = min(min(AllInfASimn));

MaxInfB_Simn = max(max(AllInfBSimn));
MinInfB_Simn = min(min(AllInfBSimn));

%Get maximum and minimum for each strain type/lineage
MaxByStrain_Simn = max(max(SimnDataAmalg),[],3);
MinByStrain_Simn = min(min(SimnDataAmalg),[],3);


%--------------------------------------------------------------------------
%%% Get axes limits, based on max/min values of simualted & bootstrap data
%--------------------------------------------------------------------------
MaxInfA = max(MaxInfA_Simn,MaxInfA_Bootstrap);
MinInfA = min(MinInfA_Simn,MinInfA_Bootstrap);

MaxInfB = max(MaxInfB_Simn,MaxInfB_Bootstrap);
MinInfB = min(MinInfB_Simn,MinInfB_Bootstrap);

MaxByStrain = max(MaxByStrain_Simn,MaxByStrain_Bootstrap);
MinByStrain = min(MinByStrain_Simn,MinByStrain_Bootstrap);

%--------------------------------------------------------------------------
%%% Iterate through each season and construct plots with 4 panels. 
%%% Left half for Type A vs Type B
%%% Top right for A/H1 vs A/H3
%%% Bottom right for B/Yam vs B/Vic
%--------------------------------------------------------------------------
TitleVec = {'2012/13','2013/14','2014/15','2015/16','2016/17','2017/18'};
PlotLBscaling = 0.98; %Scaling to be applied to minimum values across all seasons
PlotUBscaling = 1.02; %Scaling to be applied to maximum values across all seasons
for ii = 1:SeasonsToPlot
    %Horizontal bar plot for each row, overlay bars for each subytpe/lineage
    figure('Color',[1 1 1]);
    clf
    position = [100, 100, 3*550, 2*450];
    set(0, 'DefaultFigurePosition', position);
    hold on
    
    %Scatter plots, with weights mapped to point size sz
    %function format: scatter(x,y,sz)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   All A vs All B  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,3,[1 2 4 5])
    hold on
    
    %If using empirical data, plot bootstrap sample regions
    if SynthDataFlag == 0
        %----------------------------------------------------------------------
        %%% For bootstrap replicates, find points within 50% and 95% regions of
        %%% observed data point
        %----------------------------------------------------------------------
        X = [AllInfA_Bootstrap(:,:,ii),AllInfB_Bootstrap(:,:,ii)];
        Y = [ObvsData(ii,1)+ObvsData(ii,2),ObvsData(ii,3)+ObvsData(ii,4)];
        Idx_50Percentile = knnsearch(X,Y,'K',Bootstrap_50PercentReps); %finds the K nearest neighbors in X for each query point in Y
        Idx_95Percentile = knnsearch(X,Y,'K',Bootstrap_95PercentReps); %finds the K nearest neighbors in X for each query point in Y
        
        %Assign samples within regions of interest to variable
        Samples_50Percentile = X(Idx_50Percentile,:);
        Samples_95Percentile = X(Idx_95Percentile,:);
        
        %Get index of points on boundary
        Boundary_50Percentile = boundary(Samples_50Percentile(:,1),Samples_50Percentile(:,2));
        Boundary_95Percentile = boundary(Samples_95Percentile(:,1),Samples_95Percentile(:,2));
        %----------------------------------------------------------------------
        
        %Plot bootstrap replicates
        f1 = fill(Samples_50Percentile(Boundary_50Percentile,1),Samples_50Percentile(Boundary_50Percentile,2),[0 0.8 0],'DisplayName','50% CI (Obvs data)');
        f2 = fill(Samples_95Percentile(Boundary_95Percentile,1),Samples_95Percentile(Boundary_95Percentile,2),[0 0.8 0],'DisplayName','95% CI (Obvs data)');

        % For filled regions, choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
        set(f1,'facealpha',1.0)
        set(f2,'facealpha',.4)
    end
    
    %Plot empirical data point estimate
    s1 = scatter(ObvsData(ii,1)+ObvsData(ii,2),ObvsData(ii,3)+ObvsData(ii,4),50,'d','MarkerEdgeColor',[0.8 0 0],...
        'MarkerFaceColor',[0.8 0 0],'DisplayName','Observed');

    %Plot simulated data
    s2 = scatter(AllInfASimn(:,:,ii),AllInfBSimn(:,:,ii),[],[0 0.4470 0.7410],'filled','DisplayName','Predicted');
     
    %Set axes labels
    xlabel('Type A')
    ylabel('Type B')
    
    %Set axes limits
    %xlim([max(0,MinInfA*PlotLBscaling) MaxInfA*PlotUBscaling])
    %ylim([max(0,MinInfB*PlotLBscaling) MaxInfB*PlotUBscaling])
    xlim([0 max(MaxInfA,ObvsData(ii,1)+ObvsData(ii,2))*PlotUBscaling])
    ylim([0 max(MaxInfB,ObvsData(ii,3)+ObvsData(ii,4))*PlotUBscaling])

    %Add legend. 
    if SynthDataFlag == 0 % If using empirical data label bootstrap sample regions
        leg1 = legend([s1,f1,f2,s2]);
    elseif SynthDataFlag == 1
        leg1 = legend([s1,s2]);
    end
    
    %Set subplot properties
    set(gca,'Fontsize',14,'LineWidth',1);
    box on
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% A(H1N1)pdm09 vs A(H3N2) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,3,3)
    hold on
        
    %If using empirical data, plot bootstrap sample regions
    if SynthDataFlag == 0
        %----------------------------------------------------------------------
        %%% For bootstrap replicates, find points within 50% and 95% regions of
        %%% observed data point
        %----------------------------------------------------------------------
        X = [BootstrapDataAmalg(:,1,ii),BootstrapDataAmalg(:,2,ii)];
        Y = [ObvsData(ii,1),ObvsData(ii,2)];
        Idx_50Percentile = knnsearch(X,Y,'K',Bootstrap_50PercentReps); %finds the K nearest neighbors in X for each query point in Y
        Idx_95Percentile = knnsearch(X,Y,'K',Bootstrap_95PercentReps); %finds the K nearest neighbors in X for each query point in Y
        
        %Assign samples within regions of interest to variable
        Samples_50Percentile = X(Idx_50Percentile,:);
        Samples_95Percentile = X(Idx_95Percentile,:);
        
        %Get index of points on boundary
        Boundary_50Percentile = boundary(Samples_50Percentile(:,1),Samples_50Percentile(:,2));
        Boundary_95Percentile = boundary(Samples_95Percentile(:,1),Samples_95Percentile(:,2));
        %----------------------------------------------------------------------
        
        %Plot bootstrap replicates
        f3 = fill(Samples_50Percentile(Boundary_50Percentile,1),Samples_50Percentile(Boundary_50Percentile,2),[0 0.8 0],'DisplayName','50% CI (Obvs data)');
        f4 = fill(Samples_95Percentile(Boundary_95Percentile,1),Samples_95Percentile(Boundary_95Percentile,2),[0 0.8 0],'DisplayName','95% CI (Obvs data)');
    
        % For filled regions, choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
        set(f3,'facealpha',1.0)
        set(f4,'facealpha',.4)
    end

    %Plot empirical data point estimate
    scatter(ObvsData(ii,1),ObvsData(ii,2),75,'d','MarkerEdgeColor',[0.8 0 0],...
        'MarkerFaceColor',[0.8 0 0],'DisplayName','Observed')
    
    %Plot simulated data
    scatter(SimnDataAmalg(:,1,ii),SimnDataAmalg(:,2,ii),[],[0 0.4470 0.7410],'filled','DisplayName','Predicted')

          
    %Set axes labels
    xlabel('A(H1N1)pdm09')
    ylabel('A(H3N2)')     
    
    %Set axes limits
    %xlim([max(0,MinByStrain(1)*PlotLBscaling) MaxByStrain(1)*PlotUBscaling])
    %ylim([max(0,MinByStrain(2)*PlotLBscaling) MaxByStrain(2)*PlotUBscaling])
    xlim([0 max(MaxByStrain(1),ObvsData(ii,1))*PlotUBscaling])
    ylim([0 max(MaxByStrain(2),ObvsData(ii,2))*PlotUBscaling])

    %Set subplot properties
    set(gca,'Fontsize',14,'LineWidth',1);
    box on
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% B/Yam vs B/Vic %%%
    %%%%%%%%%%%%%%%%%%%%%%   
    subplot(2,3,6)
    hold on
        
    %If using empirical data, plot bootstrap sample regions
    if SynthDataFlag == 0
        %----------------------------------------------------------------------
        %%% For bootstrap replicates, find points within 50% and 95% regions of
        %%% observed data point
        %----------------------------------------------------------------------
        X = [BootstrapDataAmalg(:,3,ii),BootstrapDataAmalg(:,4,ii)];
        Y = [ObvsData(ii,3),ObvsData(ii,4)];
        Idx_50Percentile = knnsearch(X,Y,'K',Bootstrap_50PercentReps); %finds the K nearest neighbors in X for each query point in Y
        Idx_95Percentile = knnsearch(X,Y,'K',Bootstrap_95PercentReps); %finds the K nearest neighbors in X for each query point in Y
        
        %Assign samples within regions of interest to variable
        Samples_50Percentile = X(Idx_50Percentile,:);
        Samples_95Percentile = X(Idx_95Percentile,:);
        
        %Get index of points on boundary
        Boundary_50Percentile = boundary(Samples_50Percentile(:,1),Samples_50Percentile(:,2));
        Boundary_95Percentile = boundary(Samples_95Percentile(:,1),Samples_95Percentile(:,2));
        %----------------------------------------------------------------------
        
        %Plot bootstrap replicates
        f5 = fill(Samples_50Percentile(Boundary_50Percentile,1),Samples_50Percentile(Boundary_50Percentile,2),[0 0.8 0],'DisplayName','50% CI (Obvs data)');
        f6 = fill(Samples_95Percentile(Boundary_95Percentile,1),Samples_95Percentile(Boundary_95Percentile,2),[0 0.8 0],'DisplayName','95% CI (Obvs data)');
        
        % For filled regions, choose a number between 0 (invisible) and 1 (opaque) for facealpha.
        set(f5,'facealpha',1.0)
        set(f6,'facealpha',.4)
    end
    %Plot empirical data point estimate
    scatter(ObvsData(ii,3),ObvsData(ii,4),75,'d','MarkerEdgeColor',[0.8 0 0],...
        'MarkerFaceColor',[0.8 0 0],'DisplayName','Observed')

    %Plot simulated data
    scatter(SimnDataAmalg(:,3,ii),SimnDataAmalg(:,4,ii),[],[0 0.4470 0.7410],'filled')

    %Set axes labels
    xlabel('B/Yamagata')
    ylabel('B/Victoria')
    
    %Set axes limits
    %xlim([max(0,MinByStrain(3)*PlotLBscaling) MaxByStrain(3)*PlotUBscaling])
    %ylim([max(0,MinByStrain(4)*PlotLBscaling) MaxByStrain(4)*PlotUBscaling])
    xlim([0 max(MaxByStrain(3),ObvsData(ii,3))*PlotUBscaling])
    ylim([0 max(MaxByStrain(4),ObvsData(ii,4))*PlotUBscaling])
    
    %Set subplot properties
    set(gca,'Fontsize',14,'LineWidth',1);
    box on
    
    %Set figure title
    FullTitle = [TitleVec{ii},' influenza season'];
    sgtitle(FullTitle,'Fontsize',16,'FontWeight','bold')   
    
end

