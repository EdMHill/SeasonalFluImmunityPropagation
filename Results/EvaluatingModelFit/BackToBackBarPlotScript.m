%Purpose:
%Back-to-back bar plot for Model M3 simulation output 
%Plot multiple replicates at once.
%Produces two figures
%(i) Posterior predictive distributions for influenza positive GP consultations per 100,000 population.
%(ii) Posterior predictive distributions for circulating influenza virus subtype/lineage composition.

%Authors: Matt Keeling & Ed Hill
%--------------------------------------------------------------------------

clear variables

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
ConstructBarPlot(SynthDataFlag,MultiSeasonImmFlag,SeasonsToPlot)

function ConstructBarPlot(SynthDataFlag,MultiSeasonImmFlag,SeasonsToPlot)
%   SynthDataFlag - (indicator variable) 0: empirical data; 1: synthetic data
%   MultiSeasonImmFlag - (indicator variable) Model type used. 
%                                             0: Original, immunity propagation lasts a single season; 
%                                             1: Extended, allowing immunity to be propagated for more than a single season
%   SeasonsToPlot - (scalar, integer) Number of influenza seasons worth of data fit to in the inference scheme 



%% Load required simulation data
if SynthDataFlag == 0 %use outputs from fitting to empirical data
    if SeasonsToPlot == 4
        InputDataFileName = 'ModelSimnData/ModelSimns_FourSeasonFit.mat';
    elseif SeasonsToPlot == 5
        InputDataFileName = 'ModelSimnData/ModelSimns_FiveSeasonFit.mat';
    elseif SeasonsToPlot == 6        
        if MultiSeasonImmFlag == 0
            InputDataFileName = 'ModelSimnData/ModelSimns_SixSeasonFit.mat';
        else
            InputDataFileName = 'ModelSimnData/ExtendedModelSimns_SixSeasonFit.mat';
        end
    else
        error('Misspecified SeasonsToPlot value of %f. SeasonsToPlot must take a value of 4,5 or 6',SeasonsToPlot)
    end
    
    %Empirical data
    EmpDataOffset = 6 - SeasonsToPlot;
    Truth=load('../../Data/ILIData/EmpData_InfluenzaPositiveGPConsultRateByStrain_2009to2018.csv');
    TruthRetained = Truth(4:end-EmpDataOffset,:);
elseif SynthDataFlag == 1 %use outputs from fitting to synthetic data
    
    %Load required simulation data and synthetic data
    if MultiSeasonImmFlag == 0
        InputDataFileName = 'ModelSimnData/ModelSimns_SynthDataFit.mat';
        Truth=load('../../Data/ILIData/SynthData_InfluenzaPositiveGPConsultRateByStrain_2009to2018.csv');
    else
        InputDataFileName = 'ModelSimnData/ExtendedModelSimns_SynthDataFit.mat';
        Truth=load('../../Data/ILIData/ExtendedModel_SynthData_InfluenzaPositiveGPConsultRateByStrain_2009to2018.csv');
    end
    
    %Load synthetic data
    TruthRetained = Truth(4:end,:);
    
end

load(InputDataFileName,'SimnData')
Simn_CellOutput = SimnData(:,1); %From data file, get array outputs from each season

%Put into 3D array
NumOfStrains = 4;
SimnNum = 1000;
ModelM3Simn_GPFluConsultOutput = zeros(SeasonsToPlot,NumOfStrains,SimnNum);

for ii=1:SimnNum
    ModelM3Simn_GPFluConsultOutput(:,:,ii) = Simn_CellOutput{ii};
end


%% Influenza attributable GP consultations: Stacked bar plot with empirical data plotted as smaller bar above the simulated distribution figure elements

fig = figure(); clf;
set(fig,'Color', [1 1 1])
position = [100, 100, 2.5*550, 1.7*450];
set(0, 'DefaultFigurePosition', position);
hold on;

YLimMax = SeasonsToPlot+1.5;
plot([0 0],[0.8 YLimMax],'-k'); %Have central dividing line

for Year=1:SeasonsToPlot;
    clear R;
    R(1:NumOfStrains,1:SimnNum)=ModelM3Simn_GPFluConsultOutput(Year,:,:);
    AH1=R(1,:); AH3=R(2,:); BYam=R(3,:); BVic=R(4,:);
    AH1_Truth = TruthRetained(Year,1); AH3_Truth = TruthRetained(Year,2); %Empirical type A data
    BYam_Truth = TruthRetained(Year,3); BVic_Truth = TruthRetained(Year,4); %Empirical type B data
    y = 0.1 + 0.7*[0:SimnNum-1]/(SimnNum-1);
    o1=[1:SimnNum]; o2=[1:SimnNum];
    % if you include the next line the two sides of the plot are independent,
    % if you comment it out the figure is messier, but a more acurate
    % representation.
    %[x,o]=sort(BYam+BVic); o2=o(q);
    
    h1=fill([0 -AH1_Truth-AH3_Truth -AH1_Truth-AH3_Truth 0],Year+[0.85 0.85 0.95 0.95],'r','DisplayName','A(H1N1)pdm09'); set(h1,'FaceColor',[1 0 0]);
    h2=fill([0 -AH3_Truth -AH3_Truth 0],Year+[0.85 0.85 0.95 0.95],'r','DisplayName','A(H3N2)'); set(h2,'FaceColor',[1 0.6 0]);
    h3=fill([0 BYam_Truth+BVic_Truth BYam_Truth+BVic_Truth 0],Year+[0.85 0.85 0.95 0.95],'b','DisplayName','B/Yamagata'); set(h3,'FaceColor',[0 0 1]);
    h4=fill([0 BVic_Truth BVic_Truth 0],Year+[0.85 0.85 0.95 0.95],'b','DisplayName','B/Victoria'); set(h4,'FaceColor','c');

    h=fill([0 -AH1(o1)-AH3(o1) 0],Year+[0.1 y 0.8],'r'); set(h,'FaceColor',[1 0 0]); set(h,'EdgeColor',[0.99 0.99 0.99],'LineStyle','none');
    h=fill([0 -AH3(o1) 0],Year+[0.1 y 0.8],'r'); set(h,'FaceColor',[1 0.7 0],'LineStyle','none');
    h=fill([0 BYam(o2)+BVic(o2) 0],Year+[0.1 y 0.8],'b'); set(h,'LineStyle','none'); set(h,'FaceColor',[0.4 0.4 1]);
    h=fill([0 BVic(o2) 0],Year+[0.1 y 0.8],'c'); set(h,'LineStyle','none'); 
end
hold off;
set(gca,'YLim',[0.8 YLimMax]);

%Add in legend
leg1 = legend([h1;h2;h4;h3]);
set(leg1,'Location','eastoutside')
%set(leg1,'Visible','off') %Make invisible while retaining whitespace beside the plot

%Set y-axis tick length to zero
ax = gca;
ax.YAxis.TickLength = [0 0];

%Set y-axis limits and label
if SeasonsToPlot == 3
    yticks([1.5 2.5 3.5])
    yticklabels({'2012/13','2013/14','2014/15'});
elseif SeasonsToPlot == 4
    yticks([1.5 2.5 3.5 4.5])
    yticklabels({'2012/13','2013/14','2014/15','2015/16'});
elseif SeasonsToPlot == 5
    yticks([1.5 2.5 3.5 4.5 5.5])
    yticklabels({'2012/13','2013/14','2014/15','2015/16','2016/17'});
elseif SeasonsToPlot == 6
    yticks([1.5 2.5 3.5 4.5 5.5 6.5])
    yticklabels({'2012/13','2013/14','2014/15','2015/16','2016/17','2017/18'});
end

%Set x-axis limits and tick labels
if SeasonsToPlot == 4 || SeasonsToPlot == 5
    xlim([-150 150])
    xticks([-150,-100,-50,0,50,100,150])
    xticklabels({'150','100','50','0','50','100','150'})
    
    %Add influenza type descriptive text x-position
    TextLabelXPos1 = -75;
    TextLabelXPos2 = 75;
elseif SeasonsToPlot == 6  

    xlim([-200 200])
    xticks([-200,-150,-100,-50,0,50,100,150,200])
    xticklabels({'200','150','100','50','0','50','100','150','200'})
    
    %Add influenza type descriptive text x-position
    TextLabelXPos1 = -130;
    TextLabelXPos2 = 70;
end

%Add influenza type descriptive text
TextLabelYPos = YLimMax - 0.3;
txt1 = 'Influenza A';
text(TextLabelXPos1,TextLabelYPos,txt1,'FontWeight','bold','FontSize',18)

txt2 = 'Influenza B';
text(TextLabelXPos2,TextLabelYPos,txt2,'FontWeight','bold','FontSize',18)


%Axes labels
xlabel('GP consultations attributable to influenza (per 100,000)')
ylabel('Influenza Season')

%Title
title('Comparison to GP consultation data')

%Specify general axis properties
set(gca,'FontSize',18)
set(gca,'LineWidth',1)
box on

%% Influenza subtype/lineage composition: Stacked bar plot with empirical data plotted as smaller bar above the simulated distribution figure elements

fig = figure(); clf;
set(fig,'Color', [1 1 1])
position = [100, 100, 2.5*550, 1.7*450];
set(0, 'DefaultFigurePosition', position);
hold on;


YLimMax = SeasonsToPlot+1.5;
plot([0 0],[0.8 YLimMax],'-k'); %Have central dividing line

%For each season, sum across all strains Truth data and Simn data
%Divide individual strain contribution by overall sum to get the proportion
for Year=1:SeasonsToPlot;
    clear R;
    
    %Get GP consultation data per run
    R(1:NumOfStrains,1:SimnNum)=ModelM3Simn_GPFluConsultOutput(Year,:,:);
    AH1=R(1,:); AH3=R(2,:); BYam=R(3,:); BVic=R(4,:);
    AH1_Truth = TruthRetained(Year,1); AH3_Truth = TruthRetained(Year,2); %Empirical type A data
    BYam_Truth = TruthRetained(Year,3); BVic_Truth = TruthRetained(Year,4); %Empirical type B data
    y = 0.1 + 0.7*[0:SimnNum-1]/(SimnNum-1);
    o1=[1:SimnNum]; o2=[1:SimnNum];
    
    %Scale Truth data so it becomes a proportion
    TruthSum = AH1_Truth + AH3_Truth + BYam_Truth + BVic_Truth;
    AH1_TruthPropn = AH1_Truth/TruthSum; AH3_TruthPropn = AH3_Truth/TruthSum; %Empirical type A strain composition
    BYam_TruthPropn = BYam_Truth/TruthSum; BVic_TruthPropn = BVic_Truth/TruthSum; %Empirical type B strain composition
    
    %Get denominator for scaling model simulated data to proportions
    PredSum = AH1(o1) + AH3(o1) + BYam(o2) + BVic(o2);
    
    %Produce bar plots
    h1=fill([0 -AH1_TruthPropn-AH3_TruthPropn -AH1_TruthPropn-AH3_TruthPropn 0],Year+[0.85 0.85 0.95 0.95],'r','DisplayName','A(H1N1)pdm09'); set(h1,'FaceColor',[1 0 0]);
    h2=fill([0 -AH3_TruthPropn -AH3_TruthPropn 0],Year+[0.85 0.85 0.95 0.95],'r','DisplayName','A(H3N2)'); set(h2,'FaceColor',[1 0.6 0]);
    h3=fill([0 BYam_TruthPropn+BVic_TruthPropn BYam_TruthPropn+BVic_TruthPropn 0],Year+[0.85 0.85 0.95 0.95],'b','DisplayName','B/Yamagata'); set(h3,'FaceColor',[0 0 1]);
    h4=fill([0 BVic_TruthPropn BVic_TruthPropn 0],Year+[0.85 0.85 0.95 0.95],'b','DisplayName','B/Victoria'); set(h4,'FaceColor','c');

    h=fill([0 (-AH1(o1)-AH3(o1))./PredSum 0],Year+[0.1 y 0.8],'r'); set(h,'FaceColor',[1 0 0]); set(h,'EdgeColor',[0.99 0.99 0.99],'LineStyle','none');
    h=fill([0 -AH3(o1)./PredSum 0],Year+[0.1 y 0.8],'r'); set(h,'FaceColor',[1 0.7 0],'LineStyle','none');
    h=fill([0 (BYam(o2)+BVic(o2))./PredSum 0],Year+[0.1 y 0.8],'b'); set(h,'LineStyle','none'); set(h,'FaceColor',[0.4 0.4 1]);
    h=fill([0 BVic(o2)./PredSum 0],Year+[0.1 y 0.8],'c'); set(h,'LineStyle','none'); 
end
hold off;
set(gca,'YLim',[0.8 YLimMax]);

%Set y-axis tick length to zero
ax = gca;
ax.YAxis.TickLength = [0 0];

%Set y-axis limits and label
if SeasonsToPlot == 3
    yticks([1.5 2.5 3.5])
    yticklabels({'2012/13','2013/14','2014/15'});
elseif SeasonsToPlot == 4
    yticks([1.5 2.5 3.5 4.5])
    yticklabels({'2012/13','2013/14','2014/15','2015/16'});
elseif SeasonsToPlot == 5
    yticks([1.5 2.5 3.5 4.5 5.5])
    yticklabels({'2012/13','2013/14','2014/15','2015/16','2016/17'});
elseif SeasonsToPlot == 6
    yticks([1.5 2.5 3.5 4.5 5.5 6.5])
    yticklabels({'2012/13','2013/14','2014/15','2015/16','2016/17','2017/18'});
end

%Set x-axis limits and tick labels
xlim([-1 1])
xticks([-1,-0.75,-0.5,-0.25,0,0.25,0.50,0.75,1])
xticklabels({'1','0.75','0.50','0.25','0','0.25','0.50','0.75','1'})

%Add influenza type descriptive text x-position
TextLabelXPos1 = -0.65;
TextLabelXPos2 = 0.35;

%Add influenza type descriptive text
TextLabelYPos = YLimMax - 0.3;
txt1 = 'Influenza A';
text(TextLabelXPos1,TextLabelYPos,txt1,'FontWeight','bold','FontSize',20)

txt2 = 'Influenza B';
text(TextLabelXPos2,TextLabelYPos,txt2,'FontWeight','bold','FontSize',20)

%Add a placeholder legend
leg1 = legend([h1;h2;h4;h3]);
set(leg1,'Location','eastoutside')

%Axes labels
xlabel('Proportion per influenza type')
ylabel('Influenza Season')

%Plot title
title('Comparison to empirical strain distribution')

%Specify general axis properties
set(gca,'FontSize',20)
set(gca,'LineWidth',1)
box on

end

