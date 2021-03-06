%Purpose:
%Script to genereate Fig. 7 from the publication "Seasonal influenza: Modelling 
%approaches to capture immunity propagation"


%Plot predicted influenza positive GP consultation counts for next decade
%Compare high efficacy scenario (bar topping each season) to replicates
%with vaccine efficacy sampled from the empirical distribution

%Author: Ed Hill
%--------------------------------------------------------------------------

clear variables

%% Load median efficacy scenario data
InputData1 = load('ForwardProjSimns_MedVaccEff.mat');
MedVaccEffData = InputData1.SimnData; %From data file, get array outputs from each season

%Scale proportions from final 12 seasons (covering 2018/19 to 2029/30) to a 
% rateper 100,000 population
MedVaccEffEstCases = MedVaccEffData(end-11:end,:)*100000; 

%% Load low efficacy scenario data
InputData2 = load('ForwardProjSimns_MinVaccEff.mat');
LowVaccEffData = InputData2.SimnData; %From data file, get array outputs from each season

%Scale proportions from final 12 seasons (covering 2018/19 to 2029/30) to a 
% rateper 100,000 population
LowVaccEffEstCases = LowVaccEffData(end-11:end,:)*100000; 

%% Load high efficacy scenario data
InputData3 = load('ForwardProjSimns_MaxVaccEff.mat');
HighVaccEffData = InputData3.SimnData; %From data file, get array outputs from each season

%Scale proportions from final 12 seasons (covering 2018/19 to 2029/30) to a 
% rateper 100,000 population
HighVaccEffEstCases = HighVaccEffData(end-11:end,:)*100000; 


%% Load randomly sampled vaccine efficacy scenario data
InputData4 = load('ForwardProjSimns_SampledVaccEff.mat');
RandomVaccEff_CellOutput = InputData4.SimnData; %From data file, get array outputs from each season


%% Retain forward projected seasons through to 2029/2030.
%Put into 3D array
SeasonsToPlot = 12;
NumOfStrains = 4;
SimnNum = 1000;
ModelM3Simn_ForwardProjFluCaseOutput = zeros(SeasonsToPlot,NumOfStrains,SimnNum);
for ii=1:SimnNum
    %Pick out final 12 rows, corrsponding to the period 2018/2019 to 2029/2030
    %(inclusive)
    ModelM3Simn_ForwardProjFluCaseOutput(:,:,ii) = RandomVaccEff_CellOutput{ii}(end-11:end,:)*100000;
   
    %Scale by 100,000, to give a standardised reporting measure
end


%% Produce top panel plot
%% Model simulation replicates, for both influenza types, ordered by case size

for i=1:(SimnNum/2)
    q(i)=(i-1)*2+1;
    q(SimnNum+1-i)=i*2;
end

fig = figure(); clf;
set(fig,'Color', [1 1 1])
position = [100, 100, 2.5*550, 2.0*450];
set(0, 'DefaultFigurePosition', position);
hold on

%Have central dividing line
YLimMax = SeasonsToPlot+1.5;
plot([0 0],[0.8 YLimMax],'-k'); 


for Year=1:SeasonsToPlot
    clear R;
    
    %Get high vacc. eff. scenario data
    AH1_MedVaccEff = MedVaccEffEstCases(Year,1);
    AH3_MedVaccEff = MedVaccEffEstCases(Year,2);
    BYam_MedVaccEff = MedVaccEffEstCases(Year,3);
    BVic_MedVaccEff = MedVaccEffEstCases(Year,4);

    %Get randomly sampled vacc. eff. scenario data
    R(1:NumOfStrains,1:SimnNum) = ModelM3Simn_ForwardProjFluCaseOutput(Year,:,:);
    AH1=R(1,:); AH3=R(2,:); BYam=R(3,:); BVic=R(4,:);
    
    %Set up plot variables
    y=0.1 + 0.7*[0:SimnNum-1]/(SimnNum-1);
    %o1=[1:SimnNum]; o2=[1:SimnNum];

    % if you include the next line the two sides of the plot are independent,
    % if you comment it out the figure is messier, but a more acurate
    % representation.
    [x,o]=sort(AH1+AH3); o1=o(q);
    [x,o]=sort(BYam+BVic); o2=o(q);
    
    %Plot outcome from single replicate, particle set that gave lowest
    %error
    h1=fill([0 -AH1_MedVaccEff-AH3_MedVaccEff -AH1_MedVaccEff-AH3_MedVaccEff 0],Year+[0.85 0.85 0.95 0.95],'r','DisplayName','A(H1N1)pdm09'); set(h1,'FaceColor',[1 0 0]);
    h2=fill([0 -AH3_MedVaccEff -AH3_MedVaccEff 0],Year+[0.85 0.85 0.95 0.95],'r','DisplayName','A(H3N2)'); set(h2,'FaceColor',[1 0.6 0]);
    h3=fill([0 BYam_MedVaccEff+BVic_MedVaccEff BYam_MedVaccEff+BVic_MedVaccEff 0],Year+[0.85 0.85 0.95 0.95],'b','DisplayName','B/Yamagata'); set(h3,'FaceColor',[0 0 1]);
    h4=fill([0 BVic_MedVaccEff BVic_MedVaccEff 0],Year+[0.85 0.85 0.95 0.95],'b','DisplayName','B/Victoria'); set(h4,'FaceColor','c');

    %Plot outcomes from all 1,000 particle set runs
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

%Add influenza type descriptive text
TextLabelYPos = YLimMax - 0.3;
TextLabelXPos1 = -4e4;
txt1 = 'Influenza A';
text(TextLabelXPos1,TextLabelYPos,txt1,'FontWeight','bold','FontSize',18)

TextLabelXPos2 = 2e4;
txt2 = 'Influenza B';
text(TextLabelXPos2,TextLabelYPos,txt2,'FontWeight','bold','FontSize',18)

%Set axes tick lengths
set(gca,'TickLength',[0.003, 0])

%Set y-axis limits and label
yticks([1.5:1:12.5])
yticklabels({'2018/19','2019/20','2020/21','2021/22','2022/23','2023/24',...
             '2024/25','2025/26','2026/27','2027/28','2028/29','2029/30'});

% %Set x-axis limits and tick labels
xlim([-6.1e4 6.1e4])
xticks([-6e4,-4e4,-2e4,0,2e4,4e4,6e4])
xticklabels({'60,000','40,000','20,000','0','20,000','40,000','60,000'})

%Axes labels
xlabel('Projected influenza seasonal incidence (per 100,000)')
ylabel('Influenza Season')

%Plot title
title('Forward projections: ''Expected'' scenario & empirically sampled efficacy')

%Specify general axis properties
set(gca,'FontSize',18)
set(gca,'LineWidth',1)
box on

%% Produce panels (b) and (c)

% Aggregate case data into tuple variable
ScenVaccEffEstCasesTuple = {LowVaccEffEstCases,HighVaccEffEstCases};
ScenNum = numel(ScenVaccEffEstCasesTuple);

% Specify scenario-specific plot variables
TitleString = {'Pessimistic scenario', 'Optimistic scenario'};

% Generate bar plot displaying incidence per strain type for each vaccine efficacy scenario)
for ScenIdx = 1:ScenNum
    
    %Extract relevant data from case tuple
    ScenVaccEffEstCases = ScenVaccEffEstCasesTuple{ScenIdx};
    
    %Initialise figure
    fig = figure(); clf;
    set(fig,'Color', [1 1 1])
    position = [100, 100, 1.5*550, 1.5*450];
    set(0, 'DefaultFigurePosition', position);
    hold on
    
    %Have central dividing line
    YLimMax = SeasonsToPlot+1.5;
    plot([0 0],[0.8 YLimMax],'-k');
    
    for Year=1:SeasonsToPlot
        clear R;
        
        %Get fixed vacc. eff. scenario data
        AH1_VaccEff = ScenVaccEffEstCases(Year,1);
        AH3_VaccEff = ScenVaccEffEstCases(Year,2);
        BYam_VaccEff = ScenVaccEffEstCases(Year,3);
        BVic_VaccEff = ScenVaccEffEstCases(Year,4);
 
        %Plot outcome from single replicate, particle set that gave lowest
        %error
        h1=fill([0 -AH1_VaccEff-AH3_VaccEff -AH1_VaccEff-AH3_VaccEff 0],Year+[0.2 0.2 0.8 0.8],'r','DisplayName','A(H1N1)pdm09'); set(h1,'FaceColor',[1 0 0]);
        h2=fill([0 -AH3_VaccEff -AH3_VaccEff 0],Year+[0.2 0.2 0.8 0.8],'r','DisplayName','A(H3N2)'); set(h2,'FaceColor',[1 0.6 0]);
        h3=fill([0 BYam_VaccEff+BVic_VaccEff BYam_VaccEff+BVic_VaccEff 0],Year+[0.2 0.2 0.8 0.8],'b','DisplayName','B/Yamagata'); set(h3,'FaceColor',[0 0 1]);
        h4=fill([0 BVic_VaccEff BVic_VaccEff 0],Year+[0.2 0.2 0.8 0.8],'b','DisplayName','B/Victoria'); set(h4,'FaceColor','c');
        
    end
    hold off;
    set(gca,'YLim',[0.8 YLimMax]);
    
%     %Add in legend
%     leg1 = legend([h1;h2;h4;h3]);
%     set(leg1,'Location','eastoutside')
%     
    %Add influenza type descriptive text
    TextLabelYPos = YLimMax - 0.3;
    TextLabelXPos1 = -4e4;
    txt1 = 'Influenza A';
    text(TextLabelXPos1,TextLabelYPos,txt1,'FontWeight','bold','FontSize',18)
    
    TextLabelXPos2 = 2e4;
    txt2 = 'Influenza B';
    text(TextLabelXPos2,TextLabelYPos,txt2,'FontWeight','bold','FontSize',18)
    
    %Set axes tick lengths
    set(gca,'TickLength',[0.003, 0])
    
    %Set y-axis limits and label
    yticks([1.5:1:12.5])
    yticklabels({'2018/19','2019/20','2020/21','2021/22','2022/23','2023/24',...
        '2024/25','2025/26','2026/27','2027/28','2028/29','2029/30'});
    
    % %Set x-axis limits and tick labels
    xlim([-6.1e4 6.1e4])
    xticks([-6e4,-4e4,-2e4,0,2e4,4e4,6e4])
    xticklabels({'60,000','40,000','20,000','0','20,000','40,000','60,000'})
    
    %Axes labels
    xlabel('Projected influenza seasonal incidence (per 100,000)')
    ylabel('Influenza Season')
    
    %Plot title
    title(TitleString{ScenIdx})
    %title 'Forward projection (randomly sampled vaccine efficacy)'
    
    %Specify general axis properties
    set(gca,'FontSize',18)
    set(gca,'LineWidth',1)
    box on
end