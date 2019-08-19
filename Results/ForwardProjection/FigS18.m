%Purpose:
%Script to genereate Fig. S18 from the publication "Seasonal influenza: Modelling 
%approaches to capture immunity propagation"

%Forward simulation
%Plot to visualise the variability per season (amongst runs with vaccine
%efficacy sampled from empirical distribution)
%with the distribution/percentage of runs where one subtype dominates

%Author: Ed Hill
%--------------------------------------------------------------------------

%% Load randomly sampled vaccine efficacy scenario data
InputData2 = load('ForwardProjSimns_SampledVaccEff.mat');
RandomVaccEff_CellOutput = InputData2.SimnData; %From data file, get array outputs from each season


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

%% Randomly sampled vaccine efficacy scenario analysis

%Influenza A strain compositon (by season)
AllFluA = sum(ModelM3Simn_ForwardProjFluCaseOutput(:,1:2,:),2); %Flu A total, by season and simulation run
H1Propn = ModelM3Simn_ForwardProjFluCaseOutput(:,1,:)./AllFluA;
H3Propn = ModelM3Simn_ForwardProjFluCaseOutput(:,2,:)./AllFluA;

%Iterate through each season
%Get number of simulations where A(H1N1) dominated A(H3N2)
SeasonNum = size(AllFluA,1);
SimnNum = size(AllFluA,3);
H1DomFrac = zeros(1,SeasonNum);
for ii = 1:SeasonNum
    H1DomCount = sum(H1Propn(ii,:,:) > 0.5);
    H1DomFrac(ii) = H1DomCount/SimnNum;
end

H3DomFrac = 1 - H1DomFrac; %H3 dominates in remaining fraction of runs 

%Influenza B strain compositon (by season)
AllFluB = sum(ModelM3Simn_ForwardProjFluCaseOutput(:,3:4,:),2); %Flu B total, by season and simulation run
YamPropn = ModelM3Simn_ForwardProjFluCaseOutput(:,3,:)./AllFluB;
VicPropn = ModelM3Simn_ForwardProjFluCaseOutput(:,4,:)./AllFluB;

%Iterate through each season
%Get number of simulations where B/Yam dominated B/Vic
SeasonNum = size(AllFluB,1);
SimnNum = size(AllFluB,3);
YamDomFrac = zeros(1,SeasonNum);
for ii = 1:SeasonNum
    YamDomCount = sum(YamPropn(ii,:,:) > 0.5);
    YamDomFrac(ii) = YamDomCount/SimnNum;
end

VicDomFrac = 1 - YamDomFrac; %Victoria lineage dominates in remaining fraction of runs

%% Construct back-to-back bar plot
%% Per season, bar showing proportion of runs where that subtype/lineage was dominant 
%% (>50% of the overall incidence for that type of influenza)

fig = figure(); clf;
set(fig,'Color', [1 1 1])
position = [100, 100, 2.5*550, 2.0*450];
set(0, 'DefaultFigurePosition', position);
hold on;

YLimMax = SeasonsToPlot+1.5;
plot([0 0],[0.8 YLimMax],'-k'); %Have central dividing line

%Plot bars
for Year=1:SeasonsToPlot
    
    %Assign stats for current season index to variables
    AH1_DomFracInSeason = H1DomFrac(Year);
    AH3_DomFracInSeason = H3DomFrac(Year);

    Vic_DomFracInSeason = VicDomFrac(Year);
    Yam_DomFracInSeason = YamDomFrac(Year);

    %Plot outcome from single replicate, particle set that gave lowest
    %error
    h1=fill([0 -AH1_DomFracInSeason -AH1_DomFracInSeason 0],Year+[0.55 0.55 0.85 0.85],'r','DisplayName','A(H1N1)pdm09'); set(h1,'FaceColor',[1 0 0]);
    h2=fill([0 -AH3_DomFracInSeason -AH3_DomFracInSeason 0],Year+[0.15 0.15 0.45 0.45],'r','DisplayName','A(H3N2)'); set(h2,'FaceColor',[1 0.6 0]);
    h3=fill([0 Yam_DomFracInSeason Yam_DomFracInSeason 0],Year+[0.15 0.15 0.45 0.45],'b','DisplayName','B/Yamagata'); set(h3,'FaceColor',[0 0 1]);
    h4=fill([0 Vic_DomFracInSeason Vic_DomFracInSeason 0],Year+[0.55 0.55 0.85 0.85],'b','DisplayName','B/Victoria'); set(h4,'FaceColor','c');

end

%Add in legend
leg1 = legend([h1;h2;h4;h3]);
set(leg1,'Location','eastoutside')

%Set x-axis limits and tick labels
xlim([-1 1])
xticks([-1,-0.75,-0.5,-0.25,0,0.25,0.50,0.75,1])
xticklabels({'1','0.75','0.50','0.25','0','0.25','0.50','0.75','1'})

%Set y-axis limits and label
ylim([1 13])
yticks([1.5:1:12.5])
yticklabels({'2018/19','2019/20','2020/21','2021/22','2022/23','2023/24',...
             '2024/25','2025/26','2026/27','2027/28','2028/29','2029/30'});

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
%set(leg1,'Visible','off') %Make invisible while retaining whitespace beside the plot

%Axes labels
xlabel('Dominating subtype/lineage (proportion of total runs)')
ylabel('Influenza Season')

%Specify general axis properties
set(gca,'FontSize',20)
set(gca,'LineWidth',1)
box on
