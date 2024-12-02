function SlidingTimeDecodeAction(AllTrialData,Labels,DecodeLabels,TimeWindow,TimeStep,StepSize,OutDir,Tag,varargin)
% Plots confusion matrix (decode analysis) for decoding each (relevant)
% condition.

% [~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh',0.001);
[~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh','fdr');
[~,NIter] = Utilities.ProcVarargin(varargin,'NIter',100);
[~,Phase] = Utilities.ProcVarargin(varargin,'Phase','Go');
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[~,reorder]=Utilities.ProcVarargin(varargin,'reorder',[]);
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');

% ASF = Analyze.FaceScratch3.FindSignificantUnits(Tag,'PThresh',PThresh,'Labels',Labels,'Labels2Incl',DecodeLabels);
ASF = Analyze.FaceScratch3.FindSignificantUnitsAction(Tag,'PThresh',PThresh,'Labels',Labels,'Labels2Incl',DecodeLabels,'All');
% ASF = Analyze.returnUniqueFieldValues(AllTrialData,'UnitIDsSorted');

Conditions={{'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek'},...
    {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'},...
    {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'},...
    {'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'}};

counter = 1;

for i=TimeWindow(1):TimeStep:TimeWindow(end)
    
    cTW = [i i+StepSize];
    TimeCenters(1,counter) = i+TimeStep/2;
    
    for jj = 1:length(Conditions)
        
        cLabels = Conditions{jj};
        cDecodeLabels = cLabels;
        
        cvAcc{jj}(:,counter) = Decode(AllTrialData,cLabels,cDecodeLabels,cTW,OutDir,Tag,ASF);
    end
    counter = counter+1;
    
end

%% to correct the offset for a timestep of 0.15 and stepsize of 0.3

tmpcv;
tmpTime;

cvAcc=tmpcv;
TimeCenters=tmpTime;
%
cvAcc{1}(:,1)=[];
cvAcc{4}(:,1)=[];
cvAcc{3}(:,[1 2])=[];
% tmpcv{4}(:,1)=[];
TimeCenters(:,[end end-1])=[];

cvAcc{1}(:,end)=[];
cvAcc{2}(:,[end end-1])=[]; 
% cvAcc{3}(:,1)=[]; 
cvAcc{4}(:,[end])=[];

cvAcc{1}(:,[end])=[];
cvAcc{2}(:,end)=[]; 
% cvAcc{3}(:,end)=[]; 
cvAcc{4}(:,end)=[];

%% plot the cvAccuracy over time

StyleArgs = {'FontName','Arial','FontSize',9,'FontWeight','bold'};

figH = plt.fig('units','inches','width',6,'height',3.5,'font','Arial','fontsize',9);
figure(figH);
clf; pnl = panel(); pnl.margin=15; pnl.pack(1,2);

% clr=[133 133 133]./255;
clr=[];
clr(1,:) = repmat([237 0 38]./255,1,1);
clr(2,:) = repmat([31 120 180]./255,1,1);
clr(3,:) = repmat([240 228 66]./255,1,1);
clr(4,:) = repmat([189 189 189]./255,1,1);

cLabel = {{'Felt Cheek'},{'Felt Shoulder'},{'Obs Cheek'},{'Obs Shoulder'}};

% plot the discriminability of touch from baseline
pnl(1,1).select()
hold on;

for i=1:4
    
    tmp{1} = cvAcc{i};
    Analyze.plotEventRelatedAverage(tmp,cLabel{i},'TimeVec',TimeCenters,...
        'Colors',clr(i,:),'useBootStrap');
end
% plot the legend, remove the patches, plot the text in the corresponding
% color
[l,icons,plots] = legend([cLabel{:}],'Location','South','box','off','linewidth',0.01,'FontName','Arial','FontWeight','bold','FontSize',9);
for i=1:4
    icons(i+4).Visible = 'off';
    icons(i).Color = clr(i,:);
end

ax1 = gca;
set(ax1,'box','off',StyleArgs{:});
ax1.XRuler.Axle.LineStyle = 'none';
set(gca,'TickDir','out');
set(ax1,'TickLength',[0 0]);
axis tight;
xlabel('Time (s)',StyleArgs{:});
ylabel('Classification Accuracy (%)',StyleArgs{:});
set(ax1,'XTick',[TimeWindow(1):1:TimeWindow(end)],'XTickLabel',[TimeWindow(1):1:TimeWindow(end)]);
% NC = 5; xlim([ax1.XLim(1) NC+.25]);
ylim([0 105])

tEvents = [0 3];
plt.vline(tEvents,{'-','LineWidth',1.25,'Color',[133 133 133]./255});
ylims = ylim();
ytext = ylims(2) + .05*diff(ylims);
text(tEvents(1),ytext-2,'Go',StyleArgs{:});

line([ax1.XLim(1) ax1.XLim(end)],[25 25],'Color',[115 115 115]./255,'LineStyle','--','LineWidth',1.25);
text(ax1.XLim(end),24,'Chance','HorizontalAlignment','right','VerticalAlignment','top',...
    'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[115 115 115]./255);

FigsDir=fullfile(OutDir, 'PopData', 'Mine');
plt.SaveFigure(saveFig,FigsDir,['SlidingTimeDecode'],'SVG','PNG');
% plt.SaveFigure(saveFig,FigsDir,['SlidingTimeDecode'],'SVGI');

end

function cvAcc = Decode(AllTrialData,Labels,DecodeLabels,TimeWindow,OutDir,Tag,ASF,varargin)

[~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh','fdr');
[~,NIter] = Utilities.ProcVarargin(varargin,'NIter',100);
[~,Phase] = Utilities.ProcVarargin(varargin,'Phase','Go');
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[~,reorder]=Utilities.ProcVarargin(varargin,'reorder',[]);
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');


NConds = length(DecodeLabels);

%% Only include units with tunign to at least one of the relevant variables
% TODO: Dropping curve version.

NObsPerCond = 10;

if ~isempty(reorder)
    DecodeLabels = DecodeLabels(reorder);
end

UDates = Analyze.returnUniqueFieldValues(AllTrialData,'Date');
for d = 1:length(UDates)
    fprintf('Day %d/%d\n',d,length(UDates));
    cTrialData = Analyze.SubSelectTrials(AllTrialData,'Date',UDates{d},'Phase',Phase,'Condition',DecodeLabels);
    FR = squeeze(Analyze.getNeuralData(cTrialData,ASF,TimeWindow));
    FR(:,any(isnan(FR))) = [];
    Label = Analyze.returnFieldValues(cTrialData,'Condition');
    
    % Use ismember instead of unique so that label order is same as DecodeLabels
    [~,LabelId] = ismember(Label,DecodeLabels);
    %     [ULabels,~,LabelId]= unique(Label);
    
    ConfMat{d} = zeros(NConds,NConds);
    [FR2,Labels2] = Analyze.SampleFeaturePopulation(FR,LabelId,'Type','Null','NumObservationsPerCondition',NObsPerCond);
    
    % Leave one out cross-val
    %     for k = 1:size(FR2,1)
    %         testIdx = k;
    %         trainIdx = setdiff(1:size(FR2,1),k);
    %         mdl = fitcdiscr(FR2(trainIdx,:),Labels2(trainIdx),'discrimType','diagLinear');
    %         predLabels = predict(mdl,FR2(testIdx,:));
    %         ConfMat{d}(Labels2(testIdx),predLabels) = ConfMat{d}(Labels2(testIdx),predLabels) + 1;
    %     end
    for k = 1:NObsPerCond
        testIdx = k:NObsPerCond:size(FR2,1);
        trainIdx = setdiff(1:size(FR2,1),testIdx);
        mdl = fitcdiscr(FR2(trainIdx,:),Labels2(trainIdx),'discrimType','diagLinear');
        predLabels = predict(mdl,FR2(testIdx,:));
        for i = 1:length(predLabels)
            ConfMat{d}(Labels2(testIdx(i)),predLabels(i)) = ConfMat{d}(Labels2(testIdx(i)),predLabels(i)) + 1;
        end
    end
    ConfMat{d} = ConfMat{d}/(size(FR2,1)/NConds);
end

%% gather the cvAcc to send back to main function

cvAcc = [];
tmp = cellfun(@(x)diag(x),ConfMat,'UniformOutput',false);
cvAcc = vertcat(tmp{:})*100;

end