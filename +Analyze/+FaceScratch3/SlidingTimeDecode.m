function SlidingTimeDecode(AllTrialData,Labels,DecodeLabels,TimeWindow,TimeStep,StepSize,OutDir,Tag,varargin)
% Plots confusion matrix (decode analysis) for decoding each (relevant)
% condition.

% [~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh',0.001);
[~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh','fdr');
[~,NIter] = Utilities.ProcVarargin(varargin,'NIter',100);
[~,Phase] = Utilities.ProcVarargin(varargin,'Phase','Go');
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[~,reorder]=Utilities.ProcVarargin(varargin,'reorder',[]);
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');

ASF = Analyze.FaceScratch3.FindSignificantUnits(Tag,'PThresh',PThresh,'Labels',Labels,'Labels2Incl',DecodeLabels);
% ASF = Analyze.returnUniqueFieldValues(AllTrialData,'UnitIDsSorted');

counter = 1;

for i=TimeWindow(1):TimeStep:TimeWindow(end)
    
    cTW = [i i+StepSize];
    TimeCenters(1,counter) = i+TimeStep/2;
    
    cvAcc{1}(:,counter) = Decode(AllTrialData,Labels,DecodeLabels,cTW,OutDir,Tag,ASF);
    counter = counter+1;
end

%% plot the cvAccuracy over time

StyleArgs = {'FontName','Arial','FontSize',9,'FontWeight','bold'};

figH = plt.fig('units','inches','width',10,'height',5,'font','Arial','fontsize',9);
figure(figH);
clf; pnl = panel(); pnl.margin=15; pnl.pack(2,4);

clr=[82 82 82]./255;

cLabel = {'Confusion'};

% plot the discriminability of touch from baseline
pnl(1,1).select()

Analyze.plotEventRelatedAverage(cvAcc,cLabel,'TimeVec',TimeCenters,...
    'Colors',clr,'useBootStrap');

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

tEvents = [0 4.5];
plt.vline(tEvents,'k-');
ylims = ylim();
ytext = ylims(2) + .05*diff(ylims);
text(tEvents(1),ytext,'Go',StyleArgs{:});

line([ax1.XLim(1) ax1.XLim(end)],[25 25],'Color',[115 115 115]./255,'LineStyle','--');
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