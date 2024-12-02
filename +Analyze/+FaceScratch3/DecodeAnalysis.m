function DecodeAnalysis(AllTrialData,Labels,DecodeLabels,TimeWindow,OutDir,Tag,varargin)
% Plots confusion matrix (decode analysis) for decoding each (relevant)
% condition.

% [~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh',0.001);
[~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh','fdr');
[~,NIter] = Utilities.ProcVarargin(varargin,'NIter',100);
[~,Phase] = Utilities.ProcVarargin(varargin,'Phase','Go');
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[~,reorder]=Utilities.ProcVarargin(varargin,'reorder',[]);
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');

StyleArgs = {'FontSize',14,'FontName','Arial'};

NConds = length(DecodeLabels);

%% Only include units with tunign to at least one of the relevant variables
% TODO: Dropping curve version.
ASF = Analyze.FaceScratch3.FindSignificantUnits(Tag,'PThresh',PThresh,'Labels',Labels,'Labels2Incl',DecodeLabels);
% ASF = Analyze.FaceScratch3.FindSignificantUnits(Tag,'PThresh',PThresh,'Labels',Labels,'Labels2Incl',DecodeLabels,'Mirror','AllTrialData',AllTrialData);
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

%% Rearrange Confusion matrix label ordering to match DecodeLabels
% ConfMatOld = ConfMat;
% for d = 1:length(UDates)
%     for
%     ConfMat{d};
% end
%
% for k = 1:NIter
%     if mod(k,10) == 0
%         fprintf('.');
%     end
%     if mod(k,100)==0
%         fprintf('\n');
%     end
%     [testFR,testLabels,remFR,~,remLabels] = Analyze.SampleFeaturePopulation(FR,LabelId,'Type','RandSample','NumObservationsPerCondition',1);
%     [trainFR,trainLabels] = Analyze.SampleFeaturePopulation(remFR,remLabels,'Type','RandSample','NumObservationsPerCondition',10);
%     mdl = fitcdiscr(trainFR,trainLabels,'discrimType','diagLinear');
%     predLabels = predict(mdl,testFR);
%     ind = sub2ind(size(ConfMat),testLabels,predLabels);
%     ConfMat(ind) = ConfMat(ind)+1;
% %     for i = 1:length(testLabels)
% %         ConfMat(testLabels(i),predLabels(i)) = ConfMat(testLabels(i),predLabels(i)) + 1;
% %     end
% end
%% Draw per day
% plt.fig('units','inches','width',20,'height',5,'font','Arial','fontsize',14);
% pnl2 = panel(); pnl2.margin=25;
% pnl2.pack(1,length(UDates));
%
% for d = 1:length(UDates)
%     cMat = ConfMat{d}*100;
%
%     pnl2(1,d).select();
%     imagesc(cMat);
%     caxis([0 100]);
%     colorbar;
%     axis([0 NConds 0 NConds]+.5);
%     %     set(gca);
%     set(gca,'XTick',1:NConds,'XTickLabel',DecodeLabels,StyleArgs{:},'XTickLabelRotation',30);
%     %     set(gca);
%     set(gca,'YTick',1:NConds,'YTickLabel',DecodeLabels,StyleArgs{:});
%     title(UDates{d});
%     axis square;
%     colormap(jet);
% end
% suptitle(sprintf('%s Confusion Matrix',Tag));
% plt.SaveFigure(2,fullfile(OutDir,'PopData','Models'),['Model-ConfusionMatrixPerDay'],'PNG','SVGI')
%% Draw matrix
NC=16;
nSplits = length(unique(cellfun(@(x)x(1),DecodeLabels,'UniformOutput',false)));
plt.fig('units','inches','width',14,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(1,3); pnl(1,1).pack('v',{3/4,1/4},'h',{1/4,3/4})
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1,1,2).select();

ConfMatMu = zeros(NConds,NConds);
for d = 1:length(UDates)
    ConfMatMu = ConfMatMu + ConfMat{d};
end
ConfMatMu = ConfMatMu/length(UDates)*100;

H = imagesc(ConfMatMu);

axis image;

if ~isempty(reorder)
    Labels2Use = LabelsAbbrev{1};
else
    Labels2Use = LabelsAbbrev{2};
end

title('Confusion Matrix','fontsize',9);
set(gca,'XTick',1:NConds,'XTickLabel',Labels2Use);
set(gca,'YTick',1:NConds,'YTickLabel',Labels2Use);
set(gca,'YDir','normal')
colorbar;
colormap(jet);

xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',9);
axis image;
colormap(cool);

vline(0.5:nSplits:NConds+0.5,'k')
hline(0.5:nSplits:NConds+0.5,'k')

axis tight
ax1 = gca;
set(ax1,'box','off','FontWeight','bold','FontSize',9);
ax1.XRuler.Axle.LineStyle = 'none';
set(gca,'TickDir','out');
set(ax1,'TickLength',[0.00 0.0]);

if ~isempty(reorder)
    CondNames = LabelsAbbrev{2};
else
    CondNames = LabelsAbbrev{1};
end

count=1;
for i=[nSplits/2:nSplits:NConds]
    text(i+0.5,-1.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255);
    text(-1.5,i+0.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255,'Rotation',90);
    count = count+1;
end

%%
FigsDir = fullfile(env.get('results'),'FaceScratch3','SUAnal',[Tag '-Go'],'PopData','Mine');
filename = 'ConfusionMatrix';
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');
% plt.SaveFigure(saveFig,FigsDir,[filename],'SVGI');
