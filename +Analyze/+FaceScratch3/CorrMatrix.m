function CorrMatrix(Labels,DecodeLabels,OutDir,Tag,varargin)
% Plots confusion matrix (decode analysis) for decoding each (relevant)
% condition.

% [~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh',0.001);
[~,Pmethod]  = Utilities.ProcVarargin(varargin,'Pmethod','fdr');
[~,Phase] = Utilities.ProcVarargin(varargin,'Phase','Go');
[~,reorder]=Utilities.ProcVarargin(varargin,'reorder',[]);
[~,linkagetype]=Utilities.ProcVarargin(varargin,'linkagetype','weighted');
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');
[varargin,trialdata]   = Utilities.ProcVarargin(varargin,'AllTrialData',{});
[varargin,PhaseStart]   = Utilities.ProcVarargin(varargin,'PhaseStart',[]);
[varargin,PhaseDur]   = Utilities.ProcVarargin(varargin,'PhaseDur',[]);
[varargin,BaselineWindow]   = Utilities.ProcVarargin(varargin,'BaselineWindow',[]);
[varargin,Phase]   = Utilities.ProcVarargin(varargin,'Phase','Go');


StyleArgs = {'FontSize',9,'FontName','Arial','FontWeight','bold'};

NConds = length(DecodeLabels);

%% Only include units with tunign to at least one of the relevant variables
datafile = fullfile(env.get('results'),'FaceScratch3','SUAnal',[Tag '-Go'],'PopData','p1-AllDays.mat');
load(datafile); % Loads PlotData
PVals = horzcat(PlotData{1}.Pvals{:})';
[H,PThresh] = Utilities.MultipleComparisonsCorrection(PVals,'method', Pmethod);
IsSig = PVals < PThresh;
ASF = cell2mat(vertcat(Analyze.returnFieldValues(PlotData,'Unit')));

if isempty(Labels) | isempty(DecodeLabels)
    idx2Incl = 1:length(Labels);
else
    idx2Incl = find(ismember(Labels,DecodeLabels));
end
GoodIdx = any(IsSig(:,idx2Incl),2);
ASF = ASF(GoodIdx,:);

Coef = horzcat(PlotData{1}.Coef{:})';
Coef = Coef(GoodIdx,idx2Incl);
CoefCI = horzcat(PlotData{1}.CoefCI{:});
CoefCI = diff(CoefCI,[],2);
CoefCI = CoefCI(:,1:2:end)';
CoefCI = CoefCI(GoodIdx,idx2Incl);

MEAS = Coef./CoefCI;

if isempty(reorder)
    reorder = 1:NConds;
end
% reorder = [1 3 2 4];

FigsDir = fullfile(env.get('results'),'FaceScratch3','SUAnal',[Tag '-Go'],'PopData','Mine');


%%
%% reorganize splits
for i=1:size(PlotData{1}.testSplit{1},1)
    for j=1:length(PlotData{1}.testSplit)
        TestSet{i}(j,:)=PlotData{1}.testSplit{j}(i,:);
        TrainSet{i}(j,:)=PlotData{1}.trainSplit{j}(i,:);
    end
end
%%

for i=1:length(TestSet)
    
    d1=TestSet{i}(:,idx2Incl);
    d2=TrainSet{i}(:,idx2Incl);
    
    d1 = d1(:,reorder);
    d2 = d2(:,reorder);
    
    rhoRep(:,:,i) = corr(d1,d2);
    
end

rho=mean(rhoRep,3);

%% Draw matrix
NC=16;
nSplits = length(unique(cellfun(@(x)x(1),DecodeLabels,'UniformOutput',false)));
plt.fig('units','inches','width',14,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(1,3); pnl(1,1).pack('v',{3/4,1/4},'h',{1/4,3/4})
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1,1,2).select();

rhomod = rho;
rhomod(tril(rhomod,-1)==rhomod)=nan;
imAlpha=ones(size(rhomod));imAlpha(isnan(rhomod))=0;
H = imagesc(rhomod,'AlphaData',imAlpha);

axis image;



title('Correlation between conditions','fontsize',9);
set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev{1});
set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev{1});
set(gca,'YDir','normal')
colorbar;
colormap(jet);

xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',9);
axis image;
colormap(cool);
for i=[0:nSplits:NConds]
    line([i+0.5 i+0.5],[0.5 i+.5],'Color','k');
    line([i+.5 NConds+0.5],[i+0.5 i+0.5],'Color','k');
end

axis tight
ax1 = gca;
set(ax1,'box','off','FontWeight','bold','FontSize',9);
ax1.XRuler.Axle.LineStyle = 'none';
set(gca,'TickDir','out');
set(ax1,'TickLength',[0.00 0.0]);

CondNames = LabelsAbbrev{2};

count=1;
for i=[nSplits/2:nSplits:NConds]
    text(i+0.5,-1.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255);
    text(-1.5,i+0.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255,'Rotation',90);
    count = count+1;
end

%%
filename = 'CorrMat';
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');
% plt.SaveFigure(saveFig,FigsDir,[filename],'SVGI');

%% Compute distance matrices

UDates = Analyze.returnUniqueFieldValues(trialdata,'Date');

for i=1:length(Labels)
    for d=1:length(UDates)
        cLabel = Labels{i};
        ctrialdata = Analyze.SubSelectTrials(trialdata, 'Date', UDates{d}, 'Phase', Phase, 'Condition', cLabel);
        btrialdata = Analyze.SubSelectTrials(trialdata, 'Date', UDates{d}, 'Phase', 'Delay', 'Condition', cLabel);
        cFR = squeeze(Analyze.getNeuralData(ctrialdata, ASF, [PhaseStart PhaseStart+PhaseDur]));
        cFR(:,any(isnan(cFR)))=[];
        FR{i,d} = cFR;
        cbFR = squeeze(Analyze.getNeuralData(btrialdata, ASF, BaselineWindow));
        cbFR(:,any(isnan(cbFR)))=[];
        bFR{i,d} = cbFR;
    end
end

for i=1:size(FR,1)
    tmpFR = cat(2,FR{i,:});
    tmpBase = cat(2,bFR{i,:});
    newFR{i} = tmpFR-tmpBase;
end

newFR = newFR(idx2Incl);
newFR = newFR(reorder);

[dEuc,dMah,dMahDiag,dCor] = Analyze.RFHandEye.ComputeDistances(newFR);

%% Plot Distance Measures

NC =16;
nSplits = length(unique(cellfun(@(x)x(1),DecodeLabels,'UniformOutput',false)));

cmap = bone;
% cmap = cbrewer('div','RdBu',150);

plt.fig('units','inches','width',16,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,4);
pnl.fontsize=11; pnl.fontname='Arial';

pnl(1,1).select();
h1 = imagesc(dEuc);
title('Euclidean','fontsize',6);
set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev{1});
set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev{1});
set(gca,'YDir','normal')
colorbar; colormap(cmap);
xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',9');
axis image;
vline([0.5:nSplits:NConds+0.5],'k-');
hline([0.5:nSplits:NConds+0.5],'k-');
CondNames = LabelsAbbrev{2};
count=1;
for i=[nSplits/2:nSplits:NConds]
    text(i+0.5,-1.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255);
    text(-1.5,i+0.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255,'Rotation',90);
    count = count+1;
end

pnl(1,2).select();
h2 = imagesc(dMah);
title('Mahalanobis','fontsize',6);
set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev{1});
set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev{1});
set(gca,'YDir','normal')
colorbar; colormap(cmap);
xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',9');
axis image;
vline([0.5:nSplits:NConds+0.5],'k-');
hline([0.5:nSplits:NConds+0.5],'k-');
CondNames = LabelsAbbrev{2};
count=1;
for i=[nSplits/2:nSplits:NConds]
    text(i+0.5,-1.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255);
    text(-1.5,i+0.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255,'Rotation',90);
    count = count+1;
end

pnl(1,3).select();
h3 = imagesc(dMahDiag);
title('Mahalanobis Diag','fontsize',6);
set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev{1});
set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev{1});
set(gca,'YDir','normal')
colorbar; colormap(cmap);
xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',9');
axis image;
vline([0.5:nSplits:NConds+0.5],'k-');
hline([0.5:nSplits:NConds+0.5],'k-');
CondNames = LabelsAbbrev{2};
count=1;
for i=[nSplits/2:nSplits:NConds]
    text(i+0.5,-1.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255);
    text(-1.5,i+0.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255,'Rotation',90);
    count = count+1;
end

pnl(1,4).select();
h3 = imagesc(dCor);
title('Corr','fontsize',6);
set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev{1});
set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev{1});
set(gca,'YDir','normal')
colorbar; colormap(cmap);
xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',9');
axis image;
vline([0.5:nSplits:NConds+0.5],'k-');
hline([0.5:nSplits:NConds+0.5],'k-');CondNames = LabelsAbbrev{2};
count=1;
for i=[nSplits/2:nSplits:NConds]
    text(i+0.5,-1.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255);
    text(-1.5,i+0.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255,'Rotation',90);
    count = count+1;
end


%%
filename = 'Distance';
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');
% plt.SaveFigure(saveFig,FigsDir,[filename],'SVGI');


%%

Z = linkage(MEAS',linkagetype,'correlation');

%%
plt.fig('units','inches','width',12,'height',8,'font','Arial','fontsize',14);
pnl2 = panel(); pnl2.margin=25;
pnl2.pack(2,3);

currLabels = {'Felt Cheek','Felt Shoulder','Observed Cheek','Observed Shoulder'};

pnl2(1,1).select();
H = dendrogram(Z,'Labels',currLabels,'ColorThreshold',0.8);

set(H,'LineWidth',2);
yticklabels = strsplit(sprintf('%.2f,',1-get(gca,'ytick')),',');
set(gca,'yticklabels',yticklabels(1:end),StyleArgs{:});
ylabel('Correlation coefficient (r)',StyleArgs{:});
title([Tag ' Hierarchical Clustering (' linkagetype ')'],StyleArgs{:});
xtickangle(45);

%% save figure
filename = 'Clustering';
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');
% plt.SaveFigure(saveFig,FigsDir,[filename],'SVGI');

%% computing the correlation by condition (person, body part)

clear corrVals;
DiffVarNames = {'Person','Body Part'};
pairs{1} = [1 3;2 4];
pairs{2} = [1 2;3 4];
for i = 1:length(DiffVarNames)
    corrVals{i} = [];
    for j = 1:length(pairs{i})
        corrVals{i} = [corrVals{i} rhoRep(pairs{i}(j,1),pairs{i}(j,2),:)];
    end
    corrVals{i} = squeeze(corrVals{i});
    corrVals{i} = corrVals{i}(:);
end

%% for plotting

mu = cellfun(@mean,corrVals);
ci = [bootci(2000,{@mean, corrVals{1}},'type','per') bootci(2000,{@mean, corrVals{2}},'type','per')];
data = {mu; [ci]};
groupidx = 1:2;

clr = lines(7);
%%
NC=16;
plt.fig('units','inches','width',12,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(1,2); pnl(1,1).pack('v',{3/4,1/4},'h',{1/4,3/4})
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1,1,2).select();

yl = 'Correlation';
t = 'Correlation by Condition';
clr = lines(7);
h  = plt.barpatchGroup(data, 'groupidx', groupidx, ...
    'yl', yl,             ...
    'fontsize',6,         ...
    'groupspace', 2,      ...
    'patchbar', 0,        ...
    'barname', DiffVarNames,      ...
    'barcmap',clr([4 2],:));

ax1=gca;

xlim([.25 NC+.25]); xtickangle(45); 
set(gca, 'fontweight','bold','FontName','Arial','fontsize',9');
set(h.leg,'FontSize',9, 'FontName','Arial','FontWeight', 'bold');
set(h.leg,'String',{'Same person, different body part','Same body part, different person'})
set(h.xline,'Visible','off');
ax1=gca;
xlim([.25 NC+.25]); xtickangle(45); 
set(gca, 'fontweight','bold','FontName','Arial','fontsize',9');
set(ax1,'box','off','FontWeight','bold','FontSize',9);
ax1.XRuler.Axle.LineStyle = 'none';
set(gca,'TickDir','out');
set(ax1,'TickLength',[0.00 0.0]);
ylim([0.1 .9])
set(gca,'XTick',1:2,'XTickLabel',{'Person','Body Part'});
ylabel('Correlation Coefficient');
set(gca, 'fontweight','bold','FontName','Arial','fontsize',9');

%%
filename = 'CorrelationByCondition';
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');
% plt.SaveFigure(saveFig,FigsDir,[filename],'SVGI');


end