function CorrMatrixAction(Labels,DecodeLabels,OutDir,Tag,varargin)
% Plots confusion matrix (decode analysis) for decoding each (relevant)
% condition.

% [~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh',0.001);
[~,Pmethod]  = Utilities.ProcVarargin(varargin,'Pmethod','fdr');
[~,Phase] = Utilities.ProcVarargin(varargin,'Phase','Go');
[~,reorder]=Utilities.ProcVarargin(varargin,'reorder',[]);
[~,linkagetype]=Utilities.ProcVarargin(varargin,'linkagetype','weighted');
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[varargin,LabelsNames]=Utilities.ProcVarargin(varargin,'LabelsNames',{});
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

%% reorganize label names

newLabelsAbbrev = {{'Felt Pinch','Felt Press','Felt Rub','Felt Tap',...
    'Obs Pinch','Obs Press','Obs Rub','Obs Tap',...
    'Felt Pinch','Felt Press','Felt Rub','Felt Tap',...
    'Obs Pinch','Obs Press','Obs Rub','Obs Tap'},...
    {'Cheek','Shoulder'}};

%% Draw matrix
NC=16;
nSplits = length(unique(cellfun(@(x)x(1),DecodeLabels,'UniformOutput',false)));
plt.fig('units','inches','width',13,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=10; pnl.pack(1,2);
pnl(1,1).pack('v',{3/4,1/4},'h',{1/4,3/4})
pnl(1,2).pack('v',{3/4,1/4},'h',{1/4,3/4})
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1,1,2).select();

rhomod = rho;
rhomod(tril(rhomod,-1)==rhomod)=nan;
imAlpha=ones(size(rhomod));imAlpha(isnan(rhomod))=0;
H = imagesc(rhomod,'AlphaData',imAlpha);
% H = imagesc(rho);

axis image;

newLabelsAbbrev{1};

title('Correlation between conditions','fontsize',9);
set(gca,'XTick',1:NConds,'XTickLabel',newLabelsAbbrev{1});
set(gca,'YTick',1:NConds,'YTickLabel',newLabelsAbbrev{1});
set(gca,'YDir','normal')
colorbar;
colormap(jet);

xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',9);
axis image;
% colormap(cool);
for i=[0:4:NConds]
    
    for i=[0:4:NConds]
    line([i+0.5 i+0.5],[0.5 i+.5],'Color','k');
    line([i+.5 NConds+0.5],[i+0.5 i+0.5],'Color','k');
end
%     if i==0
%         line([i+0.5 i+0.5],[0.5 i+1.5],'Color','k');
%     else
%         line([i+0.5 i+0.5],[0.5 i-4+1.5],'Color','k');
%     end
    line([i+.5 NConds+0.5],[i+0.5 i+0.5],'Color','k');
%     if i~=max(NConds)
%         line([i+4.5 i+4.5],[i+1.5 i+4.5],'Color','k', 'LineStyle','--');
%         line([i+1.5 i+4.5],[i+1.5 i+1.5],'Color','k','LineStyle','--');
%     end
    
    %     line([0.5 NConds+0.5],[i+0.5 i+.5],'Color','k');
    %     line([i+.5 i+0.5],[0.5 NConds+0.5],'Color','k');
    %     if i~=max(NConds)
    %         line([i+1.5 i+4.5],[i+1.5 i+1.5],'Color','k','LineStyle','--');
    %         line([i+1.5 i+1.5],[i+1.5 i+4.5],'Color','k','LineStyle','--');
    %     end
    
end

axis tight
ax1 = gca;
set(ax1,'box','off','FontWeight','bold','FontSize',9);
ax1.XRuler.Axle.LineStyle = 'none';
set(gca,'TickDir','out');
set(ax1,'TickLength',[0.00 0.0]);

CondNames = newLabelsAbbrev{2};

count=1;
for i=[NConds/(nSplits*2):NConds/nSplits:NConds]
    text(i+0.5,-5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255);
    text(-5,i+0.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255,'Rotation',90);
    count = count+1;
end

%%

% %reorder based on body part
% 
% reorder = [1 9 2 10 3 11 4 12 5 13 6 14 7 15 8 16];
% 
% for i=1:length(TestSet)
%     
%     d1=TestSet{i}(:,idx2Incl);
%     d2=TrainSet{i}(:,idx2Incl);
%     
%     d1 = d1(:,reorder);
%     d2 = d2(:,reorder);
%     
%     rhoRep(:,:,i) = corr(d1,d2);
%     
% end
% 
% rho=mean(rhoRep,3);
% 
% newLabelsNames = {{'Felt Pinch','Obs Pinch','Felt Press','Obs Press',...
%     'Felt Rub','Obs Rub','Felt Tap','Obs Tap',...
%     'Felt Pinch','Obs Pinch','Felt Press','Obs Press',...
%     'Felt Rub','Obs Rub','Felt Tap','Obs Tap'},...
%     {'Cheek','Shoulder'}};
% 
% 
% %% Draw matrix
% 
% pnl(1,2,1,2).select();
% 
% rhomod = rho;
% rhomod(tril(rhomod,-1)==rhomod)=nan;
% imAlpha=ones(size(rhomod));imAlpha(isnan(rhomod))=0;
% H = imagesc(rhomod,'AlphaData',imAlpha);
% 
% axis image;
% 
% newLabelsNames{1};
% 
% title('Correlation between conditions','fontsize',9);
% set(gca,'XTick',1:NConds,'XTickLabel',newLabelsNames{1});
% set(gca,'YTick',1:NConds,'YTickLabel',newLabelsNames{1});
% set(gca,'YDir','normal')
% colorbar;
% colormap(jet);
% 
% xlim([.25 NC+.25]); xtickangle(45); set(gca, 'fontweight','bold','fontsize',9);
% axis image;
% % colormap(cool);
% % for i=[0:4:NConds]
% %     line([i+0.5 i+0.5],[0.5 i+.5],'Color','k');
% %     line([i+.5 NConds+0.5],[i+0.5 i+0.5],'Color','k');
% % end
% 
% axis tight
% ax1 = gca;
% set(ax1,'box','off','FontWeight','bold','FontSize',9);
% ax1.XRuler.Axle.LineStyle = 'none';
% set(gca,'TickDir','out');
% set(ax1,'TickLength',[0.00 0.0]);
% 
% CondNames = newLabelsNames{2};
% 
% count=1;
% for i=[NConds/(nSplits*2):NConds/nSplits:NConds]
%     text(i+0.5,-5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
%         'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255);
%     text(-5,i+0.5,CondNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
%         'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255,'Rotation',90);
%     count = count+1;
% end

%%
filename = 'CorrMat';
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');
% plt.SaveFigure(saveFig,FigsDir,[filename],'SVGI');


%% Also compute for pairs differing by only person, bp, and action

LabelsAbbrev = {'NPiC','NPrC','NRbC','NTpC',...
    'TPiC','TPrC','TRbC','TTpC',...
    'NPiS','NPrS','NRbS','NTpS',...
    'TPiS','TPrS','TRbS','TTpS'};

reorder = [1:NConds];

for i=1:length(TestSet)
    
    d1=TestSet{i}(:,idx2Incl);
    d2=TrainSet{i}(:,idx2Incl);
    
    d1 = d1(:,reorder);
    d2 = d2(:,reorder);
    
    rhoRep(:,:,i) = corr(d1,d2);
    
end

rho=mean(rhoRep,3);

corrVals = {[], [], []};
for i = 1:NConds
    for j = (i+1):NConds
        c1 = LabelsAbbrev{i};
        c2 = LabelsAbbrev{j};
        [n,diffId] = numdiffFn(c1,c2);
        if n ~= 1
            continue;
        end
        ii = find(diffId);
        if length(ii) > 1
            error('??');
        end
        corrVals{ii} = [corrVals{ii} rho(i,j)];
    end
end

%% Plot
DiffVarNames = {'Person','Action','BodyPart'};

reorder = [3 2 1];
corrVals = corrVals(reorder);
DiffVarNames = DiffVarNames(reorder);

% StyleArgs = {'FontName','Arial','FontSize',14};
% StyleArgsSmall = {'FontName','Arial','FontSize',12};
NC=16;
plt.fig('units','inches','width',12,'height',7,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(1,2); pnl(1,1).pack('v',{3/4,1/4},'h',{1/4,3/4})
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1,1,2).select();

for i = 1:length(corrVals)
    [mu(i),lb(i),ub(i)] = MeanAndCI(corrVals{i});
end

data = {mu'; [lb; ub]};
groupidx = 1:3;
clr = lines(7);
h  = plt.barpatchGroup(data, 'groupidx', groupidx, ...
    'fontsize',6,         ...
    'groupspace', 2,      ...
    'patchbar', 0,        ...
    'barname', DiffVarNames,      ...
    'barcmap',clr([3 4 2],:));
set(h.leg,'FontSize',9,'location', 'southoutside');

xlim([.25 NC+.25]); xtickangle(45);
set(gca, 'fontweight','bold','FontName','Arial','fontsize',9');
set(h.leg,'FontSize',9, 'FontName','Arial','FontWeight', 'bold');
set(h.leg,'String',{'Same action different person and/or body part','Same body part, different person and/or action','Same person, different action and/or body part'})
set(h.xline,'Visible','off');
ax1=gca;
xlim([.25 NC+.25]); xtickangle(45);
set(gca, 'fontweight','bold','FontName','Arial','fontsize',9');
set(ax1,'box','off','FontWeight','bold','FontSize',9);
ax1.XRuler.Axle.LineStyle = 'none';
set(gca,'TickDir','out');
set(ax1,'TickLength',[0.00 0.0]);
ylim([0.1 .9])
set(gca,'XTick',1:2,'XTickLabel',{'Action','Person','Body Part'});
ylabel('Correlation Coefficient');
xlabel('Differing Condition');
title('Correlation by Condition')
set(gca, 'fontweight','bold','FontName','Arial','fontsize',9');


filename = 'CorrelationByCondition';
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');
% plt.SaveFigure(saveFig,FigsDir,[filename],'SVGI');

%%

% Z = linkage(MEAS(:,reorder)',linkagetype,'correlation');
Z = linkage(MEAS','complete', 'correlation');

%%
plt.fig('units','inches','width',12,'height',8,'font','Arial','fontsize',14);
pnl2 = panel(); pnl2.margin=25;
pnl2.pack(2,3);

% currLabels = {'Fc Pi','Fc Pr','Fc Ru','Fc Ta',...
%     'Fs Pi','Fs Pr','Fs Ru','Fs Ta',...
%     'Oc Pi','Oc Pr','Oc Ru','Oc Ta',...
%     'Os Pi','Os Pr','Os Ru','Os Ta'};

currLabels = {'Felt Pinch Cheek','Felt Press Cheek','Felt Rub Cheek','Felt Tap Cheek',...
    'Felt Pinch Shoulder','Felt Press Shoulder','Felt Rub Shoulder','Felt Tap Shoulder',...
    'Obs Pinch Cheek','Obs Press Cheek','Obs Rub Cheek','Obs Tap Cheek',...
    'Obs Pinch Shoulder','Obs Press Shoulder','Obs Rub Shoulder','Obs Tap Shoulder'};

pnl2(1,1).select();
H = dendrogram(Z,'Labels',currLabels,'ColorThreshold',0.70);

set(H,'LineWidth',2);
yticklabels = strsplit(sprintf('%.2f,',1-get(gca,'ytick')),',');
set(gca,'yticklabels',yticklabels(1:end),StyleArgs{:});
ylabel('Correlation coefficient (r)',StyleArgs{:});
title(['Hierarchical Clustering'],StyleArgs{:});
xtickangle(45);

%% save figure
filename = 'Clustering';
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');
% plt.SaveFigure(saveFig,FigsDir,[filename],'SVGI');

end

function [n,diffId] = numdiffFn(a1,a2)
% How many conditions different and which conditions (1st 2nd 3rd) diff.
diffId = ~[strcmp(a1(1),a2(1)) strcmp(a1(2:3),a2(2:3)) strcmp(a1(4),a2(4))];
n = sum(diffId);

end