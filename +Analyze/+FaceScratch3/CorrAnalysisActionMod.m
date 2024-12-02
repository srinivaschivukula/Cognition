function [ output_args ] = CorrAnalysisActionMod(TaskLabel,varargin )
%SPECIFICITYANALYSIS Summary of this function goes here
%   Detailed explanation goes here

% [varargin,PThresh]=Utilities.ProcVarargin(varargin,'PThresh',.001);
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');
[varargin,Labels]=Utilities.ProcVarargin(varargin,'Labels',{});
[varargin,Labels2Incl]=Utilities.ProcVarargin(varargin,'Labels2Incl',{});
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[varargin,Dates]=Utilities.ProcVarargin(varargin,'Dates',[]);
[varargin,trialdata]   = Utilities.ProcVarargin(varargin,'AllTrialData',{});

datafile = fullfile(env.get('results'),'FaceScratch3','SUAnal',[TaskLabel '-Go'],'PopData','p1-AllDays.mat');
load(datafile); % Loads PlotData

if isempty(Labels) | isempty(Labels2Incl)
    idx2incl = 1:length(Labels);
else
    idx2incl = find(ismember(Labels,Labels2Incl));
end

if ~isempty(Dates)
    if ~iscell(Dates)
        Dates = {Dates};
    end
    PlotData = Analyze.SubSelectTrials(PlotData,'TaskDate',Dates);
end


ASF = cell2mat(vertcat(Analyze.returnFieldValues(PlotData,'Unit')));
PVals = horzcat(PlotData{1}.Pvals{:})';
PVals = PVals(:,idx2incl);
[IsSig,effa] = Utilities.MultipleComparisonsCorrection(PVals,'method','fdr');
effa

Coef = horzcat(PlotData{1}.Coef{:})';
Coef = Coef(:,idx2incl);
CoefCI = horzcat(PlotData{1}.CoefCI{:})';
CoefCI = CoefCI(:,idx2incl);
tmp = diff(CoefCI);
CoefCI = tmp(1:2:end,:);

NConds = size(PVals,2);
MEAS = Coef./CoefCI;

GoodIdx = any(IsSig,2);
MEAS = MEAS(GoodIdx,:);
FigsDir = fullfile(env.get('results'),'FaceScratch3','SUAnal',[TaskLabel '-Go'],'PopData','Mine');


%% Compute correlation matrix
rho = corr(MEAS);
plt.fig('units','inches','width',6,'height',5,'font','Arial','fontsize',16);
pnl2 = panel(); pnl2.margin=20; pnl2.pack(1,1);
imagesc(rho);
axis image;
title('Correlation between conditions');
set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev);
set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev);
set(gca,'YDir','normal')
colorbar;
colormap(jet);
filename = 'CorrMat';
if ~isempty(Dates)
    filename = [filename '-' strjoin(Dates,'-')];
end
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');


if ~strcmp(TaskLabel,'XTouchRF') && ~strcmp(TaskLabel,'XActionType')
    return;
end


if strcmp(TaskLabel,'XTouchRF')
    %%
    DiffVarNames = {'Side','BodyPart'};
    pairs{1} = [3 4; 5 6; 7 8; 9 10];
    pairs{2} = [3 5; 3 7; 3 9; 5 7; 5 9; 7 9; 4 6; 4 8; 4 10; 6 8; 6 10; 8 10];
    for i = 1:length(DiffVarNames)
        corrVals{i} = [];
        for j = 1:length(pairs{i})
            corrVals{i} = [corrVals{i} rho(pairs{i}(j,1),pairs{i}(j,2))];
        end
    end
    
    plt.fig('units','inches','width',3,'height',5,'font','Arial','fontsize',16);
    pnl2 = panel(); pnl2.margin=20; pnl2.pack(1,1);
    
    for i = 1:length(corrVals)
        [mu(i),lb(i),ub(i)] = MeanAndCI(corrVals{i});
    end
    
    data = {mu'; [lb; ub]};
    groupidx = 1:2;
    xl = 'Differing Condition';
    yl = 'Correlation'
    t = 'Correlation by Condition';
    clr = lines(2);
    h  = plt.barpatchGroup(data, 'groupidx', groupidx, ...
        'xl', xl,             ...
        'yl', yl,             ...
        't', t,               ...
        'groupspace', 2,      ...
        'patchbar', 0,        ...
        'barname', DiffVarNames,      ...
        'barcmap',clr);
    set(h.leg,'FontSize',14);
    plt.SaveFigure(saveFig,FigsDir,['CorrByVarDiff'],'PNG','SVGI');

    %%
    return;
end
%% Also compute for pairs differing by only person, bp, and action
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

reorder = [2 1 3];
corrVals = corrVals(reorder);
DiffVarNames = DiffVarNames(reorder);

% StyleArgs = {'FontName','Arial','FontSize',14};
% StyleArgsSmall = {'FontName','Arial','FontSize',12};
plt.fig('units','inches','width',5,'height',5,'font','Arial','fontsize',16);
pnl2 = panel(); pnl2.margin=20; pnl2.pack(1,1);

for i = 1:length(corrVals)
    [mu(i),lb(i),ub(i)] = MeanAndCI(corrVals{i});
end

data = {mu'; [lb; ub]};
groupidx = 1:3;
xl = 'Differing Condition';
yl = 'Correlation'
t = 'Correlation by Condition';
clr = lines(3);
h  = plt.barpatchGroup(data, 'groupidx', groupidx, ...
        'xl', xl,             ...
        'yl', yl,             ...
        't', t,               ...
        'groupspace', 2,      ...
        'patchbar', 0,        ...
        'barname', DiffVarNames,      ...
        'barcmap',clr);
set(h.leg,'FontSize',14);
plt.SaveFigure(saveFig,FigsDir,['CorrByVarDiff'],'PNG','SVGI');


%% Dendrogram/clustering
% Z = linkage(MEAS','single','correlation');
Z = linkage(MEAS','complete','correlation'); % yes
% Z = linkage(MEAS','average','correlation'); % no
Z = linkage(MEAS','weighted','correlation'); % yes
% Z = linkage(MEAS','centroid','correlation'); % NO
% Z = linkage(MEAS','median','correlation'); % NO
plt.fig('units','inches','width',6,'height',5,'font','Arial','fontsize',13);
% dendrogram(Z,'Labels',LabelsAbbrev,'Reorder',[1 2 3 4 9 10 12 11 5 7 6 8 13 14 15 16]);
dendrogram(Z,'Labels',LabelsAbbrev);
% reordering to look better

% Convert y labels into correlation measures
yticklabels = strsplit(sprintf('%.2f,',1-get(gca,'ytick')),',');
set(gca,'yticklabels',yticklabels(1:end));
ylabel('Correlation coefficient (r)');
title(sprintf('Hierarchical clustering'));
xtickangle(45);

%
plt.SaveFigure(saveFig,FigsDir,['Dendrogram'],'PNG','SVGI');

%%
Z = linkage(MEAS','weighted','correlation');
plt.fig('units','inches','width',6,'height',5,'font','Arial','fontsize',13);
dendrogram(Z,'Labels',LabelsAbbrev);
% reordering to look better

% Convert y labels into correlation measures
yticklabels = strsplit(sprintf('%.2f,',1-get(gca,'ytick')),',');
set(gca,'yticklabels',yticklabels(1:end));
ylabel('Correlation coefficient (r)');
title(sprintf('Hierarchical clustering'));
xtickangle(45);

plt.SaveFigure(saveFig,FigsDir,['Dendrogram-Weighted'],'PNG','SVGI');


%% Check number of units in common between cheek, neck, and shoulder
% if strcmp(TaskLabel,'XTouchRF')
%     figure;
%     Cheek = any(IsSig(:,[5 6]),2);
%     Shoulder = any(IsSig(:,[7 8]),2);
%     Neck = any(IsSig(:,[9 10]),2);
% 
%     A = [nnz(Cheek) nnz(Shoulder) nnz(Neck)];
%     I = [nnz(Cheek & Shoulder) nnz(Cheek & Neck) nnz(Shoulder & Neck) nnz(Cheek & Shoulder & Neck)];
%     venn(A,I,'FaceAlpha',0.3);
%     legend('Cheek','Shoulder','Neck');
%     title({sprintf('%d %d %d',A),sprintf('%d %d %d %d',I)});
%     plt.SaveFigure(saveFig,FigsDir,['ChShNkVenn'],'PNG','SVGI');
% 
% end


%% Number of units invariant for different combinations of N bodyparts
%% N = 1 2 3 ...
% if strcmp(TaskLabel,'XTouchRF')
%     figure;
%     
%     
%     Cheek = any(IsSig(:,[5 6]),2);
%     Shoulder = any(IsSig(:,[7 8]),2);
%     Neck = any(IsSig(:,[9 10]),2);
% 
%     A = [nnz(Cheek) nnz(Shoulder) nnz(Neck)];
%     I = [nnz(Cheek & Shoulder) nnz(Cheek & Neck) nnz(Shoulder & Neck) nnz(Cheek & Shoulder & Neck)];
%     plt.venn(A,I,'FaceAlpha',0.3);
%     legend('Cheek','Shoulder','Neck');
%     title({sprintf('%d %d %d',A),sprintf('%d %d %d %d',I)});
%     plt.SaveFigure(saveFig,FigsDir,['InvarForNBP'],'PNG','SVGI');
% 
% end
end

function [n,diffId] = numdiffFn(a1,a2)
% How many conditions different and which conditions (1st 2nd 3rd) diff.
diffId = ~[strcmp(a1(1),a2(1)) strcmp(a1(2:3),a2(2:3)) strcmp(a1(4),a2(4))];
n = sum(diffId);

end
