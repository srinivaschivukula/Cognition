function [ output_args ] = SpecificityAnalysis(TaskLabel,varargin )
%SPECIFICITYANALYSIS Summary of this function goes here
%   Detailed explanation goes here

% [varargin,PThresh]=Utilities.ProcVarargin(varargin,'PThresh',.001);
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');
[varargin,Labels]=Utilities.ProcVarargin(varargin,'Labels',{});
[varargin,Labels2Incl]=Utilities.ProcVarargin(varargin,'Labels2Incl',{});
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[varargin,Dates]=Utilities.ProcVarargin(varargin,'Dates',[]);
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig',1);
[varargin,AllTrialData]   = Utilities.ProcVarargin(varargin,'AllTrialData', {});

datafile = fullfile(env.get('results'),'FaceScratch3','SUAnal',[TaskLabel '-Go'],'PopData','p1-AllDays.mat');
load(datafile); % Loads PlotData

if isempty(Labels) | isempty(Labels2Incl)
    idx2incl = 1:length(Labels);
else
    idx2incl = find(ismember(Labels,Labels2Incl));
end

ASF = cell2mat(vertcat(Analyze.returnFieldValues(PlotData,'Unit')));
mASF = Analyze.FaceScratch3.FindSignificantUnits(TaskLabel,'PThresh','fdr','Labels',Labels,'Labels2Incl',Labels2Incl,'Mirror','AllTrialData',AllTrialData);

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
% si=@(x,y)(x-y)./(abs((x))+abs((y)));
si=@(x,y)(abs(x)-abs(y))./(abs((x))+abs((y)));
spc=linspace(-1,1,10);
MEAS = Coef./CoefCI;

%% Comptue specificity
try
    units2Use = ismember(ASF,mASF,'rows');
    MEAS = MEAS(units2Use,:);
    IsSig = IsSig(units2Use,:);
    for i = 1:NConds
        for j = (i+1):NConds
            idx = any(IsSig(:,[i j]),2);
            x = MEAS(idx,i);
            y = MEAS(idx,j);
            D{i,j} = si(x,y);
        end
    end
catch
    for i = 1:NConds
        for j = (i+1):NConds
            idx = any(IsSig(:,[i j]),2);
            x = MEAS(idx,i);
            y = MEAS(idx,j);
            D{i,j} = si(x,y);
        end
    end
end
%% Plot all specificity plots
% StyleArgs = {'FontName','Arial','FontSize',9,'FontWeight','bold'};
% % StyleArgsSmall = {'FontName','Arial','FontSize',12};
% % plt.fig('units','inches','width',20,'height',12,'font','Arial','fontsize',16);
% % plt.fig('units','inches','width',10,'height',6,'font','Arial','fontsize',16);
% plt.fig('units','inches','width',3*(NConds-1),'height',2*(NConds-1),'font','Arial','fontsize',9);
% pnl2 = panel(); pnl2.margin=20;[10 12 10 10]; pnl2.pack(NConds-1,NConds-1);
% for i = 1:NConds
%     for j = (i+1):NConds
%         leftCond = LabelsAbbrev{j};
%         rightCond = LabelsAbbrev{i};
%         pnl2(i,j-1).select();
%         histogram(D{i,j},spc,'Normalization','pdf');
%
%         [b1,~,stats] = signtest(D{i,j});
%         b1
%         stats
%         m1 = median(D{i,j});
%         title(sprintf('med=%.3f, p=%0.3f',m1,b1),StyleArgs{:});
%         xlabel(sprintf('%s <--> %s',leftCond,rightCond),StyleArgs{:});
%         ax1 = gca;
%         set(ax1,StyleArgs{:});
%         ylim([0 1.1]);
%     end
% end

%% Plot only felt/observed cheek and shoulder
StyleArgs = {'FontName','Arial','FontSize',9,'FontWeight','bold'};
% StyleArgsSmall = {'FontName','Arial','FontSize',12};
% plt.fig('units','inches','width',20,'height',12,'font','Arial','fontsize',16);
% plt.fig('units','inches','width',10,'height',6,'font','Arial','fontsize',16);
plt.fig('units','inches','width',3*(NConds-1),'height',2*(NConds-1),'font','Arial','fontsize',9);
pnl2 = panel(); pnl2.margin=20;[10 12 10 10]; pnl2.pack(NConds-1,NConds-1);

pnl2(1,1).select()
leftCond = LabelsAbbrev{3};
rightCond = LabelsAbbrev{1};
hh1 = histogram(D{1,3},spc,'Normalization','pdf');

[b1,~,stats] = signtest(D{1,3});
b1
stats
m1 = median(D{1,3});
title(sprintf('med=%.3f, p=%0.3f',m1,b1),StyleArgs{:});
xlabel(sprintf('%s <--> %s',leftCond,rightCond),StyleArgs{:});
ax1 = gca;
set(ax1,StyleArgs{:});
ylim([0 1.1]);

pnl2(1,2).select();
leftCond = LabelsAbbrev{4};
rightCond = LabelsAbbrev{2};
hh2 = histogram(D{2,4},spc,'Normalization','pdf');

[b1,~,stats] = signtest(D{2,4});
b1
stats
m1 = median(D{2,4});
title(sprintf('med=%.3f, p=%0.3f',m1,b1),StyleArgs{:});
xlabel(sprintf('%s <--> %s',leftCond,rightCond),StyleArgs{:});
ax1 = gca;
set(ax1,StyleArgs{:});
ylim([0 1.1]);

%%

FigsDir = fullfile(env.get('results'),'FaceScratch3','SUAnal',[TaskLabel '-Go'],'PopData','Mine');

filename = ['Specificity'];
plt.SaveFigure(saveFig,FigsDir,filename,'SVG','PNG');
% plt.SaveFigure(saveFig,FigsDir,filename,'SVG');

end

