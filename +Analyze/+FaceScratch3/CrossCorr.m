function CrossCorr(OutDir,Tag,varargin)
% Plots confusion matrix (decode analysis) for decoding each (relevant)
% condition.

% [~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh',0.001);
[~,Pmethod]  = Utilities.ProcVarargin(varargin,'Pmethod','fdr');
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');
[varargin,AllTrialData]   = Utilities.ProcVarargin(varargin,'AllTrialData',{});
[varargin,PhaseStart]   = Utilities.ProcVarargin(varargin,'PhaseStart',[]);
[varargin,PhaseDur]   = Utilities.ProcVarargin(varargin,'PhaseDur',[]);
[varargin,BaselineWindow]   = Utilities.ProcVarargin(varargin,'BaselineWindow',[]);
[varargin,Phase]   = Utilities.ProcVarargin(varargin,'Phase','Go');


StyleArgs = {'FontSize',9,'FontName','Arial','FontWeight','bold'};

%% Only include units with tunign to at least one of the relevant variables
datafile = fullfile(env.get('results'),'FaceScratch3','SUAnal',[Tag '-Go'],'PopData','p1-AllDays.mat');
load(datafile); % Loads PlotData

ASF = vertcat(PlotData{1}.Unit{:});
% timeWindow = [2.5 6];
timeWindow = [0.5 3];
% [1 4]

clear FR
cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase);
Condition={ {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek'},...
    {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'},...
    {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'},...
    {'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'}};

Action=Analyze.returnUniqueFieldValues(cTrialData,'Action');
Action=Action(cellfun(@isempty,strfind(Action,'Null')));

Dates = Analyze.returnUniqueFieldValues(PlotData,'TaskDate');
dateCodes = Blackrock.Helper.date2unitId(Dates);

% shuffles = perms(1:4);
% shuffles(ismember(shuffles,[1 2 3 4],'rows'),:)=[];

A=[2 3 4];B=[1 3 4]; C=[1 2 3]; D=[1 2 4];
[A,B,C,D]=ndgrid(A,B,C,D);
tmp=[A(:),B(:),C(:),D(:)];
for i=1:size(tmp,1)
    N(i)=length(unique(tmp(i,:)));
end

shuffles=tmp(find(N==4),:);

FR=cell(1,length(Dates)); FRdm=cell(1,length(Dates));
FRshuff=cell(1,length(Dates));

for sessionIdx = 1:length(Dates)
    cASF = ASF(ASF(:,4)==dateCodes(sessionIdx),:);
    
    clear sessionFR;
    sessionFR = cell(1,size(cASF,1));
    for unitIdx = 1:size(cASF,1)
        cUnit = cASF(unitIdx,:);
        
        unitFR=[];unitFRdm = [];
        for cIDX=1:length(Condition)
            cCond=Condition{cIDX};
            
            tmpFR=[];
            for actionIdx = 1:length(Action)
                cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond,'Action',Action{actionIdx});
                tmpVal = squeeze(Analyze.getNeuralData(cTrialData,cUnit,timeWindow));
                tmpVal = tmpVal(~isnan(tmpVal));
                tmpFR = cat(1,tmpFR,mean(tmpVal));
            end
            unitFR = cat(2,unitFR,tmpFR);
            unitFRdm = cat(2,unitFRdm,tmpFR-mean(tmpFR));
        end
        sessionFR{unitIdx} = unitFR;
        
        FR{sessionIdx}=cat(1,FR{sessionIdx},unitFR);
        FRdm{sessionIdx}=cat(1,FRdm{sessionIdx},unitFRdm);
    end
    
    %% do the shuffling stuff
    shuffCorr=[];
    for i=1:4
        for j=1:4
            
            tmp1 = cellfun(@(x)x(:,i),sessionFR,'UniformOutput',false);
            %             tmp1 = cellfun(@(x)x-mean(x),tmp1,'UniformOutput',false);
            tmp1 = cat(1,tmp1{:});
            
            for kk = 1:size(shuffles,1)
                cShuff = shuffles(kk,:);
                tmp2 = cellfun(@(x)x(cShuff',j),sessionFR,'UniformOutput',false);
                %                 tmp2 = cellfun(@(x)x-mean(x),tmp2,'UniformOutput',false);
                tmp2 = cat(1,tmp2{:});
                
                tmpVals=[];
                for tmpIdx = 1:4:length(tmp2)
                    tmpVals = cat(1,tmpVals,corr(tmp1([tmpIdx:(tmpIdx+3)]),tmp2([tmpIdx:(tmpIdx+3)])));
                end
                %                 shuffCorr(i,j,kk) = corr(tmp1,tmp2);
                shuffCorr(i,j,kk) = mean(tmpVals);
            end
        end
    end
    FRshuff{sessionIdx} = mean(shuffCorr,3);
end

%% for the shuffles

%% plotting the data

clr=[];
clr(1:6,:) = repmat([185 185 185]./255,6,1);
StyleArgs = {'fontweight','bold','fontname','helvetica','fontsize',9};

plt.fig('units','inches','width',5,'height',6,'font','Arial','fontsize',12);
pnl = panel();  pnl.margin=15; pnl.pack('v',2);

MEAS = FR;

pnl(1).select()
cvAccuracy =[];
cvAccuracyShuff=[];
for i=1:length(MEAS)
    cvAccuracy(:,:,i) = corr(MEAS{i});
    cvAccuracyShuff(:,:,i) = FRshuff{i};
end

T={'Fc-Fs','Fc-Oc', 'Fc-Os','Fs-Oc', 'Fs-Os','Oc-Os'};
IDX=[2 3 4 7 8 12];


clear Vals
for i=1:length(IDX)
    [I,J] = ind2sub([4 4],IDX(i));
    Vals(:,i)= squeeze(cvAccuracy(I,J,:));
    hold on;
        ValsShuff(:,i)= squeeze(cvAccuracyShuff(I,J,:));
end
muOld = mean(Vals,1);
H = notBoxPlot(Vals,[1 2 3 5 6 8],'jitter',0.3);
set([H.data],...
    'MarkerFaceColor',[1,1,1]*0.35,...
    'markerEdgeColor',[1,1,1]*0.35,...
    'MarkerSize',2)
set([H.mu],'color','k')
J=clr;
for ii=1:length(H)
    set(H(ii).sdPtch,'FaceColor',J(ii,:),...
        'EdgeColor','none')
    set(H(ii).semPtch,'FaceColor',J(ii,:)*0.3,'Visible','off',...
        'EdgeColor','none')
end

hold on 
H1 = notBoxPlot(ValsShuff,[1 2 3 5 6 8],'jitter',0.3);
set([H1.data],...
    'MarkerFaceColor','r',...
    'markerEdgeColor','r',...
    'MarkerSize',2)
set([H1.mu],'color','r')
% J=clr;
for ii=1:length(H1)
    set(H1(ii).sdPtch,'FaceColor',[254 224 210]./255,...
        'EdgeColor','none')
    set(H1(ii).semPtch,'FaceColor',J(ii,:)*0.3,'Visible','off',...
        'EdgeColor','none')
end


axis tight
ax1 = gca;
set(ax1,'box','off',StyleArgs{:});
ax1.XRuler.Axle.LineStyle = 'none';
xticklabels(T); xtickangle(45)
set(gca,'TickDir','out');
set(ax1,'TickLength',[0.005 0.025]);
ylabel('Correlation')
set(ax1,'box','off',StyleArgs{:});
NC = 14; xlim([.25 NC+.25]);
ylim([-0.2 1])

line([0 max(H(end).sdPtch.Vertices(:,1))],[0 0],'Color',[115 115 115]./255,'LineStyle','--');
text(max(H(end).sdPtch.Vertices(:,1)),0,'Chance','HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontWeight','bold','FontSize',9,'FontName','Arial','Color','k');


pnl(2).select()
cvAccuracy =[];
for i=1:length(MEAS)
    cvAccuracy(:,:,i) = partialcorr(MEAS{i});
end

T={'Fc-Fs','Fc-Oc', 'Fc-Os','Fs-Oc', 'Fs-Os','Oc-Os'};
IDX=[2 3 4 7 8 12];


clear Vals
for i=1:length(IDX)
    [I,J] = ind2sub([4 4],IDX(i));
    Vals(:,i)= squeeze(cvAccuracy(I,J,:));
        ValsShuff(:,i)= squeeze(cvAccuracyShuff(I,J,:));
end

H = notBoxPlot(Vals,[1 2 3 5 6 8],'jitter',0.3);
set([H.data],...
    'MarkerFaceColor',[1,1,1]*0.35,...
    'markerEdgeColor',[1,1,1]*0.35,...
    'MarkerSize',2)
set([H.mu],'color','k')
J=clr;
for ii=1:length(H)
    set(H(ii).sdPtch,'FaceColor',J(ii,:),...
        'EdgeColor','none')
    set(H(ii).semPtch,'FaceColor',J(ii,:)*0.3,'Visible','off',...
        'EdgeColor','none')
end
hold on;
% plot([1 2 3 5 6 8],muOld,'rd','MarkerFaceColor',[1 0 0],'MarkerSize',2.5)
axis tight
ax1 = gca;
set(ax1,'box','off',StyleArgs{:});
ax1.XRuler.Axle.LineStyle = 'none';
xticklabels(T); xtickangle(45)
set(gca,'TickDir','out');
set(ax1,'TickLength',[0.005 0.025]);
ylabel('Partial Correlation')
set(ax1,'box','off',StyleArgs{:});
NC = 14; xlim([.25 NC+.25]);
ylim([-0.2 1])
line([0 max(H(end).sdPtch.Vertices(:,1))],[0 0],'Color',[115 115 115]./255,'LineStyle','--');
text(max(H(end).sdPtch.Vertices(:,1)),0,'Chance','HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontWeight','bold','FontSize',9,'FontName','Arial','Color','k');
hold on; 



FigsDir = fullfile(env.get('results'),'FaceScratch3','SUAnal',[Tag '-Go'],'PopData','Mine');
filename = 'CrossCorrelation';
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');



end