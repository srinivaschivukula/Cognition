function CrossMahalanobisDistance(AllTrialData,Labels,TimeWindow,OutDir,Tag,varargin)
% Train on dim 1 to classify dim 2, test on how well it generalizes.
% e.g. train on nancy to diff bp, test how well bp diff'd in Tyson.

[~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh','fdr');
% [~,NIter] = Utilities.ProcVarargin(varargin,'NReps',100);
[~,Phase] = Utilities.ProcVarargin(varargin,'Phase','Go');
% [~,BasePhase] = Utilities.ProcVarargin(varargin,'BasePhase','Delay');
% [~,BaseWindow] = Utilities.ProcVarargin(varargin,'BaseWindow',[1.5 2.5]);
% [~,Dims] = Utilities.ProcVarargin(varargin,'SplitDims',{});
[~,DimNames] = Utilities.ProcVarargin(varargin,'DimNames',{});
% [~,SplitValNames] = Utilities.ProcVarargin(varargin,'SplitValNames',{});
[~,PBPSpec] = Utilities.ProcVarargin(varargin,'PBPSpec');
[~,sigUnits] = Utilities.ProcVarargin(varargin,'sigUnits');
[~,ActSpec] = Utilities.ProcVarargin(varargin,'ActSpec');
% [~,ConditionFieldLabel] = Utilities.ProcVarargin(varargin,'ConditionFieldLabel','Condition');
% [~,ObsPerCond] = Utilities.ProcVarargin(varargin,'ObsPerCond',10);
[~,FileTag] = Utilities.ProcVarargin(varargin,'FileTag','');
[~,saveFig] = Utilities.ProcVarargin(varargin,'saveFig');
[~,withShuffle] = Utilities.ProcVarargin(varargin,'withShuffle');
[~,Conditions] = Utilities.ProcVarargin(varargin,'Conditions',{});
[~,NumRandIterations] = Utilities.ProcVarargin(varargin,'NumRandIterations',10);
[~,NumRandUnits] = Utilities.ProcVarargin(varargin,'NumRandUnits',50);


StyleArgs = {'FontSize',9,'FontName','Arial','FontWeight','bold'};

DecodeLabels = [Conditions{:}];

%% Only include units with tunign to at least one of the relevant variables
Dates = Analyze.returnUniqueFieldValues(AllTrialData,'Date');
if PBPSpec
    ASF = Analyze.FaceScratch3.FindSignificantUnitsAction(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,'Dates',Dates,'UnitGroup','PBPSpec');
elseif ActSpec
    ASF = Analyze.FaceScratch3.FindSignificantUnitsAction(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,'Dates',Dates,'UnitGroup','ActSpec');
elseif sigUnits
    ASF = Analyze.FaceScratch3.FindSignificantUnitsAction(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,'Dates',Dates,'UnitGroup','sigUnits');
else
    ASF = Analyze.FaceScratch3.FindSignificantUnitsAction(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,'Dates',Dates,'UnitGroup','All');
end

%%

opts.NumShuffles=200;
opts.cvOptions.ValidationType='ClassicCrossValidation';
% opts.cvOptions.NReps=30;
opts.cvOptions.NReps=10;
opts.cvOptions.NFolds=10;
opts.dec=Predictor.FWClassifier(@Analyze.FaceScratch3.BasicClassifier);

%% first do within and across format classification of Action Selectivity

shuffles = setdiff(perms([1 2 3 4]),[1 2 3 4],'rows');

shufcvAccuracy=[];cvAccuracy=[];numUnits=[];
for dateIdx=1:length(Dates)
    
    cDate = Dates{dateIdx};
    dateCode = Blackrock.Helper.date2unitId(cDate);
    ccASF = ASF(ASF(:,4)==dateCode,:);
    numUnits(dateIdx,:) = size(ccASF,1);
    
    randAcc=[];
    for randIdx=1:1
        
        randNum = randperm(size(ccASF,1),ceil((1)*size(ccASF,1)))';
        cASF = ASF(randNum,:);
        
        for i=1:length(Conditions)
            cCond=Conditions{i};
            cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
            L{i}=Analyze.returnFieldValues(cTrialData,'Action');
            FR{i}=squeeze(Analyze.getNeuralData(cTrialData,cASF,TimeWindow));
            [FR{i},L{i}]=Analyze.SampleFeaturePopulation(FR{i},L{i},'Type','Basic');
            [FRCell{i},LabelCell{i}]=DataConvert.list2cell(FR{i},L{i});
            [FRTrue{i},LabelTrue{i}]=DataConvert.cell2list(FRCell{i});
        end
        tmpAcc = DPA.CrossValidatedDistance(FRTrue);
        randAcc(:,:,randIdx) = cell2mat(tmpAcc);
        
        % shuffling stuff
        shufVals=[];
        for jj=1:length(FRTrue)
            gr1FR = FRTrue(1:jj);
            gr2FR = FRTrue(jj+1:length(FRTrue));
            
            abc = []; tmpShufMat=[];
            for shufIdx=1:length(shuffles)
                cShuf = shuffles(shufIdx,:);
                for kk=1:length(cShuf)
                    abc = cat(2,abc,((cShuf(kk)*10)-9):(cShuf(kk)*10));
                end
                cshufFR = [gr1FR cellfun(@(x)x(abc',:),gr2FR,'UniformOutput',0)];
                tmpShufMat(:,:,shufIdx) = cell2mat(DPA.CrossValidatedDistance(cshufFR));
            end
            meanShufMat = mean(tmpShufMat,3);
            shufVals(jj,:) = meanShufMat(jj,:);
        end
        
        shufRandAccuracy(:,:,randIdx) = shufVals;
    end
    shufcvAccuracy(:,:,dateIdx) = mean(shufRandAccuracy,3);
    cvAccuracy(:,:,dateIdx)= mean(randAcc,3);
end

%% Plotting

T={'Fc-Fs','Fc-Oc', 'Fc-Os','Fs-Oc', 'Fs-Os','Oc-Os'};
IDX=[2 3 4 7 8 12];
% IDX=[1:16];

clear Vals
for i=1:length(IDX)
    [I,J] = ind2sub([4 4],IDX(i));
    Vals(:,i)= squeeze(cvAccuracy(J,I,:));
end

Vals = abs(Vals);
% Vals = Vals./numUnits;

clr=[];
clr(1:6,:) = repmat([185 185 185]./255,6,1);
StyleArgs = {'fontweight','bold','fontname','helvetica','fontsize',9};

plt.fig('units','inches','width',5,'height',6,'font','Arial','fontsize',12);
pnl = panel();  pnl.margin=15; pnl.pack('v',2);
pnl(1).select()

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

clear Vals2
for i=1:length(IDX)
    [I,J] = ind2sub([4 4],IDX(i));
    Vals2(:,i)= squeeze(shufcvAccuracy(J,I,:));
end

Vals2 = abs(Vals2);
% Vals = Vals./numUnits;

H2 = notBoxPlot(Vals2,[1 2 3 5 6 8],'jitter',0.3);
set([H2.data],...
    'MarkerFaceColor','r',...
    'markerEdgeColor','r',...
    'MarkerSize',2)
set([H2.mu],'color','r')
J=clr;
for ii=1:length(H2)
    set(H2(ii).sdPtch,'FaceColor',[254 224 210]./255,...
        'EdgeColor','none')
    set(H2(ii).semPtch,'FaceColor',J(ii,:)*0.3,'Visible','off',...
        'EdgeColor','none')
end

% plot([1 2 3 5 6 8],muOld,'rd','MarkerFaceColor',[1 0 0],'MarkerSize',2.5)
axis tight
ax1 = gca;
set(ax1,'box','off',StyleArgs{:});
ax1.XRuler.Axle.LineStyle = 'none';
xticklabels(T); xtickangle(45)
set(gca,'TickDir','out');
set(ax1,'TickLength',[0.005 0.025]);
ylabel('cv Mahalanobis Distance')
set(ax1,'box','off',StyleArgs{:});
NC = 14; xlim([.25 NC+.25]);
ylim([0 10])

%% To plot even the cross-validated distances use this:

T={'Fc-Fc','Fc-Fs','Fc-Oc','Fc-Os','Fs-Fs','Fs-Oc','Fs-Os','Oc-Oc','Oc-Os','Os-Os'};
IDX=[1 2 3 4 6 7 8 11 12 16];

clear Vals
for i=1:length(IDX)
    [I,J] = ind2sub([4 4],IDX(i));
    Vals(:,i)= squeeze(cvAccuracy(J,I,:));
end

Vals = abs(Vals);
% Vals = Vals./numUnits;

clr=[];
clr(1:10,:) = repmat([185 185 185]./255,10,1);
StyleArgs = {'fontweight','bold','fontname','helvetica','fontsize',9};

plt.fig('units','inches','width',5,'height',6,'font','Arial','fontsize',12);
pnl = panel();  pnl.margin=15; pnl.pack('v',2);
pnl(1).select()

H = notBoxPlot(Vals,[1 2 3 4 6 7 8 10 11 13],'jitter',0.3);
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

clear Vals2
for i=1:length(IDX)
    [I,J] = ind2sub([4 4],IDX(i));
    Vals2(:,i)= squeeze(shufcvAccuracy(J,I,:));
end

Vals2 = abs(Vals2);
% Vals = Vals./numUnits;

H2 = notBoxPlot(Vals2,[1 2 3 4 6 7 8 10 11 13],'jitter',0.3);
set([H2.data],...
    'MarkerFaceColor','r',...
    'markerEdgeColor','r',...
    'MarkerSize',2)
set([H2.mu],'color','r')
J=clr;
for ii=1:length(H2)
    set(H2(ii).sdPtch,'FaceColor',[254 224 210]./255,...
        'EdgeColor','none')
    set(H2(ii).semPtch,'FaceColor',J(ii,:)*0.3,'Visible','off',...
        'EdgeColor','none')
end

% plot([1 2 3 5 6 8],muOld,'rd','MarkerFaceColor',[1 0 0],'MarkerSize',2.5)
axis tight
ax1 = gca;
set(ax1,'box','off',StyleArgs{:});
ax1.XRuler.Axle.LineStyle = 'none';
xticklabels(T); xtickangle(45)
set(gca,'TickDir','out');
set(ax1,'TickLength',[0.005 0.025]);
ylabel('Distance')
set(ax1,'box','off',StyleArgs{:});
NC = 14; xlim([.25 NC+.25]);
ylim([0 10])


%%
filename = sprintf('CrossDecode-%s-vs-%s',DimNames{1},DimNames{2});
if PBPSpec
    var = 'PBPSpec';
    filename = [filename '-PBPSpec'];
    suptitle(sprintf('PBPSpec units only (%d)',size(ASF,1)));
elseif ActSpec
    var = 'ActSpec';
    filename = [filename '-ActSpec'];
    suptitle(sprintf('ActSpec units only (%d)',size(ASF,1)));
elseif sigUnits
    var = 'sigUnits';
    filename = [filename '-sigUnits'];
    suptitle(sprintf('All significant units (%d)',size(ASF,1)));
else
    var = 'allUnits';
    filename = [filename '-allUnits'];
    suptitle(sprintf('All units (%d)',size(ASF,1)));
    
end

if ~isempty(FileTag)
    filename = [filename '-' FileTag];
end

FigsDir = fullfile(env.get('results'),'FaceScratch3','SUAnal',[Tag '-Go'],'PopData','Mine');
filename = [filename];
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');
% plt.SaveFigure(saveFig,FigsDir,[filename],'SVGI');

end

