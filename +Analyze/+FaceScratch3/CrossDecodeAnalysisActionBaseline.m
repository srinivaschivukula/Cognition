function CrossDecodeAnalysisActionBaseline(AllTrialData,Labels,TimeWindow,OutDir,Tag,varargin)
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

shufcvAccuracy=[];cvAccuracy=[];
for dateIdx=1:length(Dates)
    
    cDate = Dates{dateIdx};
    dateCode = Blackrock.Helper.date2unitId(cDate);
    cASF = ASF(ASF(:,4)==dateCode,:);
    
    %     tmpCVaccuracy=[]; tmpShufCVaccuracy=[];
    %     for randIDX = 1:NumRandIterations
    %
    %         cASF = ccASF(randperm(size(ccASF,1),NumRandUnits)',:);
    
    for i=1:length(Conditions)
        cCond=Conditions{i};
        cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
        
        L{i}=Analyze.returnFieldValues(cTrialData,'Action');
        
        if iscell(TimeWindow)
            FR{i}=squeeze(Analyze.getNeuralData(cTrialData,cASF,TimeWindow{i}));
        else
            FR{i}=squeeze(Analyze.getNeuralData(cTrialData,cASF,TimeWindow));
        end
        
        [FR{i},L{i}]=Analyze.SampleFeaturePopulation(FR{i},L{i},'Type','Basic');
        
        if 1==1
            [FRCell{i},LabelCell{i}]=DataConvert.list2cell(FR{i},L{i});
            [FRTrue{i},LabelTrue{i}]=DataConvert.cell2list(FRCell{i});
        else
            [FRCell{i},LabelCell{i}]=DataConvert.list2cell(FR{i},L{i});
            
            mu=cellfun(@mean,FRCell{i},'UniformOutput',false);
            mu=cat(1,mu{:});
            mn=min(mu,[],1);
            mx=max(mu,[],1);
            FR{i}=(FR{i}-mn)./(mx-mn);
            [FRCell{i},LabelCell{i}]=DataConvert.list2cell(FR{i},L{i});
            [FRTrue{i,1},LabelTrue{i}]=DataConvert.cell2list(FRCell{i});
        end
    end
    %%
    %         Shuffles = perms(1:4);
    
    A=[2 3 4];B=[1 3 4]; C=[1 2 3]; D=[1 2 4];
    [A,B,C,D]=ndgrid(A,B,C,D);
    tmpGrid=[A(:),B(:),C(:),D(:)];
    for i=1:size(tmpGrid,1)
        N(i)=length(unique(tmpGrid(i,:)));
    end
    
    Shuffles=tmpGrid(find(N==4),:);
    %%
    for i=1:length(Conditions)
        
        ShuffleRank=1;
        opts.dec.Train('TrainingData',...
            {FRCell{i},[1:4],LabelCell{i}},'CrossValidate','cvOptions',opts.cvOptions);
        
        for j=1:4
            disp([i j])
            if i==j
                tmp.BootCI{i,j}=opts.dec.getBootCI;
                tmp.ShuffleRank(i,i)=ShuffleRank;
                
                tmp.cvAccuracy(i,j)=opts.dec.Results.cvAccuracy;
                tmp.confusionMat{i,j}=opts.dec.Results.ConfusionMat;
                
            else
                LabelEst=opts.dec.PredictBatch(FRTrue{j});
                foo=Analyze.Semantic.CrossDecodeHelper.ComputeCrossStats(LabelEst,LabelTrue{j},LabelCell{j});
                tmp.ShuffleRank(i,j)= foo.ShuffleRank; tmp.cvAccuracy(i,j)=foo.cvAccuracy;
                tmp.confusionMat{i,j}=foo.confusionMat; tmp.BootCI{i,j}=foo.BootCI;
                
                %%
                for kk=1:size(Shuffles,1)
                    LabelShuff=Analyze.Semantic.CrossDecodeHelper.ShuffleConditionLabels(LabelTrue{j},Shuffles(kk,:));
                    foo=Analyze.Semantic.CrossDecodeHelper.ComputeCrossStats(LabelEst,LabelShuff',LabelCell{j});
                    tmp.shufShuffleRank(i,j,kk)= foo.ShuffleRank; tmp.shufcvAccuracy(i,j,kk)=foo.cvAccuracy;
                    tmp.shufconfusionMat{i,j,kk}=foo.confusionMat; tmp.shufBootCI{i,j,kk}=foo.BootCI;
                end
            end
        end
    end
    
    %         tmpCVaccuracy= cat(3,tmpCVaccuracy,tmp.cvAccuracy);
    %         tmpShufCVaccuracy = cat(3,tmpShufCVaccuracy,tmp.shufcvAccuracy);
    %     end
    
    shufcvAccuracy = cat(3,shufcvAccuracy,tmp.shufcvAccuracy);
    cvAccuracy = cat(3,cvAccuracy,tmp.cvAccuracy);
    
end

%% now do within and across format classification of condition (person-body-part)

Conditions = {{'NancyPinchCheek','NancyPinchShoulder','TysonPinchCheek','TysonPinchShoulder'},...
    {'NancyPressCheek','NancyPressShoulder','TysonPressCheek','TysonPressShoulder',},...
    {'NancyRubCheek','NancyRubShoulder','TysonRubCheek','TysonRubShoulder'},...
    {'NancyTapCheek','NancyTapShoulder','TysonTapCheek','TysonTapShoulder'}};

shufcvAccuracy2=[];cvAccuracy2=[];
for dateIdx=1:length(Dates)
    
    cDate = Dates{dateIdx};
    dateCode = Blackrock.Helper.date2unitId(cDate);
    ccASF = ASF(ASF(:,4)==dateCode,:);
    
    tmpCVaccuracy=[]; tmpShufCVaccuracy=[];
    for randIDX = 1:NumRandIterations
        
        cASF = ccASF(randperm(size(ccASF,1),NumRandUnits)',:);
        
        for i=1:length(Conditions)
            cCond=Conditions{i};
            cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
            
            L{i}=Analyze.returnFieldValues(cTrialData,'PersonBodyPart');
            if iscell(TimeWindow)
                FR{i}=squeeze(Analyze.getNeuralData(cTrialData,cASF,TimeWindow{i}));
            else
                FR{i}=squeeze(Analyze.getNeuralData(cTrialData,cASF,TimeWindow));
            end
            
            [FR{i},L{i}]=Analyze.SampleFeaturePopulation(FR{i},L{i},'Type','Basic');
            
            if 1==0
                [FRCell{i},LabelCell{i}]=DataConvert.list2cell(FR{i},L{i});
                [FRTrue{i},LabelTrue{i}]=DataConvert.cell2list(FRCell{i});
            else
                [FRCell{i},LabelCell{i}]=DataConvert.list2cell(FR{i},L{i});
                
                mu=cellfun(@mean,FRCell{i},'UniformOutput',false);
                mu=cat(1,mu{:});
                mn=min(mu,[],1);
                mx=max(mu,[],1);
                FR{i}=(FR{i}-mn)./(mx-mn);
                [FRCell{i},LabelCell{i}]=DataConvert.list2cell(FR{i},L{i});
                [FRTrue{i},LabelTrue{i}]=DataConvert.cell2list(FRCell{i});
            end
        end
        %%
        Shuffles = perms(1:4);
        %%
        % rng(123);
        for i=1:length(Conditions)
            %     opts.cvOptions.ValidationType='ResamplingCrossValidation';
            %     [ShuffleRank,AllInfo]=opts.dec.ShuffleTest('TrainingData',...
            %         {FRCell{i},[1:5],LabelCell{i}},'CrossValidate','cvOptions',opts.cvOptions);
            ShuffleRank=1;
            opts.dec.Train('TrainingData',...
                {FRCell{i},[1:4],LabelCell{i}},'CrossValidate','cvOptions',opts.cvOptions);
            
            for j=1:4
                disp([i j])
                if i==j
                    tmp.BootCI{i,j}=opts.dec.getBootCI;
                    tmp.ShuffleRank(i,i)=ShuffleRank;
                    
                    tmp.cvAccuracy(i,j)=opts.dec.Results.cvAccuracy;
                    tmp.confusionMat{i,j}=opts.dec.Results.ConfusionMat;
                    
                else
                    
                    LabelEst=opts.dec.PredictBatch(FRTrue{j});
                    foo=Analyze.Semantic.CrossDecodeHelper.ComputeCrossStats(LabelEst,LabelTrue{j},LabelCell{j});
                    tmp.ShuffleRank(i,j)= foo.ShuffleRank; tmp.cvAccuracy(i,j)=foo.cvAccuracy;
                    tmp.confusionMat{i,j}=foo.confusionMat; tmp.BootCI{i,j}=foo.BootCI;
                    
                    %%
                    for kk=1:size(Shuffles,1)
                        LabelShuff=Analyze.Semantic.CrossDecodeHelper.ShuffleConditionLabels(LabelTrue{j},Shuffles(kk,:));
                        foo=Analyze.Semantic.CrossDecodeHelper.ComputeCrossStats(LabelEst,LabelShuff',LabelCell{j});
                        tmp.shufShuffleRank(i,j,kk)= foo.ShuffleRank; tmp.shufcvAccuracy(i,j,kk)=foo.cvAccuracy;
                        tmp.shufconfusionMat{i,j,kk}=foo.confusionMat; tmp.shufBootCI{i,j,kk}=foo.BootCI;
                    end
                end
            end
        end
        tmpCVaccuracy= cat(3,tmpCVaccuracy,tmp.cvAccuracy);
        tmpShufCVaccuracy = cat(3,tmpShufCVaccuracy,tmp.shufcvAccuracy);
    end
    shufcvAccuracy2 = cat(3,shufcvAccuracy2,tmpShufCVaccuracy);
    cvAccuracy2 = cat(3,cvAccuracy2,mean(tmpCVaccuracy,3));
    
end
%% plots
% first plot cross-classification of Action selectivity

IDX=[1:4:13 2:4:14 3:4:15 4:4:16];
clear plotData
for i=1:length(IDX)
    [I,J] = ind2sub([4 4],IDX(i));
    plotData(:,i)= squeeze(cvAccuracy(I,J,:));
end

StyleArgs={'FontWeight','bold','FontSize',9,'FontName','Helvetica'};

clr=[];
clr([1:4],:) = repmat([237 0 38]./255,4,1);
clr([5:8],:) = repmat([31 120 180]./255,4,1);
clr([9:12],:) = repmat([240 228 66]./255,4,1);
clr([13:16],:) = repmat([189 189 189]./255,4,1);

plt.fig('units','inches','width',17,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();
xtickangle(45); set(gca, 'fontweight','bold');

% H = notBoxPlotMod(plotData,[1:4 7:10 13:16 19:22],'jitter',0.3,'useCI');
H = notBoxPlot(plotData,[1:4 6:9 11:14 16:19],'jitter',0.3);

% if withShuffle
hold on;
IDX=[1:4:13 2:4:14 3:4:15 4:4:16];
clear plotData
for i=1:length(IDX)
    [I,J] = ind2sub([4 4],IDX(i));
    plotData(:,i)= squeeze(shufcvAccuracy(I,J,:));
end

clr1=[];
clr1([1:4 9:12],:) = repmat([25 25 25]./255,8,1);
clr1([5:8 13:16],:) = repmat([25 25 25]./255,8,1);

H1 = notBoxPlotMod(plotData,[1:4 6:9 11:14 16:19],'jitter',0.3,'useCI');
%     H1 = notBoxPlot(plotData,[1:4 7:10 13:16 19:22],'jitter',0.3);

%     set([H1.data],...
%         'MarkerFaceColor',[1,1,1]*0.35,...
%         'markerEdgeColor',[1,1,1]*0.35,...
%         'MarkerSize',1.5)
set([H1.mu],'color',[227 74 51]./255)
J=clr1;
for ii=1:length(H)
    set(H1(ii).sdPtch,'FaceColor',J(ii,:)*0.3,...
        'EdgeColor','none')
    set(H1(ii).semPtch,'FaceColor',J(ii,:)*0.3,'Visible','off',...
        'EdgeColor','none');
end
% end

set([H.data],...
    'MarkerFaceColor',[1,1,1]*0.35,...
    'markerEdgeColor',[1,1,1]*0.35,...
    'MarkerSize',2)
set([H.mu],'color','k')
J=clr;
for ii=1:length(H)
    set(H(ii).sdPtch,'FaceColor',J(ii,:),...
        'EdgeColor','none')
    set(H(ii).semPtch,'FaceColor',J(ii,:),...
        'EdgeColor','none');
end
axis tight
ax1 = gca;
ax1.XRuler.Axle.LineStyle = 'none';
cLabels = {'Fc->Fc','Fc->Fs','Fc->Oc','Fc->Os',...
    'Fs->Fc','Fs->Fs','Fs->Oc','Fs->Os',...
    'Oc->Fc','Oc->Fs','Oc->Oc','Oc->Os',...
    'Os->Fc','Os->Fs','Os->Oc','Os->Os'};
% cLabels = {'Fc->Fc','Fc->Oc','Fc->Fs','Fc->Os',...
%     'Oc->Fc','Oc->Oc','Oc->Fs','Oc->Os',...
%     'Fs->Fc','Fs->Oc','Fs->Fs','Fs->Os',...
%     'Os->Fc','Os->Oc','Os->Fs','Os->Os'};
xticklabels(cLabels);
set(ax1,'box','off',StyleArgs{:});

set(gca,'TickDir','out');
set(ax1,'TickLength',[0.00025 0.0]);

line([ax1.XLim(1) max(H(end).sdPtch.Vertices(:,1))],[25 25],'Color',[115 115 115]./255,'LineStyle','--');
text(max(H(end).sdPtch.Vertices(:,1)),11,'Chance','HorizontalAlignment','right','VerticalAlignment','bottom',...
    StyleArgs{:},'Color',[115 115 115]./255);

ylabel('Decode Accuracy (%)')
% xlabel('Digit')
title('Classifying Action',StyleArgs{:})
NC = 22; xlim([.25 NC+.25]);
ylim([0 105])


%% next plot cross-classification of Person-BodyPart

IDX=[1:4:13 2:4:14 3:4:15 4:4:16];
clear plotDataBP
for i=1:length(IDX)
    [I,J] = ind2sub([4 4],IDX(i));
    plotDataBP(:,i)= squeeze(cvAccuracy2(I,J,:));
end


clr=[];
clr([1:4],:) = repmat([44 162 95]./255,4,1);
clr([5:8],:) = repmat([117 107 177]./255,4,1);
clr([9:12],:) = repmat([255 127 66]./255,4,1);
clr([13:16],:) = repmat([255 255 179]./255,4,1);

% clr=[];
% clr([1:4 9:12],:) = repmat([44 162 95]./255,8,1);
% clr([5:8 13:16],:) = repmat([117 107 177]./255,8,1);

pnl(1,2).select();
xtickangle(45); set(gca, 'fontweight','bold');

% H = notBoxPlotMod(plotDataBP,[1:4 7:10 13:16 19:22],'jitter',0.3,'useCI');
H = notBoxPlot(plotDataBP,[1:4 6:9 11:14 16:19],'jitter',0.3);

% if withShuffle
hold on;
IDX=[1:4:13 2:4:14 3:4:15 4:4:16];
clear plotDataBP
for i=1:length(IDX)
    [I,J] = ind2sub([4 4],IDX(i));
    plotDataBP(:,i)= squeeze(shufcvAccuracy2(I,J,:));
end

clr1=[];
clr1([1:4 9:12],:) = repmat([25 25 25]./255,8,1);
clr1([5:8 13:16],:) = repmat([25 25 25]./255,8,1);

H1 = notBoxPlotMod(plotDataBP,[1:4 6:9 11:14 16:19],'jitter',0.3,'useCI');

set([H1.mu],'color',[227 74 51]./255)
J=clr1;
for ii=1:length(H)
    set(H1(ii).sdPtch,'FaceColor',J(ii,:)*0.3,...
        'EdgeColor','none')
    set(H1(ii).semPtch,'FaceColor',J(ii,:)*0.3,...
        'EdgeColor','none');
end
% end

set([H.data],...
    'MarkerFaceColor',[1,1,1]*0.35,...
    'markerEdgeColor',[1,1,1]*0.35,...
    'MarkerSize',2)
set([H.mu],'color','k')
J=clr;
for ii=1:length(H)
    set(H(ii).sdPtch,'FaceColor',J(ii,:),...
        'EdgeColor','none')
    set(H(ii).semPtch,'FaceColor',J(ii,:)*0.3,'FaceAlpha',0,...
        'EdgeColor','none');
end
axis tight
ax1 = gca;
ax1.XRuler.Axle.LineStyle = 'none';
cLabels = {'Pi->Pi','Pi->Pr','Pi->Ru','Pi->Tap',...
    'Pr->Pi','Pr->Pr','Pr->Ru','Pr->Tap',...
    'Ru->Pi','Ru->Pr','Ru->Ru','Ru->Tap',...
    'Ta->Pi','Ta->Pr','Ta->Ru','Ta->Tap'};
xticklabels(cLabels);
set(ax1,'box','off',StyleArgs{:});

set(gca,'TickDir','out');
set(ax1,'TickLength',[0.00025 0.0]);

line([ax1.XLim(1) max(H(end).sdPtch.Vertices(:,1))],[25 25],'Color',[115 115 115]./255,'LineStyle','--');
text(max(H(end).sdPtch.Vertices(:,1)),11,'Chance','HorizontalAlignment','right','VerticalAlignment','bottom',...
    StyleArgs{:},'Color',[115 115 115]./255);

ylabel('Decode Accuracy (%)')
% xlabel('Digit')
title('Classifying Condition',StyleArgs{:})
NC = 22; xlim([.25 NC+.25]);
ylim([0 110])

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

