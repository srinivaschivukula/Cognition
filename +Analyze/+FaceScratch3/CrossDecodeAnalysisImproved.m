function CrossDecodeAnalysisImproved(AllTrialData,Labels,DecodeLabels,TimeWindow,OutDir,Tag,varargin)
% Train on dim 1 to classify dim 2, test on how well it generalizes.
% e.g. train on nancy to diff bp, test how well bp diff'd in Tyson.

[~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh','fdr');
% [~,NIter] = Utilities.ProcVarargin(varargin,'NReps',100);
[~,Phase] = Utilities.ProcVarargin(varargin,'Phase','Go');
[~,BasePhase] = Utilities.ProcVarargin(varargin,'BasePhase','Delay');
[~,BaseWindow] = Utilities.ProcVarargin(varargin,'BaseWindow',[1.5 2.5]);
[~,Dims] = Utilities.ProcVarargin(varargin,'SplitDims',{});
[~,DimNames] = Utilities.ProcVarargin(varargin,'DimNames',{});
[~,SplitValNames] = Utilities.ProcVarargin(varargin,'SplitValNames',{});
[~,BPSpecOnly] = Utilities.ProcVarargin(varargin,'BPSpecOnly');
[~,PSpecOnly] = Utilities.ProcVarargin(varargin,'PSpecOnly');
[~,All] = Utilities.ProcVarargin(varargin,'All');
[~,sigUnits] = Utilities.ProcVarargin(varargin,'sigUnits');
[~,ConditionFieldLabel] = Utilities.ProcVarargin(varargin,'ConditionFieldLabel','Condition');
[~,ObsPerCond] = Utilities.ProcVarargin(varargin,'ObsPerCond',10);
[~,FileTag] = Utilities.ProcVarargin(varargin,'FileTag','');
[~,LabelsAbbrev] = Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[~,saveFig] = Utilities.ProcVarargin(varargin,'saveFig');

withShuffle =0;

StyleArgs = {'FontSize',9,'FontName','Arial','FontWeight','bold'};

NConds = length(DecodeLabels);
NShuffles = 1;
% NShuffles = 10;

%% Only include units with tunign to at least one of the relevant variables
Dates = Analyze.returnUniqueFieldValues(AllTrialData,'Date');
if BPSpecOnly
    %     [Units,UnitLabels] = Analyze.FaceScratch.FindInterestingUnits(Tag);
    %     Keep = strcmp(UnitLabels,'BodyPart') | strcmp(UnitLabels,'Both');
    ASF = Analyze.FaceScratch3.FindSignificantUnits(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,...
        'BPSpecOnly','Dates',Dates);
elseif PSpecOnly
    ASF = Analyze.FaceScratch3.FindSignificantUnits(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,...
        'PSpecOnly','Dates',Dates);
elseif sigUnits
    ASF = Analyze.FaceScratch3.FindSignificantUnits(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,'Dates',Dates);
else
    ASF = Analyze.FaceScratch3.FindSignificantUnits(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,...
        'All','Dates',Dates);
end


opts.NumShuffles=200;
opts.cvOptions.ValidationType='ClassicCrossValidation';
opts.cvOptions.NReps=30;
opts.cvOptions.NFolds=10;
opts.dec=Predictor.FWClassifier(@Analyze.FaceScratch3.BasicClassifier);

%% first do within and across format classification of Person

Conditions = {{'NancyCheek','TysonCheek'},{'NancyShoulder','TysonShoulder'}};
shufcvAccuracy=[];cvAccuracy=[];
for dateIdx=1:length(Dates)
    
    cDate = Dates{dateIdx};
    dateCode = Blackrock.Helper.date2unitId(cDate);
    cASF = ASF(ASF(:,4)==dateCode,:);
    
    for i=1:length(Conditions)
        cCond=Conditions{i};
        cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
        
        L{i}=Analyze.returnFieldValues(cTrialData,'Person');
        FR{i}=squeeze(Analyze.getNeuralData(cTrialData,cASF,TimeWindow));
        
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
            [FRTrue{i},LabelTrue{i}]=DataConvert.cell2list(FRCell{i});
        end
    end
    %%
    Shuffles = perms(1:2);
    %%
    % rng(123);
    for i=1:length(Conditions)
        %     opts.cvOptions.ValidationType='ResamplingCrossValidation';
        %     [ShuffleRank,AllInfo]=opts.dec.ShuffleTest('TrainingData',...
        %         {FRCell{i},[1:5],LabelCell{i}},'CrossValidate','cvOptions',opts.cvOptions);
        ShuffleRank=1;
        opts.dec.Train('TrainingData',...
            {FRCell{i},[1:2],LabelCell{i}},'CrossValidate','cvOptions',opts.cvOptions);
        
        for j=1:2
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
    Results{dateIdx} = tmp;
    shufcvAccuracy = cat(3,shufcvAccuracy,tmp.shufcvAccuracy);
    cvAccuracy = cat(3,cvAccuracy,tmp.cvAccuracy);
    
end

%%
% next do within and across classification of body part
Conditions = {{'NancyCheek','NancyShoulder'},{'TysonCheek','TysonShoulder'}};

shufcvAccuracy2=[];cvAccuracy2=[];
for dateIdx=1:length(Dates)
    
    cDate = Dates{dateIdx};
    dateCode = Blackrock.Helper.date2unitId(cDate);
    cASF = ASF(ASF(:,4)==dateCode,:);
    
    for i=1:length(Conditions)
        cCond=Conditions{i};
        cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
        
        L{i}=Analyze.returnFieldValues(cTrialData,'BodyPart');
        FR{i}=squeeze(Analyze.getNeuralData(cTrialData,cASF,TimeWindow));
        
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
            [FRTrue{i},LabelTrue{i}]=DataConvert.cell2list(FRCell{i});
        end
    end
    %%
    Shuffles = perms(1:2);
    %%
    % rng(123);
    for i=1:length(Conditions)
        %     opts.cvOptions.ValidationType='ResamplingCrossValidation';
        %     [ShuffleRank,AllInfo]=opts.dec.ShuffleTest('TrainingData',...
        %         {FRCell{i},[1:5],LabelCell{i}},'CrossValidate','cvOptions',opts.cvOptions);
        ShuffleRank=1;
        opts.dec.Train('TrainingData',...
            {FRCell{i},[1:2],LabelCell{i}},'CrossValidate','cvOptions',opts.cvOptions);
        
        for j=1:2
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
    Results2{dateIdx} = tmp;
    shufcvAccuracy2 = cat(3,shufcvAccuracy2,tmp.shufcvAccuracy);
    cvAccuracy2 = cat(3,cvAccuracy2,tmp.cvAccuracy);
    
end

%%
%% Plot the data
% first plot within and across classification of person

IDX=[1 3 2 4];
clear plotData
for i=1:length(IDX)
    [I,J] = ind2sub([2 2],IDX(i));
    plotData(:,i)= squeeze(cvAccuracy(I,J,:));
end

StyleArgs={'FontWeight','bold','FontSize',9,'FontName','Helvetica'};

clr=[];
clr([1:2],:) = repmat([49 163 84]./255,2,1);
clr([3:4],:) = repmat([117 107 177]./255,2,1);

plt.fig('units','inches','width',17,'height',5,'font','Arial','fontsize',9);
pnl = panel();  pnl.margin=15; pnl.pack(2,3);
pnl.fontsize=11; pnl.fontname='Arial';
pnl(1,1).select();
xtickangle(45); set(gca, 'fontweight','bold');

% H = notBoxPlotMod(plotData,[1:4 7:10 13:16 19:22],'jitter',0.3,'useCI');
H = notBoxPlot(plotData,[1 2 4 5],'jitter',0.3);

if withShuffle
    hold on;
    IDX=[1 3 2 4];
    clear plotData
    for i=1:length(IDX)
        [I,J] = ind2sub([2 2],IDX(i));
        plotData(:,i)= squeeze(shufcvAccuracy(I,J,:));
    end
    
    clr1=[];
    clr1([1:4 9:12],:) = repmat([200 200 200]./255,8,1);
    clr1([5:8 13:16],:) = repmat([200 200 200]./255,8,1);
    
    H1 = notBoxPlotMod(plotData,[1 2 4 5],'jitter',0.3,'useCI');
    %         H1 = notBoxPlot(plotData,[1 2 4 5],'jitter',0.3);
    
    %         set([H1.data],...
    %             'MarkerFaceColor',[227 74 51]./255*0.35,...
    %             'markerEdgeColor',[227 74 51]./255*0.35,...
    %             'MarkerSize',1.5)
    set([H1.mu],'color',[227 74 51]./255)
    J=clr1;
    for ii=1:length(H1)
        set(H1(ii).sdPtch,'FaceColor',J(ii,:),'Visible','off',...
            'EdgeColor','none')
        set(H1(ii).semPtch,'FaceColor',J(ii,:),...
            'EdgeColor','none');
    end
end

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
cLabels = {'Ch->Ch','Ch->Sh','Sh->Ch','Sh->Sh'};
xticklabels(cLabels);
set(ax1,'box','off',StyleArgs{:});

set(gca,'TickDir','out');
set(ax1,'TickLength',[0.00025 0.0]);

% line([0 max(H(end).sdPtch.Vertices(:,1))+0.25],[100 100],'Color',[115 115 115]./255,'LineStyle','-.');
anArrow = annotation('arrow', [0 1],[1 1],'Color',[50 50 50]./255,'LineStyle','-.','HeadStyle','cback3','HeadLength',5,'HeadWidth',5);
anArrow.Parent = gca;
anArrow.Position = [0 100 max(H(end).sdPtch.Vertices(:,1))+0.25, 0];
text(max(H(end).sdPtch.Vertices(:,1))+0.45,100,'Ch,Sh','HorizontalAlignment','left','VerticalAlignment','middle',...
    StyleArgs{:},'Color',[50 50 50]./255);


line([0 max(H(end).sdPtch.Vertices(:,1))+0.25],[50 50],'Color',[115 115 115]./255,'LineStyle','--');
text(max(H(end).sdPtch.Vertices(:,1))+0.45,50,'Chance','HorizontalAlignment','left','VerticalAlignment','middle',...
    StyleArgs{:},'Color',[115 115 115]./255);

ylabel('Decode Accuracy (%)')
% xlabel('Digit')
title('Classifying Person',StyleArgs{:})
NC = 22; xlim([.25 NC+.25]);
ylim([0 105])


%% next plot cross-classification of BodyPart

IDX=[1 3 2 4];
clear plotDataBP
for i=1:length(IDX)
    [I,J] = ind2sub([2 2],IDX(i));
    plotDataBP(:,i)= squeeze(cvAccuracy2(I,J,:));
end

clr=[];
clr([1:2],:) = repmat([237 0 38]./255,2,1);
clr([3:4],:) = repmat([31 120 180]./255,2,1);

pnl(1,2).select();
xtickangle(45); set(gca, 'fontweight','bold');

% H = notBoxPlotMod(plotDataBP,[1:4 7:10 13:16 19:22],'jitter',0.3,'useCI');
H = notBoxPlot(plotDataBP,[1 2 4 5],'jitter',0.3);

if withShuffle
    hold on;
    IDX=[1 3 2 4];
    clear plotDataBP
    for i=1:length(IDX)
        [I,J] = ind2sub([2 2],IDX(i));
        plotDataBP(:,i)= squeeze(shufcvAccuracy2(I,J,:));
    end
    
    clr1=[];
    clr1([1:4 9:12],:) = repmat([25 25 25]./255,8,1);
    clr1([5:8 13:16],:) = repmat([25 25 25]./255,8,1);
    
    H1 = notBoxPlotMod(plotDataBP,[1 2 4 5],'jitter',0.3,'useCI');
    
    set([H1.mu],'color',[227 74 51]./255)
    J=clr1;
    for ii=1:length(H)
        set(H1(ii).sdPtch,'FaceColor',J(ii,:)*0.3,...
            'EdgeColor','none')
        set(H1(ii).semPtch,'FaceColor',J(ii,:)*0.3,...
            'EdgeColor','none');
    end
end

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
cLabels = {'Felt->Felt','Felt->Obs','Obs->Felt','Obs->Obs'};
xticklabels(cLabels);
set(ax1,'box','off',StyleArgs{:});

set(gca,'TickDir','out');
set(ax1,'TickLength',[0.00025 0.0]);

% line([0 max(H(end).sdPtch.Vertices(:,1))+0.25],[100 100],'Color',[115 115 115]./255,'LineStyle','-.');
anArrow = annotation('arrow', [0 1],[1 1],'Color',[50 50 50]./255,'LineStyle','-.','HeadStyle','cback3','HeadLength',5,'HeadWidth',5);
anArrow.Parent = gca;
anArrow.Position = [0 100 max(H(end).sdPtch.Vertices(:,1))+0.25, 0];
text(max(H(end).sdPtch.Vertices(:,1))+0.45,100,'Felt,Obs','HorizontalAlignment','left','VerticalAlignment','middle',...
    StyleArgs{:},'Color',[50 50 50]./255);


line([0 max(H(end).sdPtch.Vertices(:,1))+0.25],[50 50],'Color',[115 115 115]./255,'LineStyle','--');
text(max(H(end).sdPtch.Vertices(:,1))+0.45,50,'Chance','HorizontalAlignment','left','VerticalAlignment','middle',...
    StyleArgs{:},'Color',[115 115 115]./255);

ylabel('Decode Accuracy (%)')
% xlabel('Digit')
title('Classifying BodyPart',StyleArgs{:})
NC = 22; xlim([.25 NC+.25]);
ylim([0 105])

% tr/ts-Score are matrices with each row being the dimension
%% Plot
% plt.fig('units','inches','width',6,'height',4,'font','Arial','fontsize',14);
% clf; pnl = panel();  pnl.margin=20; pnl.pack(1,2);
%
% clr=[];
% clr([1],:) = repmat([237 0 38]./255,1,1);
% clr([2],:) = repmat([31 120 180]./255,1,1);
% clr([3],:) = repmat([49 163 84]./255,1,1);
% clr([4],:) = repmat([117 107 177]./255,1,1);
%
% for cDimIdx = 1:2 % Dim to split by
%     pnl(1,cDimIdx).select();
%
%     if cDimIdx==1
%         plotData = mean(cvAccuracy,3);
%     elseif cDimIdx==2
%         plotData = mean(cvAccuracy2,3);
%     end
%
% %     plotData = mean(scoreMat{cDimIdx},3)*100;
% %     signifData = nanmean(signifMat{cDimIdx},3);
%     %     plotData = [trScore(cDimIdx,1) tsScore(cDimIdx,2); tsScore(cDimIdx,1) trScore(cDimIdx,2)];
%     %     signifData = [signifTr(cDimIdx,1) signifTs(cDimIdx,2); signifTs(cDimIdx,1) signifTr(cDimIdx,2)];
%     hh=bar(plotData','FaceColor','flat','ShowBaseLine','off');
%     if cDimIdx==1
%         hh(1).CData(2,:) = clr(4,:);
%         hh(2).CData(2,:) = clr(4,:);
%
%         hh(1).CData(1,:) = clr(3,:);
%         hh(2).CData(1,:) = clr(3,:);
%         cLabels = {'Ch->Ch','Ch->Sh','Sh->Ch','Sh->Sh'};
%         titleLabel = 'Person';
%     elseif cDimIdx==2
%         hh(1).CData(2,:) = clr(2,:);
%         hh(2).CData(2,:) = clr(2,:);
%
%         hh(1).CData(1,:) = clr(1,:);
%         hh(2).CData(1,:) = clr(1,:);
%         cLabels = {'Felt->Felt','Felt->Obs','Obs->Felt','Obs->Obs'};
%         titleLabel = 'BodyPart';
%     end
%     %     bar([.1 .2; .3 .4]);
%     title(sprintf('Classifying %s',titleLabel),StyleArgs{:});
%     ylim([0 110]);
%     xticks([0.875 1.125 1.875 2.125])
%     xticklabels(cLabels);
%     xtickangle(45);
%     %     xlabel('Tested on',StyleArgs{:});
%     ylabel('Decode Accuracy (%)',StyleArgs{:});
%     set(gca,StyleArgs{:});
%     %     legLabels = {};
%     %     for i = 1:length(SplitValNames{cDimIdx})
%     %         legLabels{i} = sprintf('Trained on %s',LabelsAbbrev{-cDimIdx+3}{i});
%     %     end
%     %     leg=legend(legLabels,'Location','SouthOutside');legend boxoff;
%     %     set(leg,StyleArgs{:});
%
%     ax1=gca;
%     plt.hline(100,{'Color','r','LineStyle','--'});
%     plt.hline((1/2)*100,{'Color',[115 115 115]./255,'LineStyle','--'});
%
%     %     axis tight
%     set(ax1,'box','off')
%     ax1.XRuler.Axle.LineStyle = 'none';
%     set(gca,'TickDir','out');
%     set(ax1,'TickLength',[0.00025 0.00]);
%     ylim([0 110])
%
%     % Draw stars for significance
% %     for j = 1:size(signifData,1)
% %         for k = 1:size(signifData,2)
% %             if signifData(j,k) | cDimIdx==2
% %                 x = hh(k).XData(j) + hh(k).XOffset-.05;
% %                 text(x,105,'*',StyleArgs{:});
% %             end
% %         end
% %     end
% %     colormap(gca,clr);
% end

%%
DimNames = {'BodyPart','Person'};
filename = sprintf('CrossDecode-%s-vs-%s',DimNames{1},DimNames{2});
if BPSpecOnly
    var = 'BPSpecOnly';
    filename = [filename '-BPSpecOnly'];
    suptitle(sprintf('BPSpec units only (%d)',size(ASF,1)));
elseif PSpecOnly
    var = 'PSpecOnly';
    filename = [filename '-PSpecOnly'];
    suptitle(sprintf('PSpec units only (%d)',size(ASF,1)));
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


