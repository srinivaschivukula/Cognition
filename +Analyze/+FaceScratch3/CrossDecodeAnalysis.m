function CrossDecodeAnalysis(AllTrialData,Labels,DecodeLabels,TimeWindow,OutDir,Tag,varargin)
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
[~,sigUnits] = Utilities.ProcVarargin(varargin,'sigUnits');
[~,ConditionFieldLabel] = Utilities.ProcVarargin(varargin,'ConditionFieldLabel','Condition');
[~,ObsPerCond] = Utilities.ProcVarargin(varargin,'ObsPerCond',10);
[~,FileTag] = Utilities.ProcVarargin(varargin,'FileTag','');
[~,LabelsAbbrev] = Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[~,saveFig] = Utilities.ProcVarargin(varargin,'saveFig');

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
    ASF = Analyze.returnUniqueFieldValues(AllTrialData,'UnitIDsSorted');
end

cTrialData = Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,ConditionFieldLabel,DecodeLabels);
BaseTrialData = Analyze.SubSelectTrials(AllTrialData,'Phase',BasePhase,ConditionFieldLabel,DecodeLabels);
FR = squeeze(Analyze.getNeuralData(cTrialData,ASF,TimeWindow));
FRBase = squeeze(Analyze.getNeuralData(BaseTrialData,ASF,BaseWindow));
FR = FR - FRBase;
% FR = FR - nanmean(FRBase);
Condition = Analyze.returnFieldValues(cTrialData,ConditionFieldLabel);

% Balance so equal number of each
% minTrials = 100;
% for c = 1:length(DecodeLabels)
%     condidx{c} = find(strcmp(Condition,DecodeLabels{c}));
%     minTrials = min(minTrials,length(condidx{c}));
% end
%
% FR2 = [];
% Condition2 = [];
% for c = 1:length(DecodeLabels)
%     FR2 = [FR2; FR(condidx{c}(1:minTrials),:)];
%     Condition2 = [Condition2; Condition(condidx{c}(1:minTrials))];
% end
% FR = FR2;
% Condition = Condition2;

% switch Tag
%     case 'XBodySpec'
%         models2use{g} = [1:4];
%         compareIdx = 1;
%     case 'XFixation'
%         models2use{g} = [5 7];
%         compareIdx = 2;
%     case 'XVideoLiveXBodySpecC'
%         models2use{g} = [3 5];
%         compareIdx = 2;
%     case 'XActSense'
%         models2use{g} = 1:2;
%         compareIdx = 2;
% %     models2use{g} = 1:size(BIC,2);
% end

if strcmp(Tag,'XBodySpec') || strcmp(Tag,'XBodySpecC')
    BodyPart = {'Cheek','Shoulder'};
    Person = {'Nancy','Tyson'};
    Dims = {BodyPart,Person};
    DimNames = {'BodyPart','Person'};
    
    for i = 1:length(Dims)
        CondPro = ~cellfun(@isempty,strfind(Condition,Dims{i}{1}));
        CondAnti = ~CondPro;
        CondCode = CondPro + 2*CondAnti;
        CondLabels{i} = CondCode;
    end
    NSplitDims = 2;
    SplitValNames = Dims;
else
    for i = 1:size(Dims,2)-1
        CondCode = zeros(length(Condition),1);
        for j = 1:size(Dims,1)
            CondCode(strcmp(Condition,Dims{j,1})) = Dims{j,i+1};
        end
        CondLabels{i} = CondCode;
    end
    NSplitDims = size(Dims,2)-1;
end

% ChCond = ~cellfun(@isempty,strfind(Condition,'Cheek'));
% ShCond = ~ChCond;
% NaCond = ~cellfun(@isempty,strfind(Condition,'Nancy'));
% TyCond = ~NaCond;
% CondBP = ChCond + 2*ShCond;
% CondPerson = NaCond + 2*TyCond;
% CondLabels = {CondBP, CondPerson};

%%
FRtrain = FR;
FRtest = FR;
for cDimIdx = 1:NSplitDims % Dim to split by
    cOtherIdx = -cDimIdx+3;
    %     cDim = Dims{cDimIdx};
    %     cLabelDim = Dims{cOtherIdx}; % other dim is dim to categorize
    NDimVals = length(unique(CondLabels{cDimIdx})); % N vals dim can take
    fprintf('Splitting by %s...\n',DimNames{cDimIdx});
    scoreMat{cDimIdx} = nan(NDimVals,NDimVals);
    signifMat{cDimIdx} = nan(NDimVals,NDimVals);
    for cValIdx = 1:NDimVals
        %         cVal = cDim{cValIdx};
        trIdx = CondLabels{cDimIdx} == cValIdx;
        
        fprintf('\t Training on %s...\n',SplitValNames{cDimIdx}{cValIdx});
        % Train model
        Labels = CondLabels{cOtherIdx};
        [FR2,Labels2] = Analyze.SampleFeaturePopulation(FRtrain(trIdx,:),Labels(trIdx),'Type','Null','NumObservationsPerCondition',ObsPerCond);
        
        for kk=1:ObsPerCond
            testIdx = kk:(ObsPerCond-4):size(FR2,1);
            trainIdx = setdiff(1:size(FR2,1),testIdx);
            
            mdl = fitcdiscr(FR2(trainIdx,:),Labels2(trainIdx),'discrimType','diagLinear');
            %         mdl = fitcdiscr(FR2,Labels2,'discrimType','diagquadratic');
            %         mdl = fitcdiscr(FR2,Labels2,'discrimType','linear');
            
            % Cross val perf
            partitionedModel = crossval(mdl, 'KFold',ObsPerCond);
            cvScore = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
            perfDist = ShuffledPerformance(NShuffles,FR2,Labels2);
            signifTr(cDimIdx,cValIdx,kk) = (nnz(cvScore > perfDist) > .95*length(perfDist));
            scoreMat{cDimIdx}(cValIdx,cValIdx,kk) = cvScore;
            signifMat{cDimIdx}(cValIdx,cValIdx,kk) = (nnz(cvScore > perfDist) > .95*length(perfDist));
            %         diasp(scoreMat{cDimIdx})
            
            % Test on other values
            testDimVals = unique(CondLabels{cDimIdx});
            testDimVals(testDimVals==cValIdx) = [];
            
            for tt = 1:length(testDimVals)
                tsDimVal = testDimVals(tt);
                fprintf('\t\t Testing on %s...\n',SplitValNames{cDimIdx}{tsDimVal});
                tsIdx = CondLabels{cDimIdx} == tsDimVal;
                
                % Test perf on other value of dimension
                [FR3,Labels3] = Analyze.SampleFeaturePopulation(FRtest(tsIdx,:),Labels(tsIdx),'Type','Null','NumObservationsPerCondition',ObsPerCond);
                predLabels = predict(mdl,FR3(testIdx,:));
                cxScore = sum(predLabels == Labels3(testIdx))/length(Labels3(testIdx));
                perfDist = ShuffledPerformance(NShuffles,FR3,Labels3);
                signifTs(cDimIdx,cValIdx,kk) = (nnz(cxScore > perfDist) > .95*length(perfDist));
                
                scoreMat{cDimIdx}(cValIdx,tsDimVal,kk) = cxScore;
                %                 signifMat{cDimIdx}{kk}(cValIdx,tsDimVal) = (nnz(cxScore > perfDist) > .95*length(perfDist));
                %             disp(scoreMat{cDimIdx})
            end
            % Save scores
            trScore{cDimIdx,cValIdx} = cvScore;
            tsScore{cDimIdx,cValIdx} = cxScore;
        end
    end
end
% tr/ts-Score are matrices with each row being the dimension
%% Plot
plt.fig('units','inches','width',6,'height',4,'font','Arial','fontsize',14);
clf; pnl = panel();  pnl.margin=20; pnl.pack(1,2);

clr=[];
clr([1],:) = repmat([237 0 38]./255,1,1);
clr([2],:) = repmat([31 120 180]./255,1,1);
clr([3],:) = repmat([49 163 84]./255,1,1);
clr([4],:) = repmat([117 107 177]./255,1,1);

for cDimIdx = 1:NSplitDims % Dim to split by
    cSplitDim = DimNames{cDimIdx};
    cCatDim = DimNames{-cDimIdx+3};
    NCatVals = length(unique(CondLabels{-cDimIdx+3}));
    %     cLabelDim = Dims{cOtherIdx}; % other dim is dim to categorize
    pnl(1,cDimIdx).select();
    plotData = mean(scoreMat{cDimIdx},3)*100;
    signifData = nanmean(signifMat{cDimIdx},3);
    %     plotData = [trScore(cDimIdx,1) tsScore(cDimIdx,2); tsScore(cDimIdx,1) trScore(cDimIdx,2)];
    %     signifData = [signifTr(cDimIdx,1) signifTs(cDimIdx,2); signifTs(cDimIdx,1) signifTr(cDimIdx,2)];
    hh=bar(plotData','FaceColor','flat','ShowBaseLine','off');
    if cDimIdx==1
        hh(1).CData(2,:) = clr(4,:);
        hh(2).CData(2,:) = clr(4,:);
        
        hh(1).CData(1,:) = clr(3,:);
        hh(2).CData(1,:) = clr(3,:);
        cLabels = {'Ch->Ch','Ch->Sh','Sh->Ch','Sh->Sh'};
    elseif cDimIdx==2
        hh(1).CData(2,:) = clr(2,:);
        hh(2).CData(2,:) = clr(2,:);
        
        hh(1).CData(1,:) = clr(1,:);
        hh(2).CData(1,:) = clr(1,:);
         cLabels = {'Felt->Felt','Felt->Obs','Obs->Felt','Obs->Obs'};
    end
    %     bar([.1 .2; .3 .4]);
    title(sprintf('Classifying %s',cCatDim),StyleArgs{:});
    ylim([0 110]);
    xticks([0.875 1.125 1.875 2.125])
    xticklabels(cLabels);
    xtickangle(45);
    %     xlabel('Tested on',StyleArgs{:});
    ylabel('Performance',StyleArgs{:});
    set(gca,StyleArgs{:});
%     legLabels = {};
%     for i = 1:length(SplitValNames{cDimIdx})
%         legLabels{i} = sprintf('Trained on %s',LabelsAbbrev{-cDimIdx+3}{i});
%     end
%     leg=legend(legLabels,'Location','SouthOutside');legend boxoff;
%     set(leg,StyleArgs{:});
    
    ax1=gca;
    plt.hline(100,{'Color','r','LineStyle','--'});
    plt.hline((1/NCatVals)*100,{'Color',[115 115 115]./255,'LineStyle','--'});
    
    %     axis tight
    set(ax1,'box','off')
    ax1.XRuler.Axle.LineStyle = 'none';
    set(gca,'TickDir','out');
    set(ax1,'TickLength',[0.00025 0.00]);
    ylim([0 110])
    
    % Draw stars for significance
    for j = 1:size(signifData,1)
        for k = 1:size(signifData,2)
            if signifData(j,k) | cDimIdx==2
                x = hh(k).XData(j) + hh(k).XOffset-.05;
                text(x,105,'*',StyleArgs{:});
            end
        end
    end
    colormap(gca,clr);
end

%%
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


function perf = ShuffledPerformance(NShuffles,Data,Labels,mdl)
% Returns distribution of models trained on Data with labels shuffled. If
% model is specified, uses that modeled and tests on Data with labels
% shuffled instead of training a new model.

warning off;
perf = zeros(1,NShuffles);
for i = 1:NShuffles
    % Shuffle Labels
    LabelsShuff = Labels(randperm(length(Labels)));
    
    if nargin < 4 % No model specified, train new one
        mdl = fitcdiscr(Data,LabelsShuff,'discrimType','diagLinear');
        partitionedModel = crossval(mdl, 'KFold',10);
        perf(i) = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
    else
        pred = predict(mdl,Data);
        perf(i) = sum(pred == LabelsShuff)/length(LabelsShuff);
    end
end
warning on;
