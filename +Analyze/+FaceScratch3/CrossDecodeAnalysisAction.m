function CrossDecodeAnalysisAction(AllTrialData,Labels,DecodeLabels,TimeWindow,OutDir,Tag,varargin)
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
[~,ConditionFieldLabel] = Utilities.ProcVarargin(varargin,'ConditionFieldLabel','Condition');
[~,ObsPerCond] = Utilities.ProcVarargin(varargin,'ObsPerCond',10);
[~,FileTag] = Utilities.ProcVarargin(varargin,'FileTag','');

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
else
    ASF = Analyze.FaceScratch3.FindSignificantUnits(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,'Dates',Dates);
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
        [FR2,Labels2] = Analyze.SampleFeaturePopulation(FR(trIdx,:),Labels(trIdx),'Type','Null','NumObservationsPerCondition',ObsPerCond);
        mdl = fitcdiscr(FR2,Labels2,'discrimType','diagLinear');
%         mdl = fitcdiscr(FR2,Labels2,'discrimType','diagquadratic');
%         mdl = fitcdiscr(FR2,Labels2,'discrimType','linear');
        
        % Cross val perf
        partitionedModel = crossval(mdl, 'KFold',ObsPerCond);
        cvScore = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
        perfDist = ShuffledPerformance(NShuffles,FR2,Labels2);
        signifTr(cDimIdx,cValIdx) = (nnz(cvScore > perfDist) > .95*length(perfDist));
        scoreMat{cDimIdx}(cValIdx,cValIdx) = cvScore;
        signifMat{cDimIdx}(cValIdx,cValIdx) = (nnz(cvScore > perfDist) > .95*length(perfDist));
%         disp(scoreMat{cDimIdx})
        
        % Test on other values
        testDimVals = unique(CondLabels{cDimIdx});
        testDimVals(testDimVals==cValIdx) = [];

        for tt = 1:length(testDimVals)
            tsDimVal = testDimVals(tt);
            fprintf('\t\t Testing on %s...\n',SplitValNames{cDimIdx}{tsDimVal});
            tsIdx = CondLabels{cDimIdx} == tsDimVal;
        
            % Test perf on other value of dimension
            [FR3,Labels3] = Analyze.SampleFeaturePopulation(FR(tsIdx,:),Labels(tsIdx),'Type','Null','NumObservationsPerCondition',ObsPerCond);
            predLabels = predict(mdl,FR3);
            cxScore = sum(predLabels == Labels3)/length(Labels3);
            perfDist = ShuffledPerformance(NShuffles,FR3,Labels3);
            signifTs(cDimIdx,cValIdx) = (nnz(cxScore > perfDist) > .95*length(perfDist));
            
            scoreMat{cDimIdx}(cValIdx,tsDimVal) = cxScore;
            signifMat{cDimIdx}(cValIdx,tsDimVal) = (nnz(cxScore > perfDist) > .95*length(perfDist));
%             disp(scoreMat{cDimIdx})
        end
        % Save scores
        trScore{cDimIdx,cValIdx} = cvScore;
        tsScore{cDimIdx,cValIdx} = cxScore;
    end
end

% tr/ts-Score are matrices with each row being the dimension 
%% Plot
plt.fig('units','inches','width',8,'height',5,'font','Arial','fontsize',9);
clf; pnl = panel();  pnl.margin=20; pnl.pack(1,2);

for cDimIdx = 1:NSplitDims % Dim to split by
    cSplitDim = DimNames{cDimIdx};
    cCatDim = DimNames{-cDimIdx+3};
    NCatVals = length(unique(CondLabels{-cDimIdx+3}));
%     cLabelDim = Dims{cOtherIdx}; % other dim is dim to categorize
    pnl(1,cDimIdx).select();
    plotData = scoreMat{cDimIdx};
    signifData = signifMat{cDimIdx};
%     plotData = [trScore(cDimIdx,1) tsScore(cDimIdx,2); tsScore(cDimIdx,1) trScore(cDimIdx,2)];
%     signifData = [signifTr(cDimIdx,1) signifTs(cDimIdx,2); signifTs(cDimIdx,1) signifTr(cDimIdx,2)];
    hh=bar(plotData);
%     bar([.1 .2; .3 .4]);
    title(sprintf('Decoding %s split by %s',cCatDim, cSplitDim),StyleArgs{:});
    ylim([0 1.2]);
    xlabel('Tested on:',StyleArgs{:});
    ylabel('Performance',StyleArgs{:});
    set(gca,'XTickLabel',SplitValNames{cDimIdx},StyleArgs{:});
    legLabels = {};
    for i = 1:length(SplitValNames{cDimIdx})
        legLabels{i} = sprintf('Trained on %s',SplitValNames{cDimIdx}{i});
    end
    leg=legend(legLabels,'Location','SouthEast');
    set(leg,StyleArgs{:});
    
    plt.hline(1,'g');
    plt.hline(1/NCatVals,'r');
    
    % Draw stars for significance
    for j = 1:size(signifData,1)
        for k = 1:size(signifData,2)
            if signifData(j,k)
                x = hh(k).XData(j) + hh(k).XOffset-.05;
                text(x,1.05,'*',StyleArgs{:});
            end
        end
    end
    colormap(gca,lines(size(scoreMat{cDimIdx},2)));
end

%%
filename = sprintf('Model-CrossDecode-%s-vs-%s',DimNames{1},DimNames{2});
if BPSpecOnly
    filename = [filename '-BPSpecOnly'];
    suptitle(sprintf('BPSpec units only (%d)',size(ASF,1)));
elseif PSpecOnly
    filename = [filename '-PSpecOnly'];
    suptitle(sprintf('PSpec units only (%d)',size(ASF,1)));
else
    suptitle(sprintf('All significant units (%d)',size(ASF,1)));
end

if ~isempty(FileTag)
    filename = [filename '-' FileTag];
end
filename = [filename '-allunits'];
plt.SaveFigure(2,fullfile(OutDir,'PopData','Models'),filename,'PNG','SVGI')


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
