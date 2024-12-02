function CrossDecodeAnalysisTrial(AllTrialData,Labels,DecodeLabels,TimeWindow,OutDir,Tag,varargin)
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
[~,All] = Utilities.ProcVarargin(varargin,'All');
[~,PSpecOnly] = Utilities.ProcVarargin(varargin,'PSpecOnly');
[~,sigUnits] = Utilities.ProcVarargin(varargin,'sigUnits');
[~,ConditionFieldLabel] = Utilities.ProcVarargin(varargin,'ConditionFieldLabel','Condition');
[~,ObsPerCond] = Utilities.ProcVarargin(varargin,'ObsPerCond',10);
[~,FileTag] = Utilities.ProcVarargin(varargin,'FileTag','');
[~,LabelsAbbrev] = Utilities.ProcVarargin(varargin,'LabelsAbbrev','');
[~,TaskLabel] = Utilities.ProcVarargin(varargin,'TaskLabel',{});
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
    ASF = Analyze.FaceScratch3.FindSignificantUnits(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,...
        'All','Dates',Dates);
end

cTrialData = Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,ConditionFieldLabel,DecodeLabels);
BaseTrialData = Analyze.SubSelectTrials(AllTrialData,'Phase',BasePhase,ConditionFieldLabel,DecodeLabels);
FR = squeeze(Analyze.getNeuralData(cTrialData,ASF,TimeWindow));
FRBase = squeeze(Analyze.getNeuralData(BaseTrialData,ASF,BaseWindow));
FR = FR - FRBase;
% FR = FR - nanmean(FRBase);
Condition = Analyze.returnFieldValues(cTrialData,ConditionFieldLabel);

for i = 1:size(Dims,2)-1
    CondCode = zeros(length(Condition),1);
    for j = 1:size(Dims,1)
        CondCode(strcmp(Condition,Dims{j,1})) = Dims{j,i+1};
    end
    CondLabels{i} = CondCode;
end
NSplitDims = size(Dims,2)-1;


%%
for cDimIdx=1:2
    cOtherIdx = -cDimIdx+3;
    
    nSplits = length(unique(CondLabels{cDimIdx}));
    countRow=[1:nSplits:NConds];
    
    NDimVals = length(unique(CondLabels{cDimIdx})); % N vals dim can take
    oNDimVals = length(unique(CondLabels{cOtherIdx}));
    fprintf('Splitting by %s...\n',DimNames{cDimIdx});
    scoreMat{cDimIdx} = nan(NDimVals*oNDimVals,NDimVals*oNDimVals);
    signifMat{cDimIdx} = nan(NDimVals*oNDimVals,NDimVals*oNDimVals);
    for cValIdx = 1:NDimVals
        countCol=[1:nSplits:NConds];
        ccountRow = cValIdx-1+countRow;
        
        %         cVal = cDim{cValIdx};
        trIdxObj = CondLabels{cDimIdx} == cValIdx;
        
        stimLabels = CondLabels{cOtherIdx};
        
        for oValIdx=1:length(unique(stimLabels))
            
            cLabels = stimLabels==oValIdx;
            
            %             trIdxStim= stimLabels==oValIdx;
            %             trIdx = [trIdxObj & trIdxStim];
            
            fprintf('\t Training on %s...%s...\n',SplitValNames{cDimIdx}{cValIdx}, SplitValNames{cOtherIdx}{oValIdx});
            % Train model
            
            [FR2,Labels2] = Analyze.SampleFeaturePopulation(FR(trIdxObj,:),cLabels(trIdxObj),'Type','Null','NumObservationsPerCondition',ObsPerCond);
            mdl = fitcdiscr(FR2,Labels2,'discrimType','diagLinear');
            %                     mdl = fitcdiscr(FR2,Labels2,'discrimType','diagquadratic');
            %         mdl = fitcdiscr(FR2,Labels2,'discrimType','linear');
            
            % Cross val perf
            fprintf('\t CrossValidating %s...%s...\n',SplitValNames{cDimIdx}{cValIdx}, SplitValNames{cOtherIdx}{oValIdx});
            
            partitionedModel = crossval(mdl, 'KFold',ObsPerCond);
            cvScore = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
            perfDist = ShuffledPerformance(NShuffles,FR2,Labels2);
            %         signifTr(cDimIdx,cValIdx) = (nnz(cvScore > perfDist) > .95*length(perfDist));
            scoreMat{cDimIdx}(ccountRow(oValIdx),ccountRow(oValIdx)) = cvScore;
            signifMat{cDimIdx}(ccountRow(oValIdx),ccountRow(oValIdx)) = (nnz(cvScore > perfDist) > .95*length(perfDist));
            
            otestDimVals = unique(stimLabels);
            otestDimVals(otestDimVals==oValIdx) =[];
            
            for ott=1:length(otestDimVals)
                otsDimVal = otestDimVals(ott);
                ccountCol=cValIdx-1+countCol;
                fprintf('\t\t\t Testing on %s...%s...\n',SplitValNames{cDimIdx}{cValIdx},SplitValNames{cOtherIdx}{otsDimVal});
                otsLabels = stimLabels==otsDimVal;
                [FR4,Labels4] = Analyze.SampleFeaturePopulation(FR(trIdxObj,:),otsLabels(trIdxObj),'Type','Null','NumObservationsPerCondition',ObsPerCond);
                %test on other values of other dimension
                predLabels = predict(mdl,FR4);
                cxScore = sum(predLabels == Labels4)/length(Labels4);
                perfDist = ShuffledPerformance(NShuffles,FR4,Labels4);
                %                 signifTs(cDimIdx,cValIdx) = (nnz(cxScore > perfDist) > .95*length(perfDist));
                
                scoreMat{cDimIdx}(ccountRow(oValIdx),ccountCol(otsDimVal)) = cxScore;
                signifMat{cDimIdx}(ccountRow(oValIdx),ccountCol(otsDimVal)) = (nnz(cxScore > perfDist) > .95*length(perfDist));
                
            end
            
            % Test on other values
            testDimVals = unique(CondLabels{cDimIdx});
            testDimVals(testDimVals==cValIdx) = [];
            
            for tt = 1:length(testDimVals)
                tsDimVal = testDimVals(tt);
                otestIdxObj = CondLabels{cDimIdx} == tsDimVal;
                ccountCol = tsDimVal-1+countCol;
                for cott=1:length(unique(stimLabels))
                    
                    cottLabels = stimLabels==cott;
                    
                    fprintf('\t\t\t Testing on %s...%s...\n',SplitValNames{cDimIdx}{tsDimVal},SplitValNames{cOtherIdx}{cott});
                    
                    [FR3,Labels3] = Analyze.SampleFeaturePopulation(FR(otestIdxObj,:),cottLabels(otestIdxObj),'Type','Null','NumObservationsPerCondition',ObsPerCond);
                    %test on other values of other dimension
                    predLabels = predict(mdl,FR3);
                    cxScore = sum(predLabels == Labels3)/length(Labels3);
                    perfDist = ShuffledPerformance(NShuffles,FR3,Labels3);
                    %                 signifTs(cDimIdx,cValIdx) = (nnz(cxScore > perfDist) > .95*length(perfDist));
                    
                    scoreMat{cDimIdx}(ccountRow(oValIdx),ccountCol(cott)) = cxScore;
                    signifMat{cDimIdx}(ccountRow(oValIdx),ccountCol(cott)) = (nnz(cxScore > perfDist) > .95*length(perfDist));
                    
                end
                %             disp(scoreMat{cDimIdx})
            end
        end
    end
end

%% Plot
Fig1 = plt.fig('units','inches','width',14,'height',5,'font','Arial','fontsize',9);
clf; pnl = panel();  pnl.margin=15; pnl.pack(1,2);orient(Fig1,'landscape');
pnl(1,1).pack('v',{3/4,1/4},'h',{1/4,3/4});
pnl(1,2).pack('v',{3/4,1/4},'h',{3/4,1/4});

cmap = cbrewer('seq','Greys',256);

NC=16;

for cDimIdx = 1:2 % Dim to split by
    nSplits = length(unique(CondLabels{cDimIdx}));
    
    cSplitDim = DimNames{cDimIdx};
    cCatDim = DimNames{-cDimIdx+3};
    NCatVals = length(unique(CondLabels{-cDimIdx+3}));
    %     cLabelDim = Dims{cOtherIdx}; % other dim is dim to categorize
    if cDimIdx==1
        pnl(1,cDimIdx,1,2).select();
    else
        pnl(1,cDimIdx,1,1).select();
    end
    
    ax1=gca; set(ax1,'box','off','FontWeight','bold','FontSize',9);
    plotData = scoreMat{cDimIdx};
    signifData = signifMat{cDimIdx};
    %     plotData = [trScore(cDimIdx,1) tsScore(cDimIdx,2); tsScore(cDimIdx,1) trScore(cDimIdx,2)];
    %     signifData = [signifTr(cDimIdx,1) signifTs(cDimIdx,2); signifTs(cDimIdx,1) signifTr(cDimIdx,2)];
    
    hh = imagesc(plotData); alpha(0.7); cc = colorbar; colormap(cmap); axis image;
    axis([0 NConds 0 NConds]+.5);
    xlim([.25 NC+.25]);
    
    title(sprintf('Decoding %s split by %s',cCatDim, cSplitDim),StyleArgs{:});
    
    xlh = xlabel('Predicted Labels');
    ylh = ylabel('True Labels');
    ylh.Position(1) = ylh.Position(1)-1;
    %     xlh.Position(2) = xlh.Position(2);
    xticks([1:NConds]); yticks([1:NConds])
    
    if cDimIdx==2
        % for the second dimension, reorganize the Labels
        %         newLabels = cell(1,length(unique(Labels)));
        %         count=1;
        %         for i=1:length(unique(vertcat(Dims{:,2})))
        %             for j=1:length(unique(vertcat(Dims{:,3})))
        %                 cL = cellfun(@(x)x==i,Dims(:,2),'UniformOutput',true) &...
        %                     cellfun(@(x)x==j,Dims(:,3),'UniformOutput',true);
        %                 newLabels{count} = Dims{cL,1};
        %                 count = count+1;
        %             end
        %         end
        
        newLabels = cell(1,length(LabelsAbbrev));
        for i=1:length(LabelsAbbrev)
            
            cLabels = LabelsAbbrev{i};
            if i==1
                newLabels{i} = unique(cLabels,'stable');
            elseif i==2
                newLabels{i}= repmat(unique(cLabels,'stable'),1,length(unique(LabelsAbbrev{1})));
            end
        end
        % use newLabels in the figure
        set(gca,'XTickLabel',newLabels{2},StyleArgs{:});
        set(gca,'YTickLabel',newLabels{2},StyleArgs{:});
        set(gca,'YDir','normal')
        %         xlim([.25 NC+.25]);
        xtickangle(45); set(gca, 'fontweight','bold','fontsize',9');
        
        StimNames =newLabels{1};
        count=1;
        for i=[nSplits/2:nSplits:NConds]
            text(i+0.5,-1,StimNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
                'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255);
            text(-1,i+0.5,StimNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
                'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255,'Rotation',90);
            count = count+1;
        end
        axis tight;
        
    else
        set(gca,'XTick',1:NConds,'XTickLabel',LabelsAbbrev{1});
        set(gca,'YTick',1:NConds,'YTickLabel',LabelsAbbrev{1});
        set(gca,'YDir','normal')
        %         xlim([.25 NC+.25]);
        xtickangle(45); set(gca, 'fontweight','bold','fontsize',9');
        
        StimNames =LabelsAbbrev{2};
        count=1;
        for i=[nSplits/2:nSplits:NConds]
            text(i+0.5,-1,StimNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
                'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255);
            text(-1,i+0.5,StimNames{count},'HorizontalAlignment','center','VerticalAlignment','bottom',...
                'FontWeight','bold','FontSize',9,'FontName','Arial','Color',[99 99 99]./255,'Rotation',90);
            count = count+1;
        end
        axis tight;
        
    end
    xtickangle(45);
    
    ax1=gca;
    vline([0.5:nSplits:NConds+0.5],'k-')
    hline([0.5:nSplits:NConds+0.5],'k-')
    set(ax1,'box','off','FontWeight','bold','FontSize',9,'FontName','Arial');
    ax1.XRuler.Axle.LineStyle = 'none';
    set(gca,'TickDir','out');
    set(ax1,'TickLength',[0.00 0.0]);
    
end


%%
filename = sprintf('CrossClassification-%s-vs-%s',DimNames{1},DimNames{2});
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

FigsDir = fullfile(env.get('results'),'FaceScratch3','SUAnal',[TaskLabel '-Go'],'PopData','Mine');

filename = [filename var];
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
