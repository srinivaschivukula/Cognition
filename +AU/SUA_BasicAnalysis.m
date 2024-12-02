function [Results]=SUA_BasicAnalysis(AllTrialData,Phase,unit,timeWindow,opts)
%%(AllTrialData,Phase,unit,timeWindow,opts)
% unit=[1 57 1];
% Phase='Go1';
% timeWindow=[.5 1];
% opts.BaselinePhase='CueTarget'

if iscell(opts.BaselinePhase)
cTrialDataBaseLine=Analyze.SubSelectTrials(AllTrialData,opts.BaselinePhase{:});
FRBaseLine=Analyze.getNeuralData(cTrialDataBaseLine, unit,opts.BaselineWindow);    
else
cTrialDataBaseLine=Analyze.SubSelectTrials(AllTrialData,'Phase',opts.BaselinePhase);
FRBaseLine=Analyze.getNeuralData(cTrialDataBaseLine, unit,opts.BaselineWindow);
end

Labels=opts.Labels;

if isfield(opts,'ConditionFieldLabel')
   ConditionFieldLabel=opts.ConditionFieldLabel;
else
     ConditionFieldLabel='Condition';
end


cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase);

ConditionLabel=Analyze.returnFieldValues(cTrialData,ConditionFieldLabel);

if isfield(opts,'CondRemap')
ConditionLabel=feval(opts.CondRemap,ConditionLabel);
end

tmp=cellfun(@(x)strcmp(ConditionLabel,x),Labels,'UniformOutput' ,false);
isValidLabel=any(cat(2,tmp{:}),2);

ConditionLabel=ConditionLabel(isValidLabel);

Condition=DataConvert.labels2IDX(ConditionLabel,Labels);

[ClassINDX,ClassNames,N]=AU.getLabelAssociations(Condition,ConditionLabel);

FR=Analyze.getNeuralData(cTrialData, unit,timeWindow);
FR=FR(isValidLabel);

ConditionOrig=Condition;
FROrig=FR;

[p,tbl,stats]=anova1(FR,Condition,'off');
Results.pAnova=p(1);


%%

FR=[FR;FRBaseLine];
Condition=[Condition;FRBaseLine*0];

data = table(Condition,FR);
data.Condition=nominal(data.Condition);
model= 'FR~Condition';
% % warning off
[tmp1,mdl1]=Analyze.FitGLM(data,model,'CrossValidate');

Results.cvR2 = tmp1.CVR2;
Results.Pvals=mdl1.Coefficients.pValue(2:end);
Results.Coef=mdl1.Coefficients.Estimate(2:end);
Results.tStat=mdl1.Coefficients.tStat(2:end);
tmp=mdl1.coefCI;
Results.CoefCI=tmp(2:end,:);
Results.R2adj=mdl1.Rsquared.Adjusted;
Results.R2=mdl1.Rsquared.Ordinary;
Results.pFull=mdl1.coefTest;
% figure; plt.barpatch([],{Results.Coef,Results.CoefCI'})
%%
% if opts.shuffle; FR=Shuffle(FR); end
% [testTuningRank,trainTuningRank,order,RAW,AUCrank]=AU.ComputeRankOrderList(FROrig-mean(FRBaseLine),ConditionOrig,250,'ComputeAUC');
[testTuningRank,trainTuningRank,order,RAW,AUCrank]=AU.ComputeRankOrderList(FROrig-mean(FRBaseLine),ConditionOrig,250);

Results.testTuningRank=testTuningRank;
Results.trainTuningRank=trainTuningRank;
Results.trainSplit=cell2mat({RAW.train}');
Results.testSplit=cell2mat({RAW.test}');


FRNB=FROrig-mean(FRBaseLine);
out=DataConvert.list2cell(FRNB,ConditionOrig);
FRmu=cellfun(@mean,out);
FRsig=cellfun(@std,out);
FRvar=cellfun(@var,out);
meanSig=mean(FRsig);

DeMean=cellfun(@(x)x-mean(x),out,'UniformOutput',false);
Results.Variance=var(cell2mat(DeMean(:)));
Results.STD=std(cell2mat(DeMean(:)));
%%

 Results.testdPBS=(testTuningRank(1)-testTuningRank(2))/sqrt(Results.Variance);
 Results.testdPBW=(testTuningRank(1)-testTuningRank(end))/sqrt(Results.Variance);
 
 Results.trainPBS=(trainTuningRank(1)-trainTuningRank(2))/sqrt(Results.Variance);
 Results.trainPBW=(trainTuningRank(1)-trainTuningRank(end))/sqrt(Results.Variance);
 
for i=1:length(ClassINDX)
    Results.test_dPBVersus(i)=(testTuningRank(1)-testTuningRank(i))/sqrt(Results.Variance);
    Results.traind_dPBVersus(i)=(trainTuningRank(1)-trainTuningRank(i))/sqrt(Results.Variance);
end

%%
if isfield(opts,'ComputeReliability') && opts.ComputeReliability
RandSplits=252;
[FRmat,cLabel]=DataConvert.list2mat(FROrig-mean(FRBaseLine),ConditionOrig);
[Results.ReliabilityNB,Results.PrefActNB,Results.ActionCountsNB]=Analyze.ActionObs.ComputeReliability(FRmat,[],RandSplits);
[Results.ReliabilityNB_NEG,Results.PrefActNB_NEG,Results.ActionCountsNB_NEG]=Analyze.ActionObs.ComputeReliability(-FRmat,[],RandSplits);
end
%%

FR=FROrig;
Condition=ConditionOrig;
AUCPairs=nan(length(ClassINDX));
for i=1:length(ClassINDX)
    for j=i+1:length(ClassINDX)
    idx1=find(Condition==i);
    idx2=find(Condition==j);
    p(i,j) = ranksum(FR(idx1),FR(idx2));    
    [AUCPairs(i,j),ShuffleRank,NullDist]=Analyze.TwoGroup.ShuffleCompare(FR(idx1),FR(idx2),'Type','AUC','NShuffles',0);
    AUCPairs(j,i)=AUCPairs(i,j);    
    end
end
Results.AUCPairs=AUCPairs;
Results.PairsSig=p;
%%

clear p
for i=1:length(ClassINDX)
    idx=find(Condition==i);
    idxAll=find(Condition~=i);
    
    % conditional on whether baseline is ITI or reference condition.
    if length(FR)==length(FRBaseLine)
        FRBASE=FRBaseLine(idx);
    else
        FRBASE=FRBaseLine;
    end
    
    p(i) = ranksum(FR(idx),FRBASE);
    %     figure; [AUC(i),ShuffleRank,NullDist]=Analyze.TwoGroup.ShuffleCompare(FR(idx),FRBASE,'Type','AUC','plotFig');
    [AUC(i),ShuffleRank,NullDist]=Analyze.TwoGroup.ShuffleCompare(FR(idx),FRBASE,'Type','AUC','NShuffles',0);
    [dPrime(i),ShuffleRank,NullDist]=Analyze.TwoGroup.ShuffleCompare(FR(idx),FRBASE,'Type','dprime','NShuffles',0);
  
Comp(i) = ranksum(FR(idx),FRBaseLine);
    %     figure; [AUC(i),ShuffleRank,NullDist]=Analyze.TwoGroup.ShuffleCompare(FR(idx),FRBASE,'Type','AUC','plotFig');
    [ComAUC(i),ShuffleRank,NullDist]=Analyze.TwoGroup.ShuffleCompare(FR(idx),FRBaseLine,'Type','AUC','NShuffles',0);
    [ComdPrime(i),ShuffleRank,NullDist]=Analyze.TwoGroup.ShuffleCompare(FR(idx),FRBaseLine,'Type','dprime','NShuffles',0);
    
    
    CompAll(i) = ranksum(FR(idx),FR(idxAll));
    %     figure; [AUC(i),ShuffleRank,NullDist]=Analyze.TwoGroup.ShuffleCompare(FR(idx),FRBASE,'Type','AUC','plotFig');
    [AUCvAll(i),ShuffleRank,NullDist]=Analyze.TwoGroup.ShuffleCompare(FR(idx),FR(idxAll),'Type','AUC','NShuffles',0);

MeanRate(i)=mean(FR(idx));
MeanFromBase(i)=mean(FR(idx))-mean(FRBaseLine);
end


Results.CompAll=CompAll;
Results.AUCvAll=AUCvAll;

Results.MeanRate=MeanRate;
Results.MeanFromBase=MeanFromBase;
Results.Comp=Comp;
Results.ComAUC=ComAUC;
Results.ComdPrime=ComdPrime;
    Results.p=p;
Results.AUC=AUC;
Results.dPrime=dPrime;

end
%%
%% 

