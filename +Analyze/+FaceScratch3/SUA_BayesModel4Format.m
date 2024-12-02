function [Results]=SUA_BayesModel4Format(AllTrialData,Phase,unit,timeWindow,opts)
%%(AllTrialData,Phase,unit,timeWindow,opts)
% unit=[1 50 1];
% % Phase='Go1';
% timeWindow=[.5 1];
% opts.BaselinePhase='CueTarget'
model= 'FR~Condition';

LabelOrder={'Pinch'    'Press'    'Rub'    'Tap'};

cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase);
% Shuffles=Analyze.ActionObsInvarianceText.CrossDecodeHelper.ComputePossibleUniquePermutations4();

%%


% Conditions is basically FeltCheek, FeltShoulder,ObsCheek,ObsShoulder

Conditions={ {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek'},...
    {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'},...
    {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'},...
    {'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'}};


%% Within Class Stats
for i=1:length(Conditions)
    cCond=Conditions{i};
    cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
    
    FRBaseLine=Analyze.getNeuralData(Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition','Null'), unit,[.5 2.5]);
    
    L=Analyze.returnFieldValues(cTrialData,'Action');
    %     if opts.shuffle
    %         %             L=L(randperm(length(L)));
    %         L=Analyze.Semantic.CrossDecodeHelper.ShuffleConditionLabels(L,Shuffles(i,:));
    %     end
    FR=Analyze.getNeuralData(cTrialData, unit,timeWindow);
    
    [Condition,Labels]=DataConvert.labels2IDX(L,LabelOrder);
    
    
    FR=[FR];
    Condition=[Condition(:)];
    %         if opts.shuffle
    %             Condition=Condition(randperm(length(Condition)));
    %         end
    dataTable = table(Condition,FR);
    dataTable.Condition=nominal(dataTable.Condition);
    %     [~,mdlText]=Analyze.FitGLM(data,model);
    dT{i}=dataTable;
    
end

FR=[FRBaseLine];
Condition=[FRBaseLine*0];
%     if opts.shuffle
%         Condition=Condition(randperm(length(Condition)));
%     end
dataTable = table(Condition,FR);
dataTable.Condition=nominal(dataTable.Condition);

NoiseData=dataTable;
%%

% mdls are the following:
% 1. invariant
% 2. idiosyncratic
% 3. BP-Nancy Only
% 4. BP-Tyson Only
% 5. BP
% 6. Person-Cheek Only
% 7. Person-Shoulder Only
% 8. Person

% mdls = {...
%     [1 1 1 1];...
%     [1 2 3 4];...
%     [1 2 0 0];...
%     [0 0 1 2];...
%     [1 2 1 2];...
%     [1 0 2 0];...
%     [0 1 0 2];...
%     [1 1 2 2];...
%     };
%
mdls={...
    [1 1 1], [1,2,0],1;...
    
    [1 2 1], [1 2 3 0],3;...
    [1 1 2], [1 2 3 0],2;...
    [2 1 1], [1 2 3 0],4;...
    
    [1 0 1],  [1 2 0],3;...
    [1 1 0],  [1 2 0],2;... 1 2 0
    [0 1 1],  [1 2 0],4;... 1 2 0
    
    [1 0 0],[1 2],5;...;... 1 2 0
    [0 1 0],[1 2],6;...;... 1 2 0
    [0 0 1],[1 2],7;...;... 1 2 0
    
    [1 0 0],0,5;...
    [0 1 0],0,6;...
    [0 0 1],0,7;...
    [0 0 0],1,8;...
    
    
    [1 2 3], [1,2,3,4,0],8;...
    %
    [1 2 0],[1 2 3 0],8;...
    [1 0 2],[1 2 3 0],8;... 1 2 3 0
    [0 1 2],[1 2 3 0],8;... 1 2 3 0
    };
%
AllMdls=[];
for i=1:size(mdls,1)
    for j=1:length(mdls{i,2})
        AllMdls=[AllMdls;mdls{i,1} mdls{i,2}(j) mdls{i,3}(1)];
    end
end
% whos AllMdls
AllMdls=[AllMdls(:,1:3) AllMdls(:,end-1)];
mdls=AllMdls;

% AllMdls = [];
% for i=1:size(mdls,1)
%     AllMdls = [AllMdls;mdls{i}];
% end
% mdls = AllMdls;

%%
% cmpEqual=[ (mdls(:,2)==mdls(:,3) & mdls(:,2)~=0) (mdls(:,2)==mdls(:,4) & mdls(:,2)~=0) (mdls(:,3)==mdls(:,4)  & mdls(:,3)~=0)];
%
% All=find(sum(cmpEqual,2)==3 );
% FL0=find(cmpEqual(:,1)==1 & ~(sum(cmpEqual,2)==3) );
% FL1=find(cmpEqual(:,2)==1 & ~(sum(cmpEqual,2)==3) );
% L0L1=find(cmpEqual(:,3)==1 & ~(sum(cmpEqual,2)==3) );
%
% Idio=find(sum(cmpEqual,2)==0 & sum(mdls(:,2:4),2)>1   );
% F=find(sum(cmpEqual,2)==0 & sum(mdls(:,2:4),2)==1 & mdls(:,2)==1   );
% L0=find(sum(cmpEqual,2)==0 & sum(mdls(:,2:4),2)==1 & mdls(:,3)==1   );
% L1=find(sum(cmpEqual,2)==0 & sum(mdls(:,2:4),2)==1 & mdls(:,4)==1   );
% TOnly=find(sum(mdls(:,2:4),2)==0);
%
% length(F) +length(L0)+length(L1)+length(Idio);
% length(All) +length(FL0)+length(FL1)+length(L0L1);
%

%%
for mdlIDX=1:size(mdls,1)
    %     sprintf('%d, ' ,mdlIDX)
    
    %     mdlIDX=7
    cdT=dT;
    cm=mdls(mdlIDX,:);
    %     pOffset=[];
    lastVal=0;
    for idx=1:4
        
        if cm(idx)==0
            offset(idx)= lastVal+1;
            cdT{idx}.Condition= nominal(double(cdT{idx}.Condition)*0+offset(idx));
            lastVal=lastVal+1;
        else
            
            %         offset=(5*(cm(idx)-1));
            if idx>1 && any(cm(idx)==cm(1:idx-1))
                tmp=find(cm(idx)==cm(1:idx-1));
                offset(idx)=offset(tmp(1));
            else
                offset(idx)=[lastVal];
            end
            cdT{idx}.Condition= nominal(double(cdT{idx}.Condition)+offset(idx));
            lastVal=lastVal+4;
        end
        %         disp([num2str(idx) ':' num2str(cm(idx))])
        %         disp(unique(cdT{idx}.Condition)')
    end
    
    %     T1=[cdT{1};cdT{2};cdT{3};cdT{4}];
    
    T1=[NoiseData;cdT{1};cdT{2};cdT{3};cdT{4}];
    %%
    
%     M1{mdlIDX}=fitlme(T1,'FR~Condition','Weights',[repmat(1,1,20),repmat(2,1,40),repmat(3,1,40),repmat(4,1,40),repmat(4,1,40)]);
    
    M1{mdlIDX}=fitlme(T1,'FR~Condition');
    [mdlFit1{mdlIDX},mdlText{mdlIDX}]=Analyze.FitGLM(T1,'FR~Condition','CrossValidate','BaseLineInds',Condition==0);
    %          [mdlFit1{mdlIDX},mdlText{mdlIDX}]=Analyze.FitGLM(T1,'FR~Condition');
%         [mdlFit1{mdlIDX},mdlText{mdlIDX}]=Analyze.FitGLM(T1,'FR~Condition','CrossValidate');

    
    Results.(['AIC'])(mdlIDX)=M1{mdlIDX}.ModelCriterion.AIC;
    Results.(['BIC' ])(mdlIDX)=M1{mdlIDX}.ModelCriterion.BIC;
    
    Results.(['adjR2' ])(mdlIDX)=mdlFit1{mdlIDX}.adjR2;
    Results.(['R2' ])(mdlIDX)=mdlFit1{mdlIDX}.R2;
    Results.(['Pval' ])(mdlIDX)=mdlFit1{mdlIDX}.Pval;
    Results.(['ErrRedCoef' ])(mdlIDX)=mdlFit1{mdlIDX}.ErrRedCoef;
    Results.(['CVcorr' ])(mdlIDX)=mdlFit1{mdlIDX}.CVcorr;
    
    Results.(['CVR2' ])(mdlIDX)=mdlFit1{mdlIDX}.CVR2;
    
    
end

%%
% now demean
for i=1:4
    dT{i}.FR=   dT{i}.FR- mean(dT{i}.FR);
end
NoiseData.FR=NoiseData.FR-mean(NoiseData.FR);
%%
for mdlIDX=1:size(mdls,1)
    %     sprintf('%d, ' ,mdlIDX)
    
    %     mdlIDX=7
    cdT=dT;
    cm=mdls(mdlIDX,:);
    pOffset=[];
    lastVal=0;
    for idx=1:4
        
        if cm(idx)==0
            offset(idx)= lastVal+1;
            cdT{idx}.Condition= nominal(double(cdT{idx}.Condition)*0+offset(idx));
            lastVal=lastVal+1;
        else
            
            %         offset=(5*(cm(idx)-1));
            if idx>1 && any(cm(idx)==cm(1:idx-1))
                tmp=find(cm(idx)==cm(1:idx-1));
                offset(idx)=offset(tmp(1));
            else
                offset(idx)=[lastVal];
            end
            cdT{idx}.Condition= nominal(double(cdT{idx}.Condition)+offset(idx));
            lastVal=lastVal+5;
        end
        %         disp([num2str(idx) ':' num2str(cm(idx))])
        %         disp(unique(cdT{idx}.Condition)')
    end
    
    T1=[NoiseData;cdT{1};cdT{2};cdT{3};cdT{4}];
    %%
    
    
    
    
    M1{mdlIDX}=fitlme(T1,'FR~Condition');
    Results.(['nmAIC'])(mdlIDX)=M1{mdlIDX}.ModelCriterion.AIC;
    Results.(['nmBIC' ])(mdlIDX)=M1{mdlIDX}.ModelCriterion.BIC;
    
    
    [mdlFit1{mdlIDX},mdlText{mdlIDX}]=Analyze.FitGLM(T1,'FR~Condition','CrossValidate', 'BaseLineInds',Condition==0);
    Results.(['nmadjR2' ])(mdlIDX)=mdlFit1{mdlIDX}.adjR2;
    Results.(['nmR2' ])(mdlIDX)=mdlFit1{mdlIDX}.R2;
    Results.(['nmPval' ])(mdlIDX)=mdlFit1{mdlIDX}.Pval;
    Results.(['nmErrRedCoef' ])(mdlIDX)=mdlFit1{mdlIDX}.ErrRedCoef;
    Results.(['nmCVcorr' ])(mdlIDX)=mdlFit1{mdlIDX}.CVcorr;
    Results.(['nmCVR2' ])(mdlIDX)=mdlFit1{mdlIDX}.CVR2;
    
    
end

end

