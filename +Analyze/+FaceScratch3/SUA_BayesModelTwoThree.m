function [Results]=SUA_BayesModelTwoThree(AllTrialData,Phase,unit,timeWindow,opts)
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

Fc = {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek'};
Fs = {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'};
Oc = {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'};
Os = {'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'};

Conditions={Fc, Fs, Oc, Os};
Names = {'Fc','Fs','Oc';...
    'Fc', 'Oc', 'Os';...
    'Fc', 'Fs', 'Os';...
    'Fs', 'Oc', 'Os';...
    };
    
AllMDL=[];

cmp={Fc Fs Oc;...
    Fc Oc Os;...
    Fc Fs Os;...
    Fs Oc Os;...
    };


for i=1:size(cmp,1)
    Conditions=cmp(i,:);
    cNames = Names(i,:);
    
    %% Within Class Stats
    for i=1:length(Conditions)
        cCond=Conditions{i};
        cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
        
        FRBaseLine=Analyze.getNeuralData(Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition','Null'), unit,[.5 2.5]);
        
        L=Analyze.returnFieldValues(cTrialData,'Action');
%         if opts.shuffle
%             %             L=L(randperm(length(L)));
%             L=Analyze.Semantic.CrossDecodeHelper.ShuffleConditionLabels(L,Shuffles(i,:));
%         end
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
    mdls=[...
        [1 1 1];...
        [1 2 3];...
        
        [1 2 1];...
        [1 1 2];...
        [2 1 1];...
        
        [1 0 1];...
        [1 1 0];...
        [0 1 1];...
        
        [1 2 0];...
        [1 0 2];...
        [0 1 2];...
        
        [1 0 0];...
        [0 1 0];...
        [0 0 1];...
        ];
    
    for mdlIDX=1:size(mdls,1)
        %%
        cdT=dT;
        
        if isequal(mdls(mdlIDX,:),[1 1 1])
            T1=[NoiseData;cdT{1} ;cdT{2};cdT{3}];
            
        elseif isequal(mdls(mdlIDX,:),[1 2 3])
            cdT{2}.Condition= nominal(double(cdT{2}.Condition)+4);
            cdT{3}.Condition= nominal(double(cdT{3}.Condition)+8);
            T1=[NoiseData;cdT{1} ;cdT{2};cdT{3}];
            
        elseif sum(mdls(mdlIDX,:)==[1 1 1])==2 & sum(mdls(mdlIDX,:)==[2 2 2])==1
            %                 [1 2 1];... [1 1 2];... [1 2 2];...
            oddIDX=find(mdls(mdlIDX,:)==[2 2 2]);
            cdT{oddIDX}.Condition= nominal(double(cdT{oddIDX}.Condition)+4);
            T1=[NoiseData;cdT{1} ;cdT{2};cdT{3}];
            
        elseif sum(mdls(mdlIDX,:)==[1 1 1])==2 & sum(mdls(mdlIDX,:)==[0 0 0])==1
            %                 [1 0 1];... [1 1 0];... [0 1 1];...
            oddIDX=find(mdls(mdlIDX,:)==[0 0 0]);
            cdT{oddIDX}.Condition= nominal(double(cdT{oddIDX}.Condition)*0+11);
            T1=[NoiseData;cdT{1} ;cdT{2};cdT{3}];
            
            
        elseif sum(mdls(mdlIDX,:)==[1 1 1])==1 & sum(mdls(mdlIDX,:)==[0 0 0])==1 & sum(mdls(mdlIDX,:)==[2 2 2])==1
            %                 [1 0 1];...   [1 1 0];...  [0 1 1];...
            oddIDX1=find(mdls(mdlIDX,:)==[2 2 2]);
            cdT{oddIDX1}.Condition= nominal(double(cdT{oddIDX1}.Condition)+5);
            oddIDX2=find(mdls(mdlIDX,:)==[0 0 0]);
            
            cdT{oddIDX2}.Condition= nominal(double(cdT{oddIDX2}.Condition)*0+6);
            T1=[NoiseData;cdT{1} ;cdT{2};cdT{3}];
            
        elseif sum(mdls(mdlIDX,:)==[1 1 1])==1 & sum(mdls(mdlIDX,:)==[0 0 0])==2
            %                 [1 0 1];... [1 1 0];...   [0 1 1];...
            oddIDX=find(mdls(mdlIDX,:)==[0 0 0]);
            cdT{oddIDX(1)}.Condition= nominal(double(cdT{oddIDX(1)}.Condition)*0+6);
            cdT{oddIDX(2)}.Condition= nominal(double(cdT{oddIDX(2)}.Condition)*0+7);
            T1=[NoiseData;cdT{1} ;cdT{2};cdT{3}];
            
        end
        
        M1{mdlIDX}=fitlme(T1,'FR~Condition');
        Condition=double(T1.Condition)-1;
        [mdlFit1{mdlIDX},mdlText{mdlIDX}]=Analyze.FitGLM(T1,'FR~Condition','CrossValidate','BaseLineInds',Condition==0);
        
        cndName=[cNames{1} cNames{2} cNames{3}];
        Results.(['AIC_' cndName])(mdlIDX)=M1{mdlIDX}.ModelCriterion.AIC;
        Results.(['BIC_' cndName])(mdlIDX)=M1{mdlIDX}.ModelCriterion.BIC;
        Results.(['CVR2_' cndName])(mdlIDX)=mdlFit1{mdlIDX}.CVR2;
        
        Results.(['CVR2_' cndName])(mdlIDX)=mdlFit1{mdlIDX}.CVR2;
        
        Results.(['ErrRedCoef_' cndName])(mdlIDX)=mdlFit1{mdlIDX}.ErrRedCoef;
        
        
    end
    
    %%
    
    
%     [v,idx]=min(Results.(['BIC_' cndName]));
%     clear ptmp
%     for i=1:size(mdls,1)
%         p=compare(M1{i},M1{idx});
%         ptmp(i)=p.pValue(2);
%     end
%     
%     Results.(['pBIC_' cndName])=ptmp;
end
%%


Fc = {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek'};
Fs = {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'};
Oc = {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'};
Os = {'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'};

Conditions={Fc, Fs, Oc, Os};
Names = {'Fc', 'Fs';...
    'Fc', 'Oc';...
    'Fc', 'Os';...
    'Fs', 'Oc';...
    'Fs', 'Os';...
    'Oc', 'Os';...
    };
    
AllMDL=[];

cmp={Fc Fs;...
    Fc Oc;...
    Fc Os;...
    Fs Oc;...
    Fs Os;...
    Oc Os;...
    };

% cmp={'L0','L1';...
%     'L0','F';...
%     'L0','T';...
%     'L1','F';...
%     'L1','T';...
%     'F','T'};
%%

for i=1:size(cmp,1)
    Conditions=cmp(i,:);
    cNames = Names(i,:);
        
    %% Within Class Stats
    for i=1:length(Conditions)
        cCond=Conditions{i};
        cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
        
        FRBaseLine=Analyze.getNeuralData(Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition','Null'), unit,[.5 2.5]);

        L=Analyze.returnFieldValues(cTrialData,'Action');
%         if opts.shuffle
%             %             L=L(randperm(length(L)));
%             L=Analyze.Semantic.CrossDecodeHelper.ShuffleConditionLabels(L,Shuffles(i,:));
%         end
        FR=Analyze.getNeuralData(cTrialData, unit,timeWindow);
        
        [Condition,Labels]=DataConvert.labels2IDX(L,LabelOrder);
        
        
        FR=[FR;FRBaseLine];
        Condition=[Condition(:);FRBaseLine*0];
        %         if opts.shuffle
        %             Condition=Condition(randperm(length(Condition)));
        %         end
        dataTable = table(Condition,FR);
        dataTable.Condition=nominal(dataTable.Condition);
        %     [~,mdlText]=Analyze.FitGLM(data,model);
        dT{i}=dataTable;
        
    end
    %%
    
    %% compute models
    % Assume same tuning to both formats
    T1=[dT{1} ;dT{2}];
    M1=fitlme(T1,'FR~Condition');
    Condition=double(T1.Condition)-1;
    [mdlFit1,mdlText]=Analyze.FitGLM(T1,'FR~Condition','CrossValidate','BaseLineInds',Condition==0);
    
    
    % Assume tuning only to condition 2
    tmp=dT{1};idx=double(tmp.Condition)>1;
    tmp.Condition(idx)=nominal(6);
    T2=[tmp ;dT{2}];
    M2=fitlme(T2,'FR~Condition');
    Condition=double(T2.Condition)-1;
    [mdlFit2,mdlText]=Analyze.FitGLM(T2,'FR~Condition','CrossValidate','BaseLineInds',Condition==0);
    
    
    % Assume tuning only to condition 1
    tmp=dT{2};idx=double(tmp.Condition)>1;
    tmp.Condition(idx)=nominal(6);
    T3=[dT{1}; tmp];
    M3=fitlme(T3,'FR~Condition');
    Condition=double(T3.Condition)-1;
    [mdlFit3,mdlText]=Analyze.FitGLM(T3,'FR~Condition','CrossValidate','BaseLineInds',Condition==0);
    
    
    % Assume idiosyncratic tuning to both formats
    tmp=dT{2};idx=double(tmp.Condition)>1;
    tmp.Condition(idx)=nominal(double(tmp.Condition(idx)) +4);
    T4=[dT{1}; tmp];
    
    M4=fitlme(T4,'FR~Condition');
    Condition=double(T4.Condition)-1;
    [mdlFit4,mdlText]=Analyze.FitGLM(T4,'FR~Condition','CrossValidate','BaseLineInds',Condition==0);
    
    cndName=[cNames{1} cNames{2}];
    Results.(['AIC_' cndName])=[M1.ModelCriterion.AIC M2.ModelCriterion.AIC M3.ModelCriterion.AIC M4.ModelCriterion.AIC];
    Results.(['BIC_' cndName])=[M1.ModelCriterion.BIC M2.ModelCriterion.BIC M3.ModelCriterion.BIC M4.ModelCriterion.BIC];
    Results.(['LL_' cndName])=[M1.ModelCriterion.LogLikelihood M2.ModelCriterion.LogLikelihood M3.ModelCriterion.LogLikelihood M4.ModelCriterion.LogLikelihood];
    Results.(['Deviance_' cndName])=[M1.ModelCriterion.Deviance M2.ModelCriterion.Deviance M3.ModelCriterion.Deviance M4.ModelCriterion.Deviance];
    
    
    Results.(['CVR2_' cndName])=[mdlFit1.CVR2 mdlFit2.CVR2 mdlFit3.CVR2 mdlFit4.CVR2];
    Results.(['CVR2p_' cndName])=[mdlFit1.CVR2p mdlFit2.CVR2p mdlFit3.CVR2p mdlFit4.CVR2p];
    Results.(['ErrRedCoef_' cndName])=[mdlFit1.ErrRedCoef mdlFit2.ErrRedCoef mdlFit3.ErrRedCoef mdlFit4.ErrRedCoef];
    
    %%
%     MM={M1 M2 M3 M4};
%     
%     [v,idx]=min(Results.(['BIC_' cndName]));
%     clear ptmp
%     for i=1:4
%         if MM{idx}.NumCoefficients>=MM{i}.NumCoefficients
%             [p,siminfo]=compare(MM{idx},MM{i});
%         else
%             [p,siminfo]=compare(MM{i},MM{idx});
%         end
%         ptmp(i)=p.pValue(2);
%     end
%     %%
%     Results.(['pBIC_' cndName])=ptmp;
    
end

end
