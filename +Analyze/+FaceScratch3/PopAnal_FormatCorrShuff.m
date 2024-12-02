function Results=PopAnal_FormatCorrShuff(AllTrialData,Phase,unit,timeWindow,opts)
%%
% timeWindow=[-1.5 -.5];
% Task related Firing Rate
disp(timeWindow)


%%
Conditions={ 'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek',...
    'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder',...
    'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek',...
    'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'};

for i=1:16
    for unitIDX=1:size(unit,1)
        cCond=Conditions{i};
        cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
        
        
        L=Analyze.returnFieldValues(cTrialData,'Action');
        
        FR=Analyze.getNeuralData(cTrialData, unit(unitIDX,:),timeWindow);
       
    end
    
end

%% Within Class Stats

clear FR
cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase);
% Condition=Analyze.returnUniqueFieldValues(cTrialData,'PersonBodyPart');
% Condition=Condition(cellfun(@isempty,strfind(Condition,'Null')));
% Condition=[Condition(end-4:end); Condition(1:15)];
for cIDX=1:length(Conditions)
    cCond=Conditions{cIDX};
    
    cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
    
    FR{cIDX}=squeeze(Analyze.getNeuralData(cTrialData, unit,timeWindow));
    
end

FR=cellfun(@(x)x(1:10,:),FR,'UniformOutput',false);

%%
mu=cellfun(@mean,FR,'UniformOutput',false);
sigma=cellfun(@std,FR,'UniformOutput',false);
mu=cat(1,mu{:});
sigma=mean(cat(1,sigma{:}),1);


%%
mu=cellfun(@mean,FR,'UniformOutput',false);
mu=cat(1,mu{:});

gr={1:4,5:8,9:12,13:16};
clear FRdm
for i=1:length(gr)
    FRdm(gr{i})=cellfun(@(x)(x-mean(mu(gr{i},:),1)),FR(gr{i}),'UniformOutput',false);
    
end
% FRdm=cellfun(@(x)x(:,H),FRdm,'UniformOutput',false);

mu=cellfun(@mean,FR,'UniformOutput',false);
mu=cat(1,mu{:});

gr={1:4,5:8,9:12,13:16};
clear mu2
for i=1:length(gr)
    tmp=mu(gr{i},:);
    mu2(:,i)=tmp(:);
    
end

Results.corr=corr(mu2);
Results.corrPart=partialcorr(mu2);


Results.corrCI=bootci(2000,{@corr,mu2},'type','per');
%%
A=[2 3 4];B=[1 3 4]; C=[1 2 3]; D=[1 2 4];
[A,B,C,D]=ndgrid(A,B,C,D);
tmp=[A(:),B(:),C(:),D(:)];
for i=1:size(tmp,1)
    N(i)=length(unique(tmp(i,:)));
end

Shuffles=tmp(find(N==4),:);

% Shuffles=Analyze.Semantic.CrossDecodeHelper.ComputePossibleUniquePermutations();

for i=1:4
    for j=1:4
        tmp1=mu(gr{i},:);
        idx2=gr{j};
        for kk=1:size(Shuffles,1)
            shuffIDX=idx2(Shuffles(kk,:));
                    tmp2=mu(shuffIDX,:);
            
            CC(i,j,kk)=corr(tmp1(:),tmp2(:));
        end
        
    end
    
end

Results.corrShuff=mean(CC,3);