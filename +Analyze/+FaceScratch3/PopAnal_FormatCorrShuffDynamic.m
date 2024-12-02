function Results=PopAnal_FormatCorrShuffDynamic(AllTrialData,Phase,unit,timeWindow,opts)
%%
% timeWindow=[-1.5 -.5];
% Task related Firing Rate
disp(timeWindow)


%%
Conditions={ 'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek',...
    'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder',...
    'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek',...
    'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'};


useSplits=0;
nReps=opts.nReps;
winSize=opts.winSize;
timeWindow=opts.timeWindow;

kern=Kinematics.MinJerkKernel(winSize+.1, .1);
MinJerk=1;
n=[repmat(10,[1 16])];
for repIDX=1:nReps
    for nC=1:16
        rP{nC,repIDX}=randperm(n(nC));
    end
end

%% Within Class Stats
for timeIDX=1:length(timeWindow)
    clear FR
    cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase);
    %     Condition=Analyze.returnUniqueFieldValues(cTrialData,'Condition');
    %     Condition=Condition(cellfun(@isempty,strfind(Condition,'Null')));
    %     Condition=[Condition(end-4:end); Condition(1:15)];
    for cIDX=1:length(Conditions)
        cCond=Conditions{cIDX};
        
        cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
        
        if MinJerk
            tmp=squeeze(Analyze.getNeuralData(cTrialData, unit,[(timeWindow(timeIDX)-winSize/2-0.05):.1:(timeWindow(timeIDX)+winSize/2+.05)]));
            k2=repmat(kern',[size(tmp,1),1,size(tmp,3)]);
            tmp=tmp.*k2; tmp=squeeze(sum(tmp,2));
            FR{cIDX}=tmp;
        else
            
            FR{cIDX}=squeeze(Analyze.getNeuralData(cTrialData, unit,[timeWindow(timeIDX)-winSize/2 timeWindow(timeIDX)+winSize/2]));
        end
    end
    
    
    
    
    
    gr={1:4,5:8,9:12,13:16};
    
    % FRdm=cellfun(@(x)x(:,H),FRdm,'UniformOutput',false);
    
    if useSplits==0
        mu=cellfun(@mean,FR,'UniformOutput',false);
        mu=cat(1,mu{:});
        clear FRdm
        for i=1:length(gr)
            FRdm(gr{i})=cellfun(@(x)(x-mean(mu(gr{i},:),1)),FR(gr{i}),'UniformOutput',false);
            
        end
        %         mu=cellfun(@mean,FRdm,'UniformOutput',false);
        mu=cellfun(@mean,FR,'UniformOutput',false);
        mu=cat(1,mu{:});
        
        clear mu2
        for i=1:length(gr)
            tmp=mu(gr{i},:);
            mu2(:,i)=tmp(:);
            
        end
        
        
        MeanResponseA{timeIDX}{1}=mu2;
        MeanResponseB{timeIDX}{1}=mu2;
        
    else
        
        
        for repIDX=1:nReps
            
            
            clear muA muB FRA FRB
            for nC=1:16
                
                FRA{nC}=FR{nC}(rP{nC,repIDX}(1:n(nC)/2),:);
                
                FRB{nC}=FR{nC}(rP{nC,repIDX}(n(nC)/2+1:end),:);
                
            end
            
            
            clear muA muB
            for nC=1:16
                muA(nC,:)=mean(FRA{nC},1);
                muB(nC,:)=mean(FRB{nC},1);
            end
            
            
            clear FRdmA FRdmB
            for i=1:length(gr)
                FRdmA(gr{i})=cellfun(@(x)(x-mean(muA(gr{i},:),1)),FRA(gr{i}),'UniformOutput',false);
                FRdmB(gr{i})=cellfun(@(x)(x-mean(muB(gr{i},:),1)),FRB(gr{i}),'UniformOutput',false);
                
                %                 FRdmA(gr{i})=cellfun(@(x)(x),FRA(gr{i}),'UniformOutput',false);
                %                 FRdmB(gr{i})=cellfun(@(x)(x),FRB(gr{i}),'UniformOutput',false);
                
            end
            
            clear muA muB
            for nC=1:16
                muA(nC,:)=mean(FRdmA{nC},1);
                muB(nC,:)=mean(FRdmB{nC},1);
                
                %                 muA(nC,:)=mean(FRA{nC},1);
                %                 muB(nC,:)=mean(FRB{nC},1);
            end
            
            clear muCombA muCombB
            for i=1:length(gr)
                tmp=muA(gr{i},:);
                muCombA(:,i)=tmp(:);
                
                tmp=muB(gr{i},:);
                muCombB(:,i)=tmp(:);
            end
            
            
            MeanResponseA{timeIDX}{repIDX}=muCombA;
            MeanResponseB{timeIDX}{repIDX}=muCombB;
        end
    end
end
%%
for timeIDX1=1:length(timeWindow)
    for timeIDX2=1:length(timeWindow)
        
        for format1=1:4
            for format2=1:4
                
                if length(MeanResponseA{timeIDX})==1
                    d1=MeanResponseA{timeIDX1}{1}(:,format1);
                    d2=MeanResponseB{timeIDX2}{1}(:,format2);
                    [tmpcorr(1),tmpP(1)]=corr(d1,d2);
                    
                    
                    tmpA = MeanResponseA{timeIDX1}{1};
                    tmpB = MeanResponseB{timeIDX2}{1};
                    
                    if isequal(tmpA,tmpB)
                        [tmp1 tmpPpartcorr] = partialcorr(tmpA);
                        Ptmpcorr(1) = tmp1(format1,format2);
                        PtmpP(1) = tmpPpartcorr(format1,format2);
                    else
                        [tmp1 tmpPpartcorr] = partialcorr([tmpA tmpB]);
                        Ptmpcorr(1) = tmp1(format1,format2+4);
                        PtmpP(1) = tmpPpartcorr(format1,format2+4);
                    end
                else
                    
                    for rep=1:length(MeanResponseA{timeIDX})
                        d1=MeanResponseA{timeIDX1}{rep}(:,format1);
                        d2=MeanResponseB{timeIDX2}{rep}(:,format2);
                        [tmpcorr(rep),tmpP(rep)]=corr(d1,d2);
                        
                        %                     also compute the partial correlation
                        [tmp1 tmpPpartcorr] = partialcorr([MeanResponseA{timeIDX1}{rep}(:,format1) MeanResponseB{timeIDX2}{rep}(:,format2)],...
                            [MeanResponseA{timeIDX1}{rep}(:,setdiff([1 2 3 4],format1)) MeanResponseB{timeIDX2}{rep}(:,setdiff([1 2 3 4],format2))]);
                        Ptmpcorr(rep)=tmp1(1,2);
                        PtmpP(rep) = tmpPpartcorr(1,2);
                        
                    end
                end
                F1{format1,format2}(timeIDX1,timeIDX2)=mean(tmpcorr);
                P1{format1,format2}(timeIDX1,timeIDX2)=mean(tmpP);
                
                F1pc{format1,format2}(timeIDX1,timeIDX2)=mean(Ptmpcorr);
                P1pc{format1,format2}(timeIDX1,timeIDX2)=mean(PtmpP);
            end
        end
        
        
    end
end

Results.FormatCorrs=F1;
Results.FormatCorrP=P1;

Results.FormatCorrsPC=F1pc;
Results.FormatCorrPPC=P1pc;

