function Results=PopAnal_FormatDecodeDynamic(AllTrialData,Phase,unit,timeWindow,opts)
%%
% timeWindow=[-1.5 -.5];
% Task related Firing Rate
disp(timeWindow)


%%
Conditions = {{'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek'},...
    {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'},...
    {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'},...
    {'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'}};

useSplits=1;
nReps=opts.nReps;
winSize=opts.winSize;
timeWindow=opts.timeWindow;

kern=Kinematics.MinJerkKernel(winSize+.1, .1);
MinJerk=1;
% n=[repmat(10,[1 16])];
% for repIDX=1:nReps
%     for nC=1:4
%         rP{nC,repIDX}=randperm(n(nC));
%     end
% end

%% Within Class Stats
for timeIDX=1:length(timeWindow)
    clear FR
    cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase);
    
    
    for condIDX=1:length(Conditions)
        cCond=Conditions{condIDX};
        cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
        
        L{condIDX}=Analyze.returnFieldValues(cTrialData,'Action');
        
        
        
        if MinJerk
            tmp=squeeze(Analyze.getNeuralData(cTrialData, unit,[(timeWindow(timeIDX)-winSize/2-0.05):.1:(timeWindow(timeIDX)+winSize/2+.05)]));
            k2=repmat(kern',[size(tmp,1),1,size(tmp,3)]);
            tmp=tmp.*k2; tmp=squeeze(sum(tmp,2));
            FR{condIDX}=tmp;
        else
            
            FR{condIDX}=squeeze(Analyze.getNeuralData(cTrialData, unit,[timeWindow(timeIDX)-winSize/2 timeWindow(timeIDX)+winSize/2]));
        end
        
       
        
        [FR{condIDX},L{condIDX}]=Analyze.SampleFeaturePopulation(FR{condIDX},L{condIDX},'Type','Basic');
        
        
        [FRCell{condIDX,timeIDX},LabelCell{condIDX,timeIDX}]=DataConvert.list2cell(FR{condIDX},L{condIDX});
        
        
         FRCell{condIDX,timeIDX}= cellfun(@(x)x(1:10,:),FRCell{condIDX,timeIDX},'UniformOutput',false);
        
        [FRTrue{condIDX,timeIDX},LabelTrue{condIDX,timeIDX}]=DataConvert.cell2list(FRCell{condIDX,timeIDX});
    end
end

%%
for timeIDX1=1:length(timeWindow)
    for timeIDX2=1:length(timeWindow)

        disp(timeIDX2)
        for format1=1:4
            for format2=1:4
                d1=FRTrue{format1,timeIDX1};
                
                d2=FRTrue{format2,timeIDX2};
                l1=LabelTrue{format1,timeIDX1};
                l2=LabelTrue{format2,timeIDX1};
                clear l2_pred l2_test
                for rep=1:size(d1,1)
                    trainIDX=setdiff(1:size(d1,1),rep);
                    
                    d1_test=d1(rep,:);
                    d1_train=d1(trainIDX,:);
                    
                    l1_test=l1(rep,:);
                    l1_train=l1(trainIDX,:);
                    if rep>size(d2,1)
                        break
                    end
                    d2_test=d2(rep,:);
                    l2_test(rep)=l2(rep,:);
                    
                    [FR,L]=DataConvert.list2cell(d1_train,l1_train);
                    
                    opts.dec.Train('TrainingData',{FR,[1:4],{'1','2','3','4'}});
                    l2_pred(rep)=opts.dec.PredictBatch(d2_test);
                end
                %                 CompleteStats=confusionmatStats(l2_test,l2_pred);
                cvAccuracy{format1,format2}(timeIDX1,timeIDX2)=(nnz(l2_test==l2_pred)/length(l2_pred))*100;
                %                 NShuffs=1000;
                %                   for shuffIDX=1:NShuffs
                %                 ChancePercent(shuffIDX)=(nnz(l2_test==Shuffle(l2_pred))/length(l2_pred))*100;
                %             end
            end
        end
        
        
    end
end

Results.cvAccuracy=cvAccuracy;
