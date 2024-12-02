function [FRorig, FR] = getDataForDPCAAction(AllTrialData,timeWindow,Phase,varargin)

[varargin,SmoothingKernel]=Utilities.ProcVarargin(varargin,'SmoothingKernel',.55);
[varargin,meanThresh]=Utilities.ProcVarargin(varargin,'meanThresh',1.5);
[varargin,TaskLabel]=Utilities.ProcVarargin(varargin,'TaskLabel','XActionType');

datafile = fullfile(env.get('results'),'FaceScratch3','SUAnal',[TaskLabel '-Go'],'PopData','p1-AllDays.mat');
load(datafile); % Loads PlotData
ASF = vertcat(PlotData{1}.Unit{:});
UDates = Analyze.returnUniqueFieldValues(AllTrialData,'Date');

smoothArg=[SmoothingKernel mean(diff(timeWindow))];
MeanRate=Analyze.getNeuralData(AllTrialData, ASF,'MeanRate');

goodFiring=nanmean(MeanRate,2)>meanThresh;

bodyNames = {'Cheek','Shoulder'};
personNames = {'Nancy','Tyson'};
actionNames = {'Pinch','Press','Rub','Tap'};

FRorig = cell(length(bodyNames),1);
FR = cell(length(bodyNames),1);

for bIdx=1:length(bodyNames)
    cBody = bodyNames{bIdx};
    
    for pIdx=1:length(personNames)
        cPerson=personNames{pIdx};
        
        for aIdx = 1:length(actionNames)
            cAction = actionNames{aIdx};
            
            combination = [cPerson cAction cBody];
            
            tmpFR=[]; tmpbFR=[]; tmpFRorig=[];
            
            for d=1:length(UDates)
                
                dateCode = Blackrock.Helper.date2unitId(UDates{d});
                cASF = ASF(((ASF(:,4)==dateCode) & goodFiring),:);
                
                ctrialdata = Analyze.SubSelectTrials(AllTrialData, 'Date', UDates{d},'Phase',Phase,'Condition', combination);
                
                cFR= Analyze.getNeuralData(ctrialdata,cASF,timeWindow,'smooth',smoothArg);
                cFR(:,any(any(isnan(cFR),1),3),:)=[];
                
                %             cbFR = Analyze.getNeuralData(ctrialdata,cASF,[-.75 -.25]);
                %             %             cbFR(:,any(isnan(cbFR)))=[];
                %             cbFR(:,any(any(isnan(cbFR),1),3),:)=[];
                
                tmpFR = cat(3,tmpFR,cFR);
            end
            
            FR{bIdx,pIdx,aIdx} = tmpFR; clear tmpFRorig;
            %         FR{bIdx,pIdx} = tmpFR; clear tmpFR; clear tmpbFR;
        end
        
    end
    
end
