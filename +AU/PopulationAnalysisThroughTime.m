function PopulationAnalysisThroughTime(config,varargin)
%%
if strcmp(mfilename,'') || nargin==0
    opts=Analyze.Libet.Config_BasicLinearModel;
else
    opts=feval(config,varargin{:});
end
if ~isfield(opts,'WindowType') 
    opts.WindowType='windowStarts';    
end 

if strcmp(opts.WindowType,'windowStartsTR')
   opts.Phases=opts.PhasesTR; 
end

if ~exist(opts.ResultsDir,'dir'); mkdir(opts.ResultsDir); end

% ConfigureSubject(opts.Subject);

if ~iscell(opts.AnalysisFCN); opts.AnalysisFCN={opts.AnalysisFCN,''}; end
for fcnIDX=1:size(opts.AnalysisFCN,1)
    clear TrialResults
    for dateIDX=1:length(opts.Dates)
        clear TrialResults
        TaskDate=opts.Dates{dateIDX};
        opts.TaskDate=TaskDate;
        if opts.saveResults
            % set save name, allowing iteration over shuffles
            if iscell(TaskDate)
                saveDate=sprintf('%dDays',length(TaskDate));
            else
                saveDate=TaskDate;
            end
             if isfield(opts,'Array') && ~isempty(opts.Array)
            saveBase=sprintf('%s-%s-%s-Array%d',opts.AnalysisFCN{fcnIDX,2},opts.Subject,saveDate,opts.Array);
             else
                             saveBase=sprintf('%s-%s-%s-ArraysAll',opts.AnalysisFCN{fcnIDX,2},opts.Subject,saveDate);     
             end
             saveBase
            opts.saveName=saveBase;
            if opts.shuffle==1;
                r=1;
                saveName=sprintf('Shuffle-%s-rep%d',saveBase,r);
                while exist(fullfile(opts.ResultsDir,[saveName '.mat']),'file')
                    r=r+1;
                    saveName=sprintf('Shuffle-%s-rep%d',saveBase,r);
                end
            else
                saveName=saveBase
            end
            if exist(fullfile(opts.ResultsDir,[saveName '.mat']),'file') && ~opts.overwrite
                fprintf('Already Processed - skipping : %s', fullfile(opts.ResultsDir,saveName));
                continue
            end
        end
        
        AllTrialData=Analyze.LoadConvertedData(opts.TaskFileName,TaskDate);
        
        % remove all trials that don't meet specified criteria
        if isfield(opts,'GlobalTrialCriteria')
            AllTrialData=Analyze.SubSelectTrials(AllTrialData,opts.GlobalTrialCriteria{:});
        end
        %%
        % Get (A)ll (S)pike (F)eatures
        clear unitsAboveThresh notANan goodSNR
        ASF=Analyze.returnUniqueFieldValues(AllTrialData,opts.unitType);
        
        
        if isfield(opts,'Array') && ~isempty(opts.Array)
           ArrayIDX=ASF(:,1)==opts.Array;
        else
            ArrayIDX=logical(ones(size(ASF(:,1))));
        end
        
        if isfield(opts,'limitWFSNR') && opts.limitWFSNR
            
            %
            WFtmp=Analyze.getNeuralData(AllTrialData,ASF,'Waveform');
            clear WF
            if size(WFtmp,2)>1
                for i=1:size(WFtmp,1)
                    a=WFtmp(i,:);
                    a(cellfun(@length,a)==1)=[];
                    WF{i}=cat(2,a{:});
                end
            else
                WF=WFtmp;
            end
            [snrval,FR,goodSNR]=computeSNR(WF','RateThresh',.05,'SNRThresh',1);
            goodSNR=snrval>opts.limitWFSNR;
            fprintf('%d/%d units meet waveform SNR criteria\n',nnz(goodSNR),length(goodSNR));
            goodSNR=goodSNR(:);
            
        else
            goodSNR=double(ones(size(ASF,1),1));
            goodSNR=goodSNR(:);
        end
        
        MR=mean(Analyze.getNeuralData(AllTrialData,ASF,'MeanRate'),2);
        unitsAboveThresh=MR>opts.FRthreshold;
        fprintf('%d/%d units meet firing rate criteria\n',nnz(unitsAboveThresh),length(unitsAboveThresh));
        
        notANan=~isnan(MR);
        fprintf('%d/%d are not NaNs\n',nnz(notANan),length(notANan));
        ASF=ASF(unitsAboveThresh & notANan & goodSNR & ArrayIDX,:);
        
        %%
        
        %%
        
        for phaseIDX=1:length(opts.Phases)
            cPhase=opts.Phases{phaseIDX};
            if strcmp(opts.WindowType,'windowStarts')
                windowStarts=opts.windowStarts.(cPhase);
                windowDuration=opts.windowDuration.(cPhase);
            elseif strcmp(opts.WindowType,'windowStartsTR')
                windowStarts=opts.windowStartsTR.(cPhase);
                windowDuration=opts.windowDurationTR.(cPhase);
            else
                error(sprintf('Invalid WindowType defition: %s',opts.WindowType))
            end
            timeWindows=[windowStarts(:) repmat(windowStarts(:),1,length(windowDuration))+repmat(windowDuration,length(windowStarts),1)];
            
            %%
            for  timeIDX=1:size(timeWindows,1)
                %%
                timeWindow=timeWindows(timeIDX,:);
                
                
                
                Results=feval(opts.AnalysisFCN{fcnIDX,1},AllTrialData,cPhase,ASF,timeWindow,opts);
                if isempty(Results); break; end
                
                fn=fieldnames(Results);
                for i=1:length(fn)
                    tmp=Results.(fn{i});
                    if isempty(opts.AnalysisFCN{fcnIDX,2})
                        fullFN=[fn{i}];
                    else
                        fullFN=[opts.AnalysisFCN{fcnIDX,2} '_' fn{i}];
                    end
                    
                    if isnumeric(tmp) && length(tmp)==1
                        TrialResults{phaseIDX}.(fullFN)(timeIDX)=tmp;
                    else
                        TrialResults{phaseIDX}.(fullFN){timeIDX}=tmp;
                    end
                end
                
                TrialResults{phaseIDX}.('Phase'){timeIDX}=cPhase;
                TrialResults{phaseIDX}.('Date'){timeIDX}=TaskDate;
                TrialResults{phaseIDX}.('Unit'){timeIDX}=ASF;
                TrialResults{phaseIDX}.('PhaseStart')(timeIDX)=timeWindow(1);
                TrialResults{phaseIDX}.('PhaseEnd')(timeIDX)=timeWindow(2);
                %%
            end
            %%
            optsOut=opts;
            %             optsOut.UnitQuality=Analyze.getNeuralData(AllTrialData, ASF,'Quality');
            optsOut.UnitWaveform=Analyze.getNeuralData(AllTrialData{1}, ASF,'MeanWaveform');
            optsOut.UnitmeanRate= Analyze.getNeuralData(AllTrialData, ASF,'MeanRate');
            optsOut.ASF=ASF;
            optsOut.windowStarts=windowStarts;
            optsOut.windowDuration=windowDuration;
            optsOut.timeWindows=timeWindows;
            optsOut.TaskDate=TaskDate;
            optsOut.Phase=opts.Phases;
            optsOut.AnalysisFCN=opts.AnalysisFCN{fcnIDX,2};
            
            TrialResults{phaseIDX}.Info=optsOut;
        end
        
        
        
        
        if opts.saveResults
            % merge this and existing file if overwrite is not specified.
            %              if exist(fullfile(opts.ResultsDir,[saveName '.mat']),'file') && ~opts.overwrite
            %                 oldTrialResults=load(fullfile(opts.ResultsDir,saveName));
            %              end
            if isfield(opts, 'AnalysisFCN_NoTime')
                save(fullfile(opts.ResultsDir,saveName),'TrialResults','TrialResults_NoTime')
            else
                save(fullfile(opts.ResultsDir,saveName),'TrialResults')
            end
        end
    end
end
end
