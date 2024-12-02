function SinglechanBandAnalysis(config,varargin)
%%
if strcmp(mfilename,'') || nargin==0
    opts=Analyze.Libet.Config_BasicLinearModel;
else
    opts=feval(config,varargin{:});
end


if ~exist(opts.ResultsDir,'dir'); mkdir(opts.ResultsDir); end

ConfigureSubject(opts.Subject);

for dateIDX=1:length(opts.Dates)
    TaskDate=opts.Dates{dateIDX};
    
    AllTrialData=Analyze.LoadConvertedData(opts.TaskFileName,TaskDate);
    
    % remove all trials that don't meet specified criteria
    if isfield(opts,'GlobalTrialCriteria')
        AllTrialData=Analyze.SubSelectTrials(AllTrialData,opts.GlobalTrialCriteria{:});
    end
    
    %%
    
    
    idx=1;
    for j=opts.Channels
        for i=1:length(opts.FreqBands)
            ASF(idx,:)=[j opts.FreqBands{i}];
            idx=idx+1;
        end
    end
    %%
    
    for cChanINDX=1:size(ASF,1);
        chanBand=ASF(cChanINDX,:);
        Feature=sprintf('%d-%d-%d',chanBand(1),chanBand(2),chanBand(3));
        disp(Feature)
        
        % clear variables.
        TrialResults=[];
        TrialResults_NoTime=[];
        
        if opts.saveResults
            namePrefix= sprintf('%s-%s-LFP-%d-%d-%d',opts.Subject,TaskDate,chanBand(1),chanBand(2),chanBand(3));
            
            % set save name, allowing iteration over shuffles
            if opts.shuffle==1;
                r=1;
                saveName=sprintf('Shuffle-%s-rep%d',namePrefix,r);
                while exist(fullfile(opts.ResultsDir,[saveName '.mat']),'file')
                    r=r+1;
                    saveName=sprintf('Shuffle-%s-rep%d',namePrefix,r);
                end
            else
                saveName=namePrefix;
            end
            if exist(fullfile(opts.ResultsDir,[saveName '.mat']),'file') && ~opts.overwrite
                fprintf('Already Processed - skipping : %s', fullfile(opts.ResultsDir,saveName));
                continue
            end
        end
        %%
        % Specialy analysis that is only performed once, and does not
        % depend on phases or time windows.
        if isfield(opts, 'AnalysisFCN_NoTime')
            [Results,opts]=feval(opts.AnalysisFCN_NoTime,AllTrialData,chanBand,opts);
            fn=fieldnames(Results);
            for phaseIDX=1:length(opts.Phases)
                for i=1:length(fn)
                    tmp=Results.(fn{i});
                    if isnumeric(tmp) && length(tmp)==1
                        TrialResults_NoTime{phaseIDX}.(fn{i})=tmp;
                    else
                        TrialResults_NoTime{phaseIDX}.(fn{i})=tmp;
                    end
                end
            end
        end
        
        for phaseIDX=1:length(opts.Phases)
            cPhase=opts.Phases{phaseIDX};
            windowStarts=opts.windowStarts.(cPhase);
            windowDuration=opts.windowDuration.(cPhase);
            timeWindows=[windowStarts(:) repmat(windowStarts(:),1,length(windowDuration))+repmat(windowDuration,length(windowStarts),1)];
            
            %%
            for fcnIDX=1:size(opts.AnalysisFCN,1)
                opts.getFeats=1;
                FeatInfo=feval(opts.AnalysisFCN{fcnIDX,1},AllTrialData,cPhase,chanBand,timeWindows,opts);
                opts.getFeats=0;
                
                for  timeIDX=1:size(timeWindows,1)
                    %%
                    timeWindow=timeWindows(timeIDX,:);
                    if ~iscell(opts.AnalysisFCN); opts.AnalysisFCN={opts.AnalysisFCN,''}; end
                    
                    Results=feval(opts.AnalysisFCN{fcnIDX,1},AllTrialData,cPhase,FeatInfo,timeIDX,opts);
                    
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
                TrialResults{phaseIDX}.('chanBand'){timeIDX}=chanBand(1:3);
                TrialResults{phaseIDX}.('PhaseStart')(timeIDX)=timeWindow(1);
                TrialResults{phaseIDX}.('PhaseEnd')(timeIDX)=timeWindow(2);
                end
               
                %%
            end
            %%
            optsOut=opts;
            optsOut.chanBand=chanBand;
            optsOut.windowStarts=windowStarts;
            optsOut.windowDuration=windowDuration;
            optsOut.timeWindows=timeWindows;
            optsOut.TaskDate=TaskDate;
            optsOut.Phase=opts.Phases;
            
            
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
