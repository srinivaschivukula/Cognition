function BehaviorAnalysis(config,varargin)
%%
if strcmp(mfilename,'') || nargin==0
    opts=Analyze.Libet.Config_BasicLinearModel;
else
    opts=feval(config,varargin{:});
end


if ~exist(opts.ResultsDirBehav,'dir'); mkdir(opts.ResultsDirBehav); end

ConfigureSubject(opts.Subject);

if ~iscell(opts.AnalysisFCN); opts.AnalysisFCN={opts.AnalysisFCN,''}; end
for fcnIDX=1:size(opts.AnalysisFCN,1)
     opts.FCN=opts.BehaviorAnalysisFCN{fcnIDX,1};
    clear TrialResults
    for dateIDX=1:length(opts.Dates)
        TaskDate=opts.Dates{dateIDX};
        opts.TaskDate=TaskDate;
        % set save name, allowing iteration over shuffles
    if iscell(TaskDate)
                saveDate=sprintf('%dDays',length(TaskDate));
            else
                saveDate=TaskDate;
            end
           
            saveBase=sprintf('%s-%s-%s',opts.AnalysisFCN{fcnIDX,2},opts.Subject,saveDate);
             opts.saveName=saveBase;       
        if opts.shuffle==1;
            r=1;
            saveName=sprintf('Shuffle-%s-rep%d',saveBase,r);
            while exist(fullfile(opts.ResultsDirBehav,[saveName '.mat']),'file')
                r=r+1;
                saveName=sprintf('Shuffle-%s-rep%d',saveBase,r);
            end
        else
            saveName=saveBase;
        end
        if exist(fullfile(opts.ResultsDirBehav,[saveName '.mat']),'file') && ~opts.overwrite
            fprintf('Already Processed - skipping : %s', fullfile(opts.ResultsDirBehav,saveName));
            continue
        end
        
        opts.saveName=saveName;
        AllTrialData=Analyze.LoadConvertedData(opts.TaskFileName,TaskDate);
        
        % remove all trials that don't meet specified criteria
        if isfield(opts,'GlobalTrialCriteria')
            AllTrialData=Analyze.SubSelectTrials(AllTrialData,opts.GlobalTrialCriteria);
        end
        %%
        
        
        for phaseIDX=1:length(opts.Phases)
            cPhase=opts.Phases{phaseIDX};
            opts.phase=cPhase;
            if ~opts.saveResults
                feval(opts.BehaviorAnalysisFCN{fcnIDX,1},AllTrialData,cPhase,opts);
            else
                Results=feval(opts.BehaviorAnalysisFCN{fcnIDX,1},AllTrialData,cPhase,opts);
                
                fn=fieldnames(Results);
                for i=1:length(fn)
                    tmp=Results.(fn{i});
                    if isempty(opts.BehaviorAnalysisFCN{fcnIDX,2})
                        fullFN=[fn{i}];
                    else
                        fullFN=[opts.BehaviorAnalysisFCN{fcnIDX,2} '_' fn{i}];
                    end
                    
                    if isnumeric(tmp) && length(tmp)==1
                        TrialResults{phaseIDX}.(fullFN)=tmp;
                    else
                        TrialResults{phaseIDX}.(fullFN)=tmp;
                    end
                end
                
                TrialResults{phaseIDX}.('Phase')=cPhase;
                TrialResults{phaseIDX}.('Date')=TaskDate;
                TrialResults{phaseIDX}.('Unit')=ASF;
                TrialResults{phaseIDX}.('PhaseStart')=timeWindow(1);
                TrialResults{phaseIDX}.('PhaseEnd')=timeWindow(2);
                
                %%
                optsOut=opts;
                optsOut.TaskDate=TaskDate;
                optsOut.Phase=opts.Phases;
                optsOut.BehaviorAnalysisFCN=opts.BehaviorAnalysisFCN{fcnIDX,2};
                
                TrialResults{phaseIDX}.Info=optsOut;
            end
        end
        
        
        
        
        if opts.saveResults
            % merge this and existing file if overwrite is not specified.
            %              if exist(fullfile(opts.ResultsDirBehav,[saveName '.mat']),'file') && ~opts.overwrite
            %                 oldTrialResults=load(fullfile(opts.ResultsDirBehav,saveName));
            %              end
            
            save(fullfile(opts.ResultsDirBehav,saveName),'TrialResults')
            
        end
    end
end
end
