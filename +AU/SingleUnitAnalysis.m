function SingleUnitAnalysis(config,varargin)
%%
% DEtermine whether to use single windows or finescale windows



% opts..(FN{i})=opts..(FN{i});

if strcmp(mfilename,'') || nargin==0
    opts=Analyze.Libet.Config_BasicLinearModel;
else
    opts=feval(config,varargin{:});
end

WindowType=Utilities.getInputField(opts,'WindowType','windowStarts',{'windowStarts','windowStartsTR'});



if ~exist(opts.ResultsDir,'dir'); mkdir(opts.ResultsDir); end

% ConfigureSubject(opts.Subject);

for dateIDX=1:length(opts.Dates)
    TaskDate=opts.Dates{dateIDX};
    opts.TaskDate=TaskDate;
    AllTrialData=Analyze.LoadConvertedData(opts.TaskFileName,TaskDate);
    
  %% Sri Chivukula Edits 
%   if strcmp(TaskDate,'20180412') | strcmp(TaskDate, '20180416') | strcmp(TaskDate, '20180425')| strcmp(TaskDate, '20180426')
%         
%         tmp = opts.Labels';
%         tmp2 = {};count = 1;
%         for tmpIDX = 1:length(tmp)
%             
%             
%             if strcmp(tmp{tmpIDX}(1),'R') | strcmp(tmp{tmpIDX}(1),'F')
%                 tmp2{count} = tmp{tmpIDX};
%                 count = count+1;
%             end
%         end
%    opts.Labels = tmp2;      
%   end
   
  %%  
    % Remove TrialData that don't have TaskLabel field (may be from video
    % task instead)
    warning('Skipping something')
    if 1==0
    goodTrialData = logical(zeros(1,length(AllTrialData)));
    for i = 1:length(AllTrialData)
        goodTrialData(i) = isfield(AllTrialData{i},'TaskLabel');
    end
    AllTrialData=AllTrialData(goodTrialData);
    end
    % remove all trials that don't meet specified criteria
    if isfield(opts,'GlobalTrialCriteria')
        AllTrialData=Analyze.SubSelectTrials(AllTrialData,opts.GlobalTrialCriteria{:});
    end
    %%
    % Get (A)ll (S)pike (F)eatures
    if isnumeric(opts.unitType)
        ASF=opts.unitType;
    else
        ASF=Analyze.returnUniqueFieldValues(AllTrialData,opts.unitType);
        %%
        if isfield(opts,'limitWFSNR') && opts.limitWFSNR
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
            [snrval,FR,goodSNR]=computeSNR(WF','RateThresh',.05,'SNRThresh',.5);
            fprintf('%d/%d units meet waveform SNR criteria\n',nnz(goodSNR),length(goodSNR));
            goodSNR=goodSNR(:);
        else
            goodSNR=double(ones(size(ASF,1),1));
            goodSNR=goodSNR(:);
        end
        %%
        if strcmpi(opts.GlobalTrialCriteria{2},'RFMap')
            MR=mean(Analyze.getNeuralData(AllTrialData(2:3),ASF,'MeanRate'),2);
        else
            MR=mean(Analyze.getNeuralData(AllTrialData,ASF,'MeanRate'),2);
        end
        unitsAboveThresh=MR>opts.FRthreshold;
        fprintf('%d/%d units meet firing rate criteria\n',nnz(unitsAboveThresh),length(unitsAboveThresh));
        
        notANan=~isnan(MR);
        fprintf('%d/%d are not NaNs\n',nnz(notANan),length(notANan));
        ASF=ASF(unitsAboveThresh & notANan & goodSNR,:);
    end
    %%
    
    for cUnitINDX=1:size(ASF,1)
        unit=ASF(cUnitINDX,:);
        Feature=sprintf('%d-%d-%d',unit(1),unit(2),unit(3));
        disp(Feature)
        
        % clear variables.
        TrialResults=[];
        TrialResults_NoTime=[];
        
        if opts.saveResults
            % set save name, allowing iteration over shuffles
            if opts.shuffle==1
                r=1;
                saveName=sprintf('Shuffle-%s-%s-%d-%d-%d-rep%d',opts.Subject,TaskDate,unit(1),unit(2),unit(3),r);
                while exist(fullfile(opts.ResultsDir,[saveName '.mat']),'file')
                    r=r+1;
                    saveName=sprintf('Shuffle-%s-%s-%d-%d-%d-rep%d',opts.Subject,TaskDate,unit(1),unit(2),unit(3),r);
                end
            else
                saveName=sprintf('%s-%s-%d-%d-%d',opts.Subject,TaskDate,unit(1),unit(2),unit(3));
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
      
            for phaseIDX=1:length(opts.Phases)
                      [Results]=feval(opts.AnalysisFCN_NoTime{1},AllTrialData,opts.Phases{phaseIDX},unit,opts);
            fn=fieldnames(Results);
                
                for i=1:length(fn)
                    tmp=Results.(fn{i});
                    if isnumeric(tmp) && length(tmp)==1
                        TrialResults_NoTime{phaseIDX}.(fn{i})=tmp;
                    else
                        TrialResults_NoTime{phaseIDX}.(fn{i})={tmp};
                    end
                end
                %%
                %%%%%%%%%%%%
                  optsOut_NoTime=opts;
            TrialResults_NoTime{phaseIDX}.UnitQuality=Analyze.getNeuralData(AllTrialData, unit,'Quality');
            TrialResults_NoTime{phaseIDX}.UnitWaveform=Analyze.getNeuralData(AllTrialData{1}, unit,'MeanWaveform');
            TrialResults_NoTime{phaseIDX}.UnitmeanRate= mean(Analyze.getNeuralData(AllTrialData, unit,'MeanRate'));

            Waveform=Analyze.getNeuralData(AllTrialData, unit,'Waveform');
            Waveform=double(cat(2,Waveform{:}));
            try
            [~,idx]=max(mean(Waveform')); p1=idx+(-1:1);b=(0:5)+1;
            DP1=max(-(mean(nanmedian(Waveform(b,:)'))-nanmedian(Waveform(p1,:)'))./sqrt(mean(nanstd(Waveform(b,:)')).^2+nanstd(Waveform(p1,:)').^2)/2);
            TrialResults_NoTime{phaseIDX}.unitWFDP=DP1;
            catch
                optsOut_NoTime.unitWFDP=nan;
            end

            %             b=(0:5)+1; p1=(0:5)+21;p2=(0:5)+29;
            %             DP2=max(-(mean(nanmedian(Waveform(b,:)'))-nanmedian(Waveform(p2,:)'))./sqrt(mean(nanstd(Waveform(b,:)')).^2+nanstd(Waveform(p2,:)').^2)/2);
            %             optsOut.unitWFDP=max([DP1 DP2]);
            %%
            TrialResults_NoTime{phaseIDX}.unit={unit};
            TrialResults_NoTime{phaseIDX}.timeWindows={opts.timeWindow};
            TrialResults_NoTime{phaseIDX}.TaskDate={TaskDate};
            TrialResults_NoTime{phaseIDX}.Phase=opts.Phases;
            
            
            TrialResults_NoTime{phaseIDX}.Info=optsOut_NoTime;

                
                %%%%%%%
                %%
            end
        end
        
        for phaseIDX=1:length(opts.Phases)
            cPhase=opts.Phases{phaseIDX};
            switch WindowType
                case 'windowStarts'
                    windowStarts=opts.windowStarts.(cPhase);
                    windowDuration=opts.windowDuration.(cPhase);
                    timeWindows=[windowStarts(:) repmat(windowStarts(:),1,length(windowDuration))+repmat(windowDuration,length(windowStarts),1)];
                case 'windowStartsTR';
                    windowStarts=opts.windowStartsTR.(cPhase);
                    windowDuration=opts.windowDurationTR.(cPhase);
                    timeWindows=[windowStarts(:) repmat(windowStarts(:),1,length(windowDuration))+repmat(windowDuration,length(windowStarts),1)];
            end
            
            
            %%
            for  timeIDX=1:size(timeWindows,1)
                %%
                timeWindow=timeWindows(timeIDX,:);
                if ~iscell(opts.AnalysisFCN); opts.AnalysisFCN={opts.AnalysisFCN,''}; end
                
                for fcnIDX=1:size(opts.AnalysisFCN,1)
                    opts.fName=opts.AnalysisFCN{fcnIDX,2};
                    try
                        if strcmp(opts.AnalysisName,'WAStim') || strcmp(opts.AnalysisName,'StimProfiles') ||...
                                strcmpi(opts.AnalysisName,'StimProfilesFourObj') || strcmpi(opts.AnalysisName,'StimProfilesMultipleObj')
                            Results=feval(opts.AnalysisFCN{fcnIDX,1},AllTrialData,cPhase,unit,timeWindow,...
                                opts.BaselineWindow,opts);
                        elseif strcmpi(opts.AnalysisName,'RFMap') 
                            Results=feval(opts.AnalysisFCN{fcnIDX,1},AllTrialData,cPhase,unit,timeWindow,opts.BaselineWindow,opts,TaskDate);
                        elseif strcmpi(opts.AnalysisName,'Direction')
                            Results=feval(opts.AnalysisFCN{fcnIDX,1},AllTrialData,cPhase,unit,timeWindow,opts.BaselineWindow,opts,TaskDate);
                        else
                            Results=feval(opts.AnalysisFCN{fcnIDX,1},AllTrialData,cPhase,unit,timeWindow,opts);
                        end
                    catch
                        keyboard
                    end
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
                end
                TrialResults{phaseIDX}.('Phase'){timeIDX}=cPhase;
                TrialResults{phaseIDX}.('Date'){timeIDX}=TaskDate;
                TrialResults{phaseIDX}.('Unit'){timeIDX}=unit;
                TrialResults{phaseIDX}.('PhaseStart')(timeIDX)=timeWindow(1);
                TrialResults{phaseIDX}.('PhaseEnd')(timeIDX)=timeWindow(2);
                %%
            end
            %%
            optsOut=opts;
            optsOut.UnitQuality=Analyze.getNeuralData(AllTrialData, unit,'Quality');
            optsOut.UnitWaveform=Analyze.getNeuralData(AllTrialData{1}, unit,'MeanWaveform');
            optsOut.UnitmeanRate= Analyze.getNeuralData(AllTrialData, unit,'MeanRate');
%             optsOut.IsolDist= Analyze.getNeuralData(AllTrialData, unit,'CV2Pooled');

            %%
            if strcmpi(opts.AnalysisName,'RFMap') || strcmpi(opts.AnalysisName,'Direction')
                Waveform=Analyze.getNeuralData(AllTrialData(2:end), unit,'Waveform');
            else
                Waveform=Analyze.getNeuralData(AllTrialData, unit,'Waveform');
            end                        
            Waveform=double(cat(2,Waveform{:}));
            try
            [~,idx]=max(mean(Waveform')); p1=idx+(-1:1);b=(0:5)+1;
            DP1=max(-(mean(nanmedian(Waveform(b,:)'))-nanmedian(Waveform(p1,:)'))./sqrt(mean(nanstd(Waveform(b,:)')).^2+nanstd(Waveform(p1,:)').^2)/2);
            optsOut.unitWFDP=DP1;
            catch
                optsOut.unitWFDP=nan;
            end

            %             b=(0:5)+1; p1=(0:5)+21;p2=(0:5)+29;
            %             DP2=max(-(mean(nanmedian(Waveform(b,:)'))-nanmedian(Waveform(p2,:)'))./sqrt(mean(nanstd(Waveform(b,:)')).^2+nanstd(Waveform(p2,:)').^2)/2);
            %             optsOut.unitWFDP=max([DP1 DP2]);
            %%
            optsOut.unit=unit;
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
            if (isfield(opts, 'AnalysisFCN_NoTime') & isempty(opts.AnalysisFCN))
                TrialResults=TrialResults_NoTime;
                                save(fullfile(opts.ResultsDir,saveName),'TrialResults')

            elseif isfield(opts, 'AnalysisFCN_NoTime')
                save(fullfile(opts.ResultsDir,saveName),'TrialResults','TrialResults_NoTime')
            else                
                save(fullfile(opts.ResultsDir,saveName),'TrialResults')
            end
        end
    end
end
