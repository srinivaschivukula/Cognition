% Address the basic question, when x is the salient category, what
% percentage of units are tuned to x.


function PlotAnalysisThroughTime(config,varargin)

[varargin,DateName]   = Utilities.ProcVarargin(varargin,'DateName',[]);


if strcmp(mfilename,'') || nargin==0
    opts=Analyze.Libet.Config_BasicLinearModel;
elseif isa(config,'function_handle') || isstr(config)
    opts=feval(config,varargin{:});
else
    opts=config;
end

if ~isfield(opts,'PlotSaveType')
    PlotSaveType={'PNG','PDF'};
else
    PlotSaveType=opts.PlotSaveType;
end
%%

BaseDir=opts.ResultsDir;
TaskType=opts.TaskName;
TaskDate=opts.PlotDates;
Phases=opts.Phases;


if isfield(opts,'PlotShuffle') && opts.PlotShuffle==1;
    Prefix=sprintf('Shuffle-%s-', opts.PlotSubject);
else
    Prefix=sprintf('%s-', opts.PlotSubject);
end

if isfield(opts,'ArrayPlot') && ~isempty(opts.ArrayPlot)
%     array=['-Array' num2str(opts.ArrayPlot)];
    array=['-' num2str(opts.ArrayPlot)];
else
%     array='-ArraysAll';
    array='';
end
%%
fprintf('Searching %s \n',BaseDir)
if strcmp(TaskDate,'*')
    FitList=dir(fullfile(BaseDir,[Prefix '*' array '-*.mat']));
    Dates='AllDays';
else
    if iscell(TaskDate)
        FitList=[];
        for i=1:length(TaskDate)
            FitList=[FitList; dir(fullfile(BaseDir,[Prefix '' TaskDate{i} array '-*.mat']))];
            TaskDateSave{i}=[TaskDate{i}(3:end) '-'];
        end
        Dates=cat(2,TaskDateSave{:});
    else
        FitList=dir(fullfile(BaseDir,[Prefix '' TaskDate , array,'*.mat']));
        Dates=TaskDate;
    end
end
PopulationData=fullfile(BaseDir,'PopData'); if ~exist(PopulationData), mkdir(PopulationData); end

if strfind(opts.PlotSubject,'*')
    PrefixName=strrep(Prefix,'*','All');
else
    PrefixName=Prefix;
end

if isempty(DateName)
    ResultsFile=[PrefixName Dates array];
else
    ResultsFile=[PrefixName DateName  array];
    
end





%% Gather data
if exist(fullfile(PopulationData,[ResultsFile '.mat']),'file') && ~opts.PopulationOverwrite
    fprintf('Population Data Already Generated: loading...');
    load(fullfile(PopulationData,ResultsFile))
else
    % Note Out is nPhases x nTime bins
    idx=1;
    fprintf('Loading files similar to %s \n',  (FitList(1).name))
    for unitIDX=1:length(FitList)
        %     for unitIDX=1:2
        clear TrialResults_NoTime TrialResults
        load(fullfile(BaseDir,FitList(unitIDX).name));
        
        for phaseIDX=1:length(Phases)
            cTrialResults=Analyze.SubSelectTrials(TrialResults,'Phase',Phases(phaseIDX));
            
            
            %% get output field names
            if ~isfield(opts,'PlotFields') || isempty(opts.PlotFields)
                allFields=fieldnames(TrialResults{1});
                idx2remove=ismember(allFields,{'Phase','Date','Unit','PhaseStart','PhaseEnd','Info'});
                allFields(idx2remove)=[];
                opts.PlotFields=allFields;
            end
            %%
            for dataIDX=1:length(opts.PlotFields)
                PlotData{phaseIDX}.(opts.PlotFields{dataIDX})(:,unitIDX)=Analyze.returnFieldValues(cTrialResults,opts.PlotFields{dataIDX});
            end
            
            if exist('TrialResults_NoTime')
                PlotFields_NoTime=fieldnames(TrialResults_NoTime{1});
                for dataIDX=1:length(PlotFields_NoTime)
                    PlotData{phaseIDX}.(PlotFields_NoTime{dataIDX})(:,unitIDX)=Analyze.returnFieldValues(TrialResults_NoTime(1),PlotFields_NoTime{dataIDX});
                end
            end
            %%
            try
                PlotData{phaseIDX}.UnitQuality(unitIDX)=cTrialResults{phaseIDX}.Info.UnitQuality;
                PlotData{phaseIDX}.UnitWaveform(:,unitIDX)=cTrialResults{phaseIDX}.Info.UnitWaveform{1};
                PlotData{phaseIDX}.TaskDate{unitIDX}=cTrialResults{phaseIDX}.Info.TaskDate;
                PlotData{phaseIDX}.UnitmeanRate(unitIDX)=cTrialResults{phaseIDX}.Info.UnitmeanRate(1);
                PlotData{phaseIDX}.unitWFDP(unitIDX)=cTrialResults{phaseIDX}.Info.unitWFDP;
                PlotData{phaseIDX}.Unit{unitIDX}=cTrialResults{phaseIDX}.Unit{1};
                PlotData{phaseIDX}.Array{unitIDX}=cTrialResults{phaseIDX}.Unit{1}(1);
            catch
            end
            PlotData{phaseIDX}.Phase{1,unitIDX}=Phases{phaseIDX};
            %%
            %     for i=1:length(FitResults)
            %         Population(unitIDX).CoefEst(i,:)=FitResults(i).Coefficients.Estimate;
            %         Population(unitIDX).CoefSE(i,:)=FitResults(i).Coefficients.SE;
            %         Population(unitIDX).CoeftStat(i,:)=FitResults(i).Coefficients.tStat;
            %         Population(unitIDX).CoefpValue(i,:)=FitResults(i).Coefficients.pValue;
            %     end
            Time{phaseIDX}.WindowStarts=Analyze.returnFieldValues(cTrialResults,'PhaseStart');
            Time{phaseIDX}.WindowStops=Analyze.returnFieldValues(cTrialResults,'PhaseEnd');
            Time{phaseIDX}.WindowMid=[Time{phaseIDX}.WindowStarts+Time{phaseIDX}.WindowStops]/2;
        end
        
    end
    save(fullfile(PopulationData,ResultsFile),'PlotData','Time')
end


%%
try
for i=1:length(Phases)
    Dur(i)=Time{i}.WindowStops(end)-Time{i}.WindowStarts(1);
end
PhaseSizes=Dur/sum(Dur);
for i=1:length(PhaseSizes); hPack{i}=PhaseSizes(i); end

catch
   warning('WindowStops and WindowStarts not defined. do your own damn sizing.') 
end

%% If option enabled, only look at well tuned units
if isfield(opts,'UseSparse') && opts.UseSparse
    PlotData = Analyze.SubSelectTrials(PlotData,'>=testdPBW',2);
    %     dPBW = Analyze.returnFieldValues(PlotData,'testdPBW');
    %     cdf = cumsum(dPBW);
    %     figure; histogram(dPBW,'Normalization','cumcount');
end
% figure; histogram(dP);

%%
if isfield(opts,'unitWFDPThresh') && ~isempty(opts.unitWFDPThresh)
    if isnan(opts.unitWFDPThresh)
        unitWFDP=Analyze.returnFieldValues(PlotData(1),'unitWFDP');
        figure; histogram(unitWFDP,round(length(unitWFDP)/10))
        opts.unitWFDPThresh=prctile(unitWFDP,50);
        title('unitWFDP')
        plt.vline(opts.unitWFDPThresh);
    end
    PlotData=Analyze.SubSelectTrials(PlotData,'>=unitWFDP',opts.unitWFDPThresh);
    
    PrefixName=[PrefixName sprintf('WFDP%0.1f', opts.unitWFDPThresh) ,'-'];
end
disp(1)
if 1==0
    %%
SortOpts.Experiment=opts.TaskName;
                SortOpts.TaskName=opts.TaskFileName;
                SortOpts.Subject=opts.Subject;
                SortOpts.TaskDates=opts.Dates;
        
        PSS=AU.ProcSortStats(SortOpts);
        %%
        fn=fieldnames(PlotData);
    
        cTaskDate=Analyze.returnFieldValues(PlotData.(fn{1}),'TaskDate');
        cUnit=Analyze.returnFieldValues(PlotData.(fn{1}),'Unit');
        cTaskDate=Blackrock.Helper.date2unitId(cTaskDate);
        unitdateID=[cell2mat(cUnit),cTaskDate];
%         tmp=PSS.returnValidInds('PeakSNR','>',50,unit,DateID,'ValueType','Percentile');
%         [A,B,C]=PSS.returnSortStats('PeakSNR',unitdateID);
end

for plotFCNidx=1:size(opts.PlotFCN,1)
    try
        pltFCNArgs=opts.PlotFCN{plotFCNidx,3};
        feval(opts.PlotFCN{plotFCNidx,1})
        if opts.PlotResults && ~isempty(opts.PlotFCN{plotFCNidx,2})
            
            
            if ~exist('FigHandels') || length(FigHandels)==1 || isempty(FigHandels)
                plt.SaveFigure(1,PopulationData,[opts.PlotFCN{plotFCNidx,2} '-' PrefixName Dates array], PlotSaveType{:})
            else
                
                for figIDX=1:length(FigHandels)
                    figure(FigHandels(figIDX))
                    plt.SaveFigure(1,PopulationData,[opts.PlotFCN{plotFCNidx,2}{figIDX} '-' PrefixName Dates array], PlotSaveType{:})
                end
            end
        end
        
    catch
        close
        lasterr
        warning(['skipping' [opts.PlotFCN{plotFCNidx,2} '-' PrefixName Dates array]])
    end
end
end

