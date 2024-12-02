% Address the basic question, when x is the salient category, what
% percentage of units are tuned to x.


function PlotAnalysisThroughTimeLFP(config,varargin)

if strcmp(mfilename,'') || nargin==0
    opts=Analyze.Libet.Config_BasicLinearModel;
else
    opts=feval(config,varargin{:});
end
%%

BaseDir=opts.ResultsDir;
TaskType=opts.TaskName;
TaskDate=opts.PlotDates;
Phases=opts.Phases;

for bandIDX=1:length(opts.FreqBands)
    cBand=opts.FreqBands{bandIDX};
    cBandTXT=sprintf('%d-%d',cBand(1),cBand(2));
    if opts.PlotShuffle==1;
        Prefix=sprintf('Shuffle-%s-', opts.PlotSubject);
    else
        Prefix=sprintf('%s-', opts.PlotSubject);
    end
    %%
    if strcmp(TaskDate,'*')
        FitList=dir(fullfile(BaseDir,[Prefix '*' cBandTXT '.mat']));
        Dates='AllDays';
    else
        if iscell(TaskDate)
            FitList=[];
            for i=1:length(TaskDate)
                FitList=[FitList; dir(fullfile(BaseDir,[Prefix '' TaskDate{i} cBandTXT '.mat']))];
                TaskDateSave{i}=[TaskDate{i}(3:end) '-'];
            end
            Dates=cat(2,TaskDateSave{:});
        else
            FitList=dir(fullfile(BaseDir,[Prefix '' TaskDate ,'*-' cBandTXT '.mat']));
            Dates=TaskDate;
        end
    end
    PopulationData=fullfile(BaseDir,'PopData'); if ~exist(PopulationData), mkdir(PopulationData); end
    
    if strfind(opts.PlotSubject,'*')
        PrefixName=strrep(Prefix,'*','All');
    else
        PrefixName=Prefix;
    end
    ResultsFile=[PrefixName Dates '-' cBandTXT];
    
    
    
    
    
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
                
                PlotData{phaseIDX}.TaskDate{unitIDX}=cTrialResults{phaseIDX}.Info.TaskDate;
                PlotData{phaseIDX}.chanBand{unitIDX}=cTrialResults{phaseIDX}.chanBand{1};
                
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
    for i=1:length(Phases)
        Dur(i)=Time{i}.WindowStops(end)-Time{i}.WindowStarts(1);
    end
    PhaseSizes=Dur/sum(Dur);
    for i=1:length(PhaseSizes); hPack{i}=PhaseSizes(i); end
    
    %%
    for plotFCNidx=1:size(opts.PlotFCN,1)
        feval(opts.PlotFCN{plotFCNidx,1})
        plt.SaveFigure(1,PopulationData,[opts.PlotFCN{plotFCNidx,2} '-' PrefixName Dates '-' cBandTXT], 'PNG','SVGI')
    end
end

