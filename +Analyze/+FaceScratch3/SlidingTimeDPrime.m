% Address the basic question, when x is the salient category, what
% percentage of units are tuned to x.


function SlidingTimeDPrime(config,varargin)

[varargin,DateName]   = Utilities.ProcVarargin(varargin,'DateName',[]);
[varargin,Pmethod]   = Utilities.ProcVarargin(varargin,'Pmethod','fdr');
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');
[varargin,DatesAll]   = Utilities.ProcVarargin(varargin,'Dates',{});
[varargin,trialdata]   = Utilities.ProcVarargin(varargin,'AllTrialData', {});
[varargin,ImagTrialData]   = Utilities.ProcVarargin(varargin,'ImagTrialData', {});
[varargin,FigSize] = Utilities.ProcVarargin(varargin,'FigSize',[]);
[varargin,SmoothingKernel]=Utilities.ProcVarargin(varargin,'SmoothingKernel', struct('mode','MinJerk','value',.5));
[varargin,TimeStep]   = Utilities.ProcVarargin(varargin,'TimeStep',0.1);
[varargin,StepSize]   = Utilities.ProcVarargin(varargin,'StepSize',0.3);
[varargin,SubtractBaseline]   = Utilities.ProcVarargin(varargin,'SubtractBaseline',0);
[varargin,NumSamples]   = Utilities.ProcVarargin(varargin,'NumSamples',15);



if isa(config,'function_handle') || isstr(config)
    opts=feval(config,varargin{:});
else
    opts=config;
end

if ~isfield(opts,'PlotSaveType')
    PlotSaveType={'PDF','PNG','PDF'};
else
    PlotSaveType=opts.PlotSaveType;
end
%%
TaskLabel = opts.GlobalTrialCriteria{2};
BaseDir=opts.ResultsDir;
TaskType=opts.TaskName;
TaskDate=opts.PlotDates;
Phases=opts.Phases;

Prefix=sprintf('%s-', opts.PlotSubject);

array='';

%%
% fprintf('Searching %s \n',BaseDir)
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
load(fullfile(PopulationData,ResultsFile))


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
end

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

plotFCNidx=1;

pltFCNArgs=opts.PlotFCN{plotFCNidx,3};

[~,Labels]   = Utilities.ProcVarargin(pltFCNArgs,'Labels','');
[~,PThresh]   = Utilities.ProcVarargin(pltFCNArgs,'PThresh',0.05);
SavePre=[opts.PlotFCN{plotFCNidx,2} '-' PrefixName Dates array];

if isfield(opts,'UseSparse') && opts.UseSparse
    SavePre=[opts.PlotFCN{plotFCNidx,2} '-' PrefixName 'sparse-' Dates];
else
    SavePre=[opts.PlotFCN{plotFCNidx,2} '-' PrefixName Dates];
end
SavePre;
[~,Phase]   = Utilities.ProcVarargin(pltFCNArgs,'Phase','Go');
PlotData=Analyze.SubSelectTrials(PlotData,'Phase',Phase);

ActionPvals=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'Pvals'));
ActionCoef=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'Coef'));
NC=size(ActionPvals,2);

pAnova=Analyze.returnFieldValues(PlotData,'pAnova');
[H,effAlphaVal,adj_p]=Utilities.MultipleComparisonsCorrection(pAnova,'method','fdr');
[Hhb,effAlphaValHB]=Utilities.MultipleComparisonsCorrection(pAnova,'method','holm-bonferroni');
R2=Analyze.returnFieldValues(PlotData,'R2adj');

dPrime=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'ComdPrime'));
AUC=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'ComAUC'));
AUCvAll=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'AUCvAll'));


ModelP=Analyze.returnFieldValues(PlotData,'pFull');
PopTuned=ModelP<PThresh;

%% Plotting stuff

pVals=cat(2,PlotData{1}.Pvals{:})';
if strcmp(Pmethod, 'holm-bonferroni')
    [Hhb,effAlphaValHB]=Utilities.MultipleComparisonsCorrection(pVals,'method',Pmethod);
    PThresh = effAlphaValHB;
elseif strcmp(Pmethod, 'fdr')
    [H,effAlphaVal,adj_p]=Utilities.MultipleComparisonsCorrection(pVals,'method',Pmethod);
    PThresh = effAlphaVal;
else
    [H,effAlphaVal,adj_p]=Utilities.MultipleComparisonsCorrection(pVals,'method',Pmethod);
    PThresh = min([0.001 effAlphaVal]);
end

ActionTuned=ActionPvals<PThresh;

%% Calculating the time-resolved distance measures for 100 randomly sampled units, 100 times

datafile = fullfile(env.get('results'),'FaceScratch','SUAnal',['XImagery2' '-Go'],'PopData','p1-AllDays.mat');
tmpData = load(datafile); % Loads PlotData
ImagData{1} = tmpData.PlotData{1};

%% this section is for felt vs imagined data
for sampIdx=1:NumSamples
    
    %% first do all the felt/imagined stuff
    
     TimeWindows = [-2.5:TimeStep:4.5];
%     TimeWindows = [-4.5:TimeStep:5.5];
    %     TimeWindows = [-1.5:TimeStep:7.5];
    %     TimeCenters = TimeWindows+TimeStep;
    TimeCenters=[];
    for timeIdx = 1:length(TimeWindows)
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
    end
    
    %% this portion is for felt vs imagined data
    AllImagUnits = vertcat(ImagData{1}.Unit{:});
        PVals = horzcat(ImagData{1}.Pvals{:})';
    H = Utilities.MultipleComparisonsCorrection(PVals,'method','fdr');
    AllImagUnits=AllImagUnits(any(H(:,[2,3,5,6]),2),:);
%         randSampImag = randsample(length(AllImagUnits),100);
    
    cUnitsImag = AllImagUnits(randSampImag,:);
    
    condNamesImag = {'ImagCheek','ImagShoulder','CheekR','ShoulderR'};
    
    FeltTrialDataFI_cheek = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go', 'Condition', {'CheekR'});
    FeltTrialDataFI_shoulder = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go', 'Condition', {'ShoulderR'});
    
    ImaginedTrialData_cheek = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go', 'Condition', {'ImagCheek'});
    ImaginedTrialData_shoulder = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go', 'Condition', {'ImagShoulder'});
    
    BaselineDataFI_cheekFelt = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Delay','Condition',{'CheekR'});
    BaselineDataFI_shoulderFelt = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Delay','Condition',{'ShoulderR'});
    BaselineDataFI_cheekImag = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go','Condition',{'ImagCheek'});
    BaselineDataFI_shoulderImag = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go','Condition',{'ImagShoulder'});
    BaselineDataFI_null = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go','Condition',{'XX'});
    
    %% do the sliding time calculation of distances/dprimes
    
    for timeIdx = 1:length(TimeWindows)
        
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        
        %% get the felt/imag data
        
        feltFRFI = [squeeze(Analyze.getNeuralData(FeltTrialDataFI_cheek, cUnitsImag, cTW)); squeeze(Analyze.getNeuralData(FeltTrialDataFI_shoulder, cUnitsImag, cTW))];       
        imagFRFI = [squeeze(Analyze.getNeuralData(ImaginedTrialData_cheek, cUnitsImag, cTW));squeeze(Analyze.getNeuralData(ImaginedTrialData_shoulder, cUnitsImag, cTW))];
        
%         if cTW(1)>-1 || cTW(2)<-2
% %                         baseFRFI = [squeeze(Analyze.getNeuralData(BaselineDataFI_cheekFelt, cUnitsImag, [-2 -1]));squeeze(Analyze.getNeuralData(BaselineDataFI_shoulderFelt, cUnitsImag, [-2 -1]));squeeze(Analyze.getNeuralData(BaselineDataFI_cheekImag, cUnitsImag, [-2 -1]));squeeze(Analyze.getNeuralData(BaselineDataFI_shoulderImag, cUnitsImag, [-2 -1]))];
%                         baseFRFI = [squeeze(Analyze.getNeuralData(BaselineDataFI_cheekFelt, cUnitsImag, [-2 0]));squeeze(Analyze.getNeuralData(BaselineDataFI_shoulderFelt, cUnitsImag, [-2 0]))];
% %             baseFRFI = squeeze(Analyze.getNeuralData(BaselineDataFI_null, cUnitsImag, [0 3]));
%         else
%                         baseFRFI = [squeeze(Analyze.getNeuralData(BaselineDataFI_cheekFelt, cUnitsImag, [-2 0]));squeeze(Analyze.getNeuralData(BaselineDataFI_shoulderFelt, cUnitsImag, [-2 0]))];
% %             baseFRFI = squeeze(Analyze.getNeuralData(BaselineDataFI_null, cUnitsImag, [0 3]));
%         end
        
        if cTW(1)>0
            %                         baseFRFI = [squeeze(Analyze.getNeuralData(BaselineDataFI_cheekFelt, cUnitsImag, [-2 -1]));squeeze(Analyze.getNeuralData(BaselineDataFI_shoulderFelt, cUnitsImag, [-2 -1]));squeeze(Analyze.getNeuralData(BaselineDataFI_cheekImag, cUnitsImag, [-2 -1]));squeeze(Analyze.getNeuralData(BaselineDataFI_shoulderImag, cUnitsImag, [-2 -1]))];
            baseFRFI = [squeeze(Analyze.getNeuralData(BaselineDataFI_cheekFelt, cUnitsImag, [-1 0]));squeeze(Analyze.getNeuralData(BaselineDataFI_shoulderFelt, cUnitsImag, [-1 0]))];
            %             baseFRFI = squeeze(Analyze.getNeuralData(BaselineDataFI_null, cUnitsImag, [0 3]));
        else
            baseFRFI = [squeeze(Analyze.getNeuralData(BaselineDataFI_cheekFelt, cUnitsImag, cTW));squeeze(Analyze.getNeuralData(BaselineDataFI_shoulderFelt, cUnitsImag, cTW))];
            %             baseFRFI = squeeze(Analyze.getNeuralData(BaselineDataFI_null, cUnitsImag, [0 3]));
        end
        
        newFeltFRFI =[]; newImagFRFI = []; newBaseFRFI=[];
        for i=1:size(feltFRFI,2)
            tmp = feltFRFI(:,i);
            newFeltFRFI(:,i) = tmp(~isnan(tmp));
            
            tmp = imagFRFI(:,i);
            newImagFRFI(:,i) = tmp(~isnan(tmp));
            
            tmp = baseFRFI(:,i);
            newBaseFRFI(:,i) = tmp(~isnan(tmp));
        end
        feltFRFI = newFeltFRFI; imagFRFI=newImagFRFI; baseFRFI=newBaseFRFI;
        newFeltFRFI=[]; newObsFRFI=[]; newBaseFRFI=[];
        
        %% calculate all the distances/dprimes
        
        ResultsFeltFI = Analyze.FitPLSLM([feltFRFI;baseFRFI],[ones(size(feltFRFI,1),1);zeros(size(baseFRFI,1),1)],'CrossValidate',1);
        ResultsObsFI = Analyze.FitPLSLM([imagFRFI;baseFRFI],[ones(size(imagFRFI,1),1);zeros(size(baseFRFI,1),1)],'CrossValidate',1);
        
        DataFeltFI{sampIdx}{timeIdx} = ResultsFeltFI.Dist;
        DataObsFI{sampIdx}{timeIdx} = ResultsObsFI.Dist;
        
        DataFeltDPFI{sampIdx}{timeIdx} = ResultsFeltFI.dP;
        DataObsDPFI{sampIdx}{timeIdx} = ResultsObsFI.dP;
        
        
    end
    
    
    %% this section is for felt vs observed data
    
    
    TimeWindows = [-1.5:TimeStep:7.5];
    %     TimeCenters = TimeWindows+TimeStep;
    TimeCenters=[];
    for timeIdx = 1:length(TimeWindows)
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
    end
    
    %% this section is for felt vs observed data
    
    AllUnits = vertcat(PlotData{1}.Unit{:});
    randSamp = randsample(length(AllUnits),100);
    
    cUnits = AllUnits(randSamp,:);
    
    condNames = {'NancyCheek','NancyShoulder','TysonCheek','TysonShoulder'};
    
    FeltTrialData_cheek = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'NancyCheek'});
    FeltTrialData_shoulder = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'NancyShoulder'});
    
    ObsTrialData_cheek = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'TysonCheek'});
    ObsTrialData_shoulder = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'TysonShoulder'});
    
    BaselineData_cheekFelt = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go','Condition',{'NancyCheek'});
    BaselineData_shoulderFelt = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go','Condition',{'NancyShoulder'});
    BaselineData_cheekObs = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go','Condition',{'TysonCheek'});
    BaselineData_shoulderObs = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go','Condition',{'TysonShoulder'});
    
    
    for timeIdx = 1:length(TimeWindows)
        %     for timeIdx=1
        
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        
        %% get the felt/obs data
        
        feltFR = [squeeze(Analyze.getNeuralData(FeltTrialData_cheek, cUnits, cTW));squeeze(Analyze.getNeuralData(FeltTrialData_shoulder, cUnits, cTW))];
        obsFR = [squeeze(Analyze.getNeuralData(ObsTrialData_cheek, cUnits, cTW));squeeze(Analyze.getNeuralData(ObsTrialData_shoulder, cUnits, cTW))];
        
        if cTW(1)>0
            baseFR = [squeeze(Analyze.getNeuralData(BaselineData_cheekFelt, cUnits, [-1 -0]));squeeze(Analyze.getNeuralData(BaselineData_shoulderFelt, cUnits, [-1 -0]));squeeze(Analyze.getNeuralData(BaselineData_cheekObs, cUnits, [-1 -0]));squeeze(Analyze.getNeuralData(BaselineData_shoulderObs, cUnits, [-1 -0]))];
        else
            baseFR = [squeeze(Analyze.getNeuralData(BaselineData_cheekFelt, cUnits, cTW));squeeze(Analyze.getNeuralData(BaselineData_shoulderFelt, cUnits, cTW));squeeze(Analyze.getNeuralData(BaselineData_cheekObs, cUnits, cTW));squeeze(Analyze.getNeuralData(BaselineData_shoulderObs, cUnits, cTW))];
        end
        
        newFeltFR =[]; newObsFR = []; newBaseFR=[];
        for i=1:size(feltFR,2)
            tmp = feltFR(:,i);
            newFeltFR(:,i) = tmp(~isnan(tmp));
            
            tmp = obsFR(:,i);
            newObsFR(:,i) = tmp(~isnan(tmp));
            
            tmp = baseFR(:,i);
            newBaseFR(:,i) = tmp(~isnan(tmp));
        end
        feltFR = newFeltFR; obsFR=newObsFR; baseFR=newBaseFR;
        newFeltFR=[]; newObsFR=[]; newBaseFR=[];
        
        
        ResultsFelt = Analyze.FitPLSLM([feltFR;baseFR],[ones(size(feltFR,1),1);zeros(size(baseFR,1),1)],'CrossValidate',1);
        ResultsObs = Analyze.FitPLSLM([obsFR;baseFR],[ones(size(obsFR,1),1);zeros(size(baseFR,1),1)],'CrossValidate',1);
        
        DataFelt{sampIdx}{timeIdx} = ResultsFelt.Dist;
        DataObs{sampIdx}{timeIdx} = ResultsObs.Dist;
        
        DataFeltDP{sampIdx}{timeIdx} = ResultsFelt.dP;
        DataObsDP{sampIdx}{timeIdx} = ResultsObs.dP;
        
    end
    
end

%% post processing of data

Felt=[]; Obs=[]; FeltFI=[]; Imag=[];

for i=1:length(DataFelt)
    tmp = [DataFelt{i}{:}];
    Felt = cat(1,Felt,tmp);
    
    tmp = [DataObs{i}{:}];
    Obs = cat(1,Obs,tmp);
    
end

for i=1:length(DataFeltFI)
    
    tmp = [DataFeltFI{i}{:}];
    FeltFI = cat(1,FeltFI,tmp);
    
    tmp = [DataObsFI{i}{:}];
    Imag = cat(1,Imag,tmp);
    
end

if SubtractBaseline
    
    for i=1:size(Felt,2)
        
                TimeWindows = [-1.5:TimeStep:7.5];
    %     TimeCenters = TimeWindows+TimeStep;
    TimeCenters=[];
    for timeIdx = 1:length(TimeWindows)
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
    end
    BaseIdx = TimeCenters<=0;
        
        tmp = mean(Felt(:,BaseIdx),2);
        Felt = Felt-tmp;
        
        tmp = mean(Obs(:,BaseIdx),2);
        Obs = Obs-tmp;
        
    end
    
    for i=1:size(FeltFI,2)
        
                TimeWindows = [-4.5:TimeStep:5.5];
    %     TimeWindows = [-1.5:TimeStep:7.5];
    %     TimeCenters = TimeWindows+TimeStep;
    TimeCenters=[];
    for timeIdx = 1:length(TimeWindows)
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
    end

            BaseIdx = TimeCenters<=0 & TimeCenters>=-3;

        tmp = mean(FeltFI(:,BaseIdx),2);
        FeltFI = FeltFI-tmp;
        
        tmp = mean(Imag(:,BaseIdx),2);
        Imag = Imag-tmp;
        
    end
    
end

PlotData = {{Felt},{Obs},{FeltFI},{Imag}};
%% Plotting by stimulus

StyleArgs = {'FontName','Arial','FontSize',9,'FontWeight','bold'};

figH = plt.fig('units','inches','width',14,'height',8,'font','Arial','fontsize',9);
figure(figH);
clf; pnl = panel(); pnl.margin=15; pnl.pack(3,7);

condNames = {'Felt','Observed','Felt','Imagined'};

for condIdx=1:length(condNames)
    
    cData = PlotData{condIdx};
    
    clr=[];
    clr(1,:) = repmat([237 0 38]./255,1,1);
    clr(2,:) = repmat([31 120 180]./255,1,1);
    clr(3,:) = repmat([237 0 38]./255,1,1);
    clr(4,:) = repmat([1 108 89]./255,1,1);
    
    cLabel = condNames{condIdx};
    
    if condIdx==1 || condIdx==2
    
        TimeWindows = [-1.5:TimeStep:7.5];
    %     TimeCenters = TimeWindows+TimeStep;
    TimeCenters=[];
    for timeIdx = 1:length(TimeWindows)
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
    end
    
    elseif condIdx==3 || condIdx==4
        TimeWindows = [-4.5:TimeStep:5.5];
    %     TimeWindows = [-1.5:TimeStep:7.5];
    %     TimeCenters = TimeWindows+TimeStep;
    TimeCenters=[];
    for timeIdx = 1:length(TimeWindows)
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
    end
    end
    % plot the discriminability of touch from baseline
    pnl(1,condIdx).select()
    
    Analyze.plotEventRelatedAverage(cData,cLabel,'TimeVec',TimeCenters,'usePercentile',...
        'Colors',clr(condIdx,:));
    
    ax1 = gca;
    set(ax1,'box','off',StyleArgs{:});
    ax1.XRuler.Axle.LineStyle = 'none';
    set(gca,'TickDir','out');
    set(ax1,'TickLength',[0 0]);
    axis tight;
    xlabel('Time (s)',StyleArgs{:});
    ylabel('From Baseline',StyleArgs{:});
    
    if condIdx==1 || condIdx==2
        xl = xlim;
        xlim([xl(1) 6])
    elseif condIdx==3 || condIdx==4
        xl = xlim;
        xlim([-2.8 5])
    end
end

for plotIdx = 1:length(condNames)
    
    pnl(1,plotIdx).select()
    ylimtmp1(plotIdx,:) = ylim();
    
end

for plotIdx = 1:length(condNames)
    
    pnl(1,plotIdx).select()
    ylim([min([ylimtmp1(:,1);0]) max(ylimtmp1(:,2))]);
%      ylim([0 15]);
    tEvents = [0];
    plt.vline(tEvents,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEvents,ytext,'Go',StyleArgs{:});
    
    if plotIdx==3 || plotIdx==4
        
    tEventsCue = [-2];
    plt.vline(tEventsCue,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEventsCue,ytext,'Cue',StyleArgs{:});

    end
end

FigsDir=fullfile(BaseDir, 'PopData', 'Mine');
plt.SaveFigure(saveFig,FigsDir,['FeltImagFeltObs'],'SVG','PNG');
% plt.SaveFigure(saveFig,FigsDir,['FeltImagFeltObs'],'SVG');

end
