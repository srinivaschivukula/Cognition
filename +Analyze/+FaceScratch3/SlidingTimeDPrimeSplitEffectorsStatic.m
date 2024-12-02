% Address the basic question, when x is the salient category, what
% percentage of units are tuned to x.


function SlidingTimeDPrimeSplitEffectorsStatic(config,varargin)

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
    disp(sampIdx)
    %% first do all the felt/imagined stuff
    
    TimeWindows = [0 3];
    %     TimeWindows = [-4.5:TimeStep:5.5];
    %     TimeWindows = [-1.5:TimeStep:7.5];
    %     TimeCenters = TimeWindows+TimeStep;
%     TimeCenters=[];
%     for timeIdx = 1:length(TimeWindows)
%         cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
%         TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
%     end
    
    %% this portion is for felt vs imagined data
    AllImagUnits = vertcat(ImagData{1}.Unit{:});
    %         PVals = horzcat(ImagData{1}.Pvals{:})';
    %     H = Utilities.MultipleComparisonsCorrection(PVals,'method','fdr');
    %     AllImagUnits=AllImagUnits(any(H(:,[2,3,5,6]),2),:);
    randSampImag = randsample(length(AllImagUnits),100);
    
    cUnitsImag = AllImagUnits(randSampImag,:);
    
    condNamesImag = {'ImagCheek','ImagShoulder','CheekR','ShoulderR'};
    
    FeltTrialDataFI_cheek = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go', 'Condition', {'CheekR'});
    FeltTrialDataFI_shoulder = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go', 'Condition', {'ShoulderR'});
    
    ImaginedTrialData_cheek = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go', 'Condition', {'ImagCheek'});
    ImaginedTrialData_shoulder = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go', 'Condition', {'ImagShoulder'});
    
    %     BaselineDataFI_cheekFelt = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Delay','Condition',{'CheekR'});
    %     BaselineDataFI_shoulderFelt = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Delay','Condition',{'ShoulderR'});
    %     BaselineDataFI_cheekImag = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go','Condition',{'ImagCheek'});
    %     BaselineDataFI_shoulderImag = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go','Condition',{'ImagShoulder'});
%     BaselineDataFI_null = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go','Condition',{'XX'});
BaselineDataFI_null = Analyze.SubSelectTrials(ImagTrialData, 'Phase', 'Go','Condition',{'CheekR','ShoulderR'});
    
    %% do the sliding time calculation of distances/dprimes
    
    for timeIdx = 1:1
        
        cTW = [0 3];
        
        %% get the felt/imag data
        
        feltFRFI_cheek = [squeeze(Analyze.getNeuralData(FeltTrialDataFI_cheek, cUnitsImag, cTW))];
        feltFRFI_shoulder = [squeeze(Analyze.getNeuralData(FeltTrialDataFI_shoulder, cUnitsImag, cTW))];
        
        imagFRFI_cheek = [squeeze(Analyze.getNeuralData(ImaginedTrialData_cheek, cUnitsImag, cTW))];
        imagFRFI_shoulder = [squeeze(Analyze.getNeuralData(ImaginedTrialData_shoulder, cUnitsImag, cTW))];
        
%         if cTW(1)>0
%             baseFRFI = squeeze(Analyze.getNeuralData(BaselineDataFI_null, cUnitsImag, [0 3]));
%         else
            
            baseFRFI = squeeze(Analyze.getNeuralData(BaselineDataFI_null, cUnitsImag, [-1.5 -.5]));
%             baseFRFI = squeeze(Analyze.getNeuralData(BaselineDataFI_null, cUnitsImag, [0 3]));
%         end
        
        newFeltFRFI_cheek =[]; newFeltFRFI_shoulder =[]; newImagFRFI_cheek = []; newImagFRFI_shoulder = []; newBaseFRFI=[];
        for i=1:size(feltFRFI_cheek,2)
                tmp = feltFRFI_cheek(:,i);
                newFeltFRFI_cheek(:,i) = tmp(~isnan(tmp));
                
                                tmp = feltFRFI_shoulder(:,i);
                newFeltFRFI_shoulder(:,i) = tmp(~isnan(tmp));

            
            tmp = imagFRFI_cheek(:,i);
            newImagFRFI_cheek(:,i) = tmp(~isnan(tmp));
            
            tmp = imagFRFI_shoulder(:,i);
            newImagFRFI_shoulder(:,i) = tmp(~isnan(tmp));
            
            tmp = baseFRFI(:,i);
            newBaseFRFI(:,i) = tmp(~isnan(tmp));
            
        end
        
        feltFRFI_cheek = newFeltFRFI_cheek; feltFRFI_shoulder = newFeltFRFI_shoulder; imagFRFI_cheek=newImagFRFI_cheek; imagFRFI_shoulder=newImagFRFI_shoulder; baseFRFI=newBaseFRFI;
        newFeltFRFI_cheek=[]; newObsFRFI_cheek=[]; newFeltFRFI_shoulder=[]; newObsFRFI_shoulder=[]; newBaseFRFI=[];
        
        %% calculate all the distances/dprimes
        
        ResultsFeltFI_cheek = Analyze.FitPLSLM([feltFRFI_cheek;baseFRFI],[ones(size(feltFRFI_cheek,1),1);zeros(size(baseFRFI,1),1)],'CrossValidate',1);
        ResultsFeltFI_shoulder = Analyze.FitPLSLM([feltFRFI_shoulder;baseFRFI],[ones(size(feltFRFI_shoulder,1),1);zeros(size(baseFRFI,1),1)],'CrossValidate',1);
        ResultsObsFI_cheek = Analyze.FitPLSLM([imagFRFI_cheek;baseFRFI],[ones(size(imagFRFI_cheek,1),1);zeros(size(baseFRFI,1),1)],'CrossValidate',1);
        ResultsObsFI_shoulder = Analyze.FitPLSLM([imagFRFI_shoulder;baseFRFI],[ones(size(imagFRFI_shoulder,1),1);zeros(size(baseFRFI,1),1)],'CrossValidate',1);
        
        DataFeltFI_cheek{sampIdx}{timeIdx} = ResultsFeltFI_cheek.Dist;
        DataFeltFI_shoulder{sampIdx}{timeIdx} = ResultsFeltFI_shoulder.Dist;
        DataObsFI_cheek{sampIdx}{timeIdx} = ResultsObsFI_cheek.Dist;
        DataObsFI_shoulder{sampIdx}{timeIdx} = ResultsObsFI_shoulder.Dist;
                
    
    end
    
    %% this section is for felt vs observed data
    
    


%% post processing and plotting

%% plot felt vs imagined    TimeWindows = [-2.5:TimeStep:5];
    %     TimeCenters = TimeWindows+TimeStep;
%     TimeCenters=[];
%     for timeIdx = 1:length(TimeWindows)
%         cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
%         TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
%     end
    
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
    
    
    for timeIdx = 1:1
        %     for timeIdx=1
        
        cTW = [0 4];
        
        %% get the felt/obs data
        
        feltFR_cheek = [squeeze(Analyze.getNeuralData(FeltTrialData_cheek, cUnits, cTW))];
        feltFR_shoulder = [squeeze(Analyze.getNeuralData(FeltTrialData_shoulder, cUnits, cTW))];
        
        obsFR_cheek = [squeeze(Analyze.getNeuralData(ObsTrialData_cheek, cUnits, cTW))];
        obsFR_shoulder = [squeeze(Analyze.getNeuralData(ObsTrialData_shoulder, cUnits, cTW))];
        
%         if cTW(1)>0
%             baseFR = [squeeze(Analyze.getNeuralData(BaselineData_cheekFelt, cUnits, [-1 -0]));squeeze(Analyze.getNeuralData(BaselineData_shoulderFelt, cUnits, [-1 -0]));squeeze(Analyze.getNeuralData(BaselineData_cheekObs, cUnits, [-1 -0]));squeeze(Analyze.getNeuralData(BaselineData_shoulderObs, cUnits, [-1 -0]))];
%         else
            baseFR = [squeeze(Analyze.getNeuralData(BaselineData_cheekFelt, cUnits, [-3.5 -2.5]));squeeze(Analyze.getNeuralData(BaselineData_shoulderFelt, cUnits, [-3.5 -2.5]));...
                squeeze(Analyze.getNeuralData(BaselineData_cheekObs, cUnits, [-3.5 -2.5]));squeeze(Analyze.getNeuralData(BaselineData_shoulderObs, cUnits, [-3.5 -2.5]))];
%         end
        
        newFeltFR_cheek =[]; newFeltFR_shoulder =[]; newObsFR_cheek = []; newObsFR_shoulder = []; newBaseFR=[];
        for i=1:size(feltFR_cheek,2)
            tmp = feltFR_cheek(:,i);
            newFeltFR_cheek(:,i) = tmp(~isnan(tmp));
            
            tmp = feltFR_shoulder(:,i);
            newFeltFR_shoulder(:,i) = tmp(~isnan(tmp));
            
            tmp = obsFR_cheek(:,i);
            newObsFR_cheek(:,i) = tmp(~isnan(tmp));
            
            tmp = obsFR_shoulder(:,i);
            newObsFR_shoulder(:,i) = tmp(~isnan(tmp));
            
            tmp = baseFR(:,i);
            newBaseFR(:,i) = tmp(~isnan(tmp));
        end
        feltFR_cheek = newFeltFR_cheek; feltFR_shoulder = newFeltFR_shoulder; obsFR_cheek=newObsFR_cheek; obsFR_shoulder=newObsFR_shoulder; baseFR=newBaseFR;
        newFeltFR_cheek=[]; newFeltFR_shoulder=[]; newObsFR_cheek=[]; newObsFR_shoulder=[]; newBaseFR=[];
        
        
        ResultsFelt_cheek = Analyze.FitPLSLM([feltFR_cheek;baseFR],[ones(size(feltFR_cheek,1),1);zeros(size(baseFR,1),1)],'CrossValidate',1);
        ResultsFelt_shoulder = Analyze.FitPLSLM([feltFR_shoulder;baseFR],[ones(size(feltFR_shoulder,1),1);zeros(size(baseFR,1),1)],'CrossValidate',1);
        ResultsObs_cheek = Analyze.FitPLSLM([obsFR_cheek;baseFR],[ones(size(obsFR_cheek,1),1);zeros(size(baseFR,1),1)],'CrossValidate',1);
        ResultsObs_shoulder = Analyze.FitPLSLM([obsFR_shoulder;baseFR],[ones(size(obsFR_shoulder,1),1);zeros(size(baseFR,1),1)],'CrossValidate',1);
        
        DataFelt_cheek{sampIdx}{timeIdx} = ResultsFelt_cheek.Dist;
        DataFelt_shoulder{sampIdx}{timeIdx} = ResultsFelt_shoulder.Dist;
        DataObs_cheek{sampIdx}{timeIdx} = ResultsObs_cheek.Dist;
        DataObs_shoulder{sampIdx}{timeIdx} = ResultsObs_shoulder.Dist;
        
    end

end

FeltFI_cheek=[]; FeltFI_shoulder=[]; Imag_cheek=[]; Imag_shoulder=[];
for i=1:length(DataFeltFI_cheek)
    
    tmp = [DataFeltFI_cheek{i}{:}];
    FeltFI_cheek = cat(1,FeltFI_cheek,tmp);
    
    tmp = [DataFeltFI_shoulder{i}{:}];
    FeltFI_shoulder = cat(1,FeltFI_shoulder,tmp);
    
    tmp = [DataObsFI_cheek{i}{:}];
    Imag_cheek = cat(1,Imag_cheek,tmp);
    
    tmp = [DataObsFI_shoulder{i}{:}];
    Imag_shoulder = cat(1,Imag_shoulder,tmp);
end
PlotData = {{FeltFI_cheek},{Imag_cheek},{FeltFI_shoulder},{Imag_shoulder}};

StyleArgs = {'FontName','Arial','FontSize',9,'FontWeight','bold'};

figH = plt.fig('units','inches','width',14,'height',8,'font','Arial','fontsize',9);
figure(figH);
clf; pnl = panel(); pnl.margin=15; pnl.pack(3,7);

condNames = {'Felt Cheek','Imag Cheek','Felt Shoulder','Imag Shoulder'};

for condIdx=1:4
    
    cData = PlotData{condIdx};
    
    clr=[];
    clr(1,:) = repmat([237 0 38]./255,1,1);
    clr(2,:) = repmat([31 120 180]./255,1,1);
    clr(3,:) = repmat([237 0 38]./255,1,1);
    clr(4,:) = repmat([1 108 89]./255,1,1);
    
    cLabel = condNames{condIdx};
    
    TimeWindows = [-2.5:TimeStep:4];
    TimeCenters=[];
    for timeIdx = 1:length(TimeWindows)
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
    end
    
    % plot the discriminability of touch from baseline
    pnl(1,condIdx).select()
    
    Analyze.plotEventRelatedAverage(cData,cLabel,'TimeVec',TimeCenters,'usePercentile',...
        'Colors',clr(condIdx,:));
    title(cLabel);
    ax1 = gca;
    set(ax1,'box','off',StyleArgs{:});
    ax1.XRuler.Axle.LineStyle = 'none';
    set(gca,'TickDir','out');
    set(ax1,'TickLength',[0 0]);
    axis tight;
    xlabel('Time (s)',StyleArgs{:});
    ylabel('From Baseline',StyleArgs{:});
    
    xl = xlim;
    xlim([-2.5 4])
    
end

for plotIdx = 1:4
    
    pnl(1,plotIdx).select()
    ylimtmp1(plotIdx,:) = ylim();
    
end

for plotIdx = 1:4
    
    pnl(1,plotIdx).select()
    ylim([min([ylimtmp1(:,1);0]) max(ylimtmp1(:,2))]);
    %      ylim([0 10]);
    % xlim([0.2 3.5]);
    tEvents = [0];
    plt.vline(tEvents,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEvents,ytext-2,'Go',StyleArgs{:});
    
    
    tEventsCue = [-2];
    plt.vline(tEventsCue,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEventsCue,ytext-2,'Cue',StyleArgs{:});
    
end

FigsDir=fullfile(BaseDir, 'PopData', 'Mine');
plt.SaveFigure(saveFig,FigsDir,['FeltImagByBodyPart'],'SVG','PNG');
% plt.SaveFigure(saveFig,FigsDir,['FeltImagByBodyPart'],'SVG');


%% plot felt vs observed

Felt_cheek=[]; Felt_shoulder=[]; Obs_cheek=[]; Obs_shoulder=[];
for i=1:length(DataFelt_cheek)
    
    tmp = [DataFelt_cheek{i}{:}];
    Felt_cheek = cat(1,Felt_cheek,tmp);
    
    tmp = [DataFelt_shoulder{i}{:}];
    Felt_shoulder = cat(1,Felt_shoulder,tmp);
    
    tmp = [DataObs_cheek{i}{:}];
    Obs_cheek = cat(1,Obs_cheek,tmp);
    
    tmp = [DataObs_shoulder{i}{:}];
    Obs_shoulder = cat(1,Obs_shoulder,tmp);
end
PlotData = {{Felt_cheek},{Obs_cheek},{Felt_shoulder},{Obs_shoulder}};

StyleArgs = {'FontName','Arial','FontSize',9,'FontWeight','bold'};

figH = plt.fig('units','inches','width',14,'height',8,'font','Arial','fontsize',9);
figure(figH);
clf; pnl = panel(); pnl.margin=15; pnl.pack(3,7);

condNames = {'Felt Cheek','Obs Cheek','Felt Shoulder','Obs Shoulder'};

for condIdx=1:4
    
    cData = PlotData{condIdx};
    
    clr=[];
    clr(1,:) = repmat([237 0 38]./255,1,1);
    clr(2,:) = repmat([1 108 89]./255,1,1);
    clr(3,:) = repmat([237 0 38]./255,1,1);
    clr(4,:) = repmat([1 108 89]./255,1,1);
    
    cLabel = condNames{condIdx};
    
    TimeWindows = [-2.5:TimeStep:5];
    TimeCenters=[];
    for timeIdx = 1:length(TimeWindows)
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
    end
    
    % plot the discriminability of touch from baseline
    pnl(1,condIdx).select()
    
    Analyze.plotEventRelatedAverage(cData,cLabel,'TimeVec',TimeCenters,'usePercentile',...
        'Colors',clr(condIdx,:));
    title(cLabel);
    ax1 = gca;
    set(ax1,'box','off',StyleArgs{:});
    ax1.XRuler.Axle.LineStyle = 'none';
    set(gca,'TickDir','out');
    set(ax1,'TickLength',[0 0]);
    axis tight;
    xlabel('Time (s)',StyleArgs{:});
    ylabel('From Baseline',StyleArgs{:});
    
    xl = xlim;
    xlim([-2.5 5.5])
    
end

for plotIdx = 1:4
    
    pnl(1,plotIdx).select()
    ylimtmp1(plotIdx,:) = ylim();
    
end

for plotIdx = 1:4
    
    pnl(1,plotIdx).select()
    ylim([min([ylimtmp1(:,1);0]) max(ylimtmp1(:,2))]);
    %      ylim([0 10]);
    % xlim([0.2 3.5]);
    tEvents = [0];
    plt.vline(tEvents,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEvents,ytext-2,'Go',StyleArgs{:});
    
    
    tEventsCue = [-2];
    plt.vline(tEventsCue,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEventsCue,ytext-2,'Cue',StyleArgs{:});
    
end

FigsDir=fullfile(BaseDir, 'PopData', 'Mine');
plt.SaveFigure(saveFig,FigsDir,['FeltObsByBodyPart'],'SVG','PNG');
% plt.SaveFigure(saveFig,FigsDir,['FeltObsByBodyPart'],'SVG');


%% plot all felt, observed, imagined, averaged across body parts

Felt=[]; Obs=[]; FeltFI=[]; Imag=[];
Felt = [Felt_cheek+Felt_shoulder]./2; Obs = [Obs_cheek+Obs_shoulder]./2;
FeltFI = [FeltFI_cheek+FeltFI_shoulder]./2; Imag = [Imag_cheek+Imag_shoulder]./2;
PlotData = {{Felt},{Obs},{FeltFI},{Imag}};

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
        
        TimeWindows = [-2.5:TimeStep:5];
        TimeCenters=[];
        for timeIdx = 1:length(TimeWindows)
            cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
            TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
        end
        
    elseif condIdx==3 || condIdx==4
        TimeWindows = [-2.5:TimeStep:4];
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
    
    title(cLabel)
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
    text(tEvents,ytext-2,'Go',StyleArgs{:});
    
    
    tEventsCue = [-2];
    plt.vline(tEventsCue,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEventsCue,ytext-2,'Cue',StyleArgs{:});
    
end


FigsDir=fullfile(BaseDir, 'PopData', 'Mine');
plt.SaveFigure(saveFig,FigsDir,['FeltImagFeltObs'],'SVG','PNG');
% plt.SaveFigure(saveFig,FigsDir,['FeltImagFeltObs'],'SVG');

end
