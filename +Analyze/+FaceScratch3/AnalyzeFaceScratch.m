function AnalyzeFaceScratch()
Convert();

% Possible args: {'SUA','Pop','ERA','Decode'};

BodySpec();
BodySpecC();
ActionType();
end

function Convert()
%% Body spec

fn ={};
tl = {};
filename = {'20170222-120530-120650-TextTraining',...
    '20170222-120530-121309-TextTraining',...
    '20170315-104616-104804-TextTraining',...
    '20170315-104616-105552-TextTraining',...
    '20170322-105656-105955-TextTraining',...
    '20170405-102741-124026-TextTraining',...
    '20170407-102824-103455-TextTraining',...
    '20170426-110048-110746-TextTraining'};

taskLabel = {'XBodySpec','XBodySpec','XBodySpec','XBodySpec',...
    'XBodySpec','XBodySpec','XBodySpec','XBodySpec'};

fn = [fn filename];
tl = [tl taskLabel];
%
filename = fn;
taskLabel = tl;

for i = 1:length(filename)
    Analyze.ProcessRawData('FaceScratch3',...
        'TaskFile',filename{i},...
        'Run',i,...
        'overWrite',1,...
        'TaskLabel',taskLabel{i},...
        'TaskConfig',@Analyze.FaceScratch3.Config_text);
end

%% ActionType

fn ={};
tl = {};
filename = {'20171215-104431-104555-TextTraining',...
    '20171215-104431-105543-TextTraining',...
    '20171220-111903-112826-TextTraining',...
    '20171220-111903-114011-TextTraining',...
    '20180105-121305-121827-TextTraining',...
    '20180105-121305-122928-TextTraining',...
    '20180110-112402-113112-TextTraining',...
    '20180110-112402-114157-TextTraining',...
    '20180112-114135-115553-TextTraining',...
    '20180112-114135-120615-TextTraining',...
    '20180115-110243-110540-TextTraining',...
    '20180115-110243-111723-TextTraining',...
    '20180117-113257-113554-TextTraining',...
    '20180117-113257-115125-TextTraining'};

taskLabel = {'XActionType','XActionType','XActionType','XActionType',...
    'XActionType','XActionType','XActionType','XActionType',...
    'XActionType','XActionType','XActionType','XActionType',...
    'XActionType','XActionType'};

fn = [fn filename];
tl = [tl taskLabel];
%
filename = fn;
taskLabel = tl;

for i = 1:length(filename)
    Analyze.ProcessRawData('FaceScratch3',...
        'TaskFile',filename{i},...
        'Run',i,...
        'overWrite',1,...
        'TaskLabel',taskLabel{i},...
        'TaskConfig',@Analyze.FaceScratch3.Config_text);
end

end

function BodySpec(varargin)
Dates = {'20170222','20170315','20170322','20170405','20170407','20170426'};
% Dates = {'20170407'};
Phase = 'Go';
Labels = {'NancyCheek','NancyHead','NancyShoulder','TysonCheek','TysonHead','TysonShoulder'};
GlobalTrialCriteria = {'TaskLabel','XBodySpec'};

ModelGroups = {{{{'NancyCheek','NancyShoulder'},{'TysonCheek','TysonShoulder'}},...
    {{'NancyCheek','TysonCheek'},{'NancyShoulder','TysonShoulder'}}}};
ModelNames = {{'Person','BodyPart'}};

ModelCompare = {...
    {'NancyCheek',      1,1,1,0,1,1,0,1;...
    'NancyShoulder',    1,2,2,0,2,0,1,1;...
    'TysonCheek',       1,3,0,1,1,2,0,2;...
    'TysonShoulder',    1,4,0,2,2,0,2,2}};
ModelCompareNames = {{'Invariant','Idiosyncratic',...
    'BPNancyOnly','BPTysonOnly','BPSpec',...
    'PCheekOnly','PShoulderOnly','PersonSpec'}};

Tag = GlobalTrialCriteria{2};
OutDir = fullfile(env.get('result'),'FaceScratch3','SUAnal',[Tag '-Go']);

Phases = {'CueTarget','Delay','Go'};
TimeWindows = [-.5 2.5; 0 3; 0 4.5];

PhaseStart=1;
PhaseDur =3;
BaselineWindow = [1.5 2.5];

%%
AllTrialData = Analyze.LoadConvertedData('FaceScratch3',Dates);
AllTrialData = Analyze.SubSelectTrials(AllTrialData,GlobalTrialCriteria{:});

%% demixed PCA
% 
% timeWindow=[-6:0.05:4.5];
% timeEvents=[-5.5 -3 0]; % cue onset, delay onset, go onset
% 
% SmoothingKernel=.85;
% Phase='Go';
% 
% [~,FR]=Analyze.FaceScratch3.getDataForDPCA(AllTrialData,timeWindow,...
%     Phase,'meanThresh',1.5,'SmoothingKernel',SmoothingKernel,GlobalTrialCriteria{:});
% 
% [FRcond,FRall,trialNum]=Cell2dPCA2Mod(FR);
% 
% Out=Analyze.FaceScratch3.dpcaAnalysis(FRcond,FRall,trialNum,timeWindow(1:end-1),'timeEvents',timeEvents,'numComponents',12,...
%     'decode_Reg','plotPCABasic',1,'plotPCAMarg',0);

%% for the sliding time comparison of felt/obs to felt/imag
% need to load trialdata also for imag/felt task
ImagTrialData = Analyze.LoadConvertedData('FaceScratch',{'20180709','20180711',...
    '20180716','20180718','20180723','20180730','20180806','20180815'});
ImagTrialData = Analyze.SubSelectTrials(ImagTrialData,'TaskLabel','XImagery2');

%% Single unit analysis: BodySpec
tic
AU.SingleUnitAnalysis(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,'GlobalTrialCriteria',GlobalTrialCriteria,...
    'ModelGroups',ModelGroups,'ModelNames',ModelNames,...
    'ModelCompare',ModelCompare,'ModelCompareNames',ModelCompareNames);
toc
%% PlotAnalysisThroughTime: BodySpec
AU.PlotAnalysisThroughTime(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'ModelGroups',ModelGroups,'ModelNames',ModelNames,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'ModelCompare',ModelCompare,'ModelCompareNames',ModelCompareNames);

%% Tuning and Discriminability

models2use = logical([1 1 1 1 1 0 0 1]);

Analyze.FaceScratch3.TuningDPrime(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Pmethod','fdr','saveFig',...
    'AllTrialData',AllTrialData,...
    'ModelCompare',ModelCompare,'ModelCompareNames',ModelCompareNames,...
    'models2use',models2use,...
    'SelectionMode','sigUnits');

%% Sliding TimeWindow Discriminability Index for Imag/Felt and for Obs/Felt

% step size is how big of a window to calculate the metric on
% time step is how spaced out to make the new time windows
% so a StepSize<TimeStep is not overlapping. StepSize=TimeStep is continous
% non-overlapping. StepSize>TimeStep is overlapping windows.
TimeStep = 0.1;
StepSize = 0.5;
NumSamples = 100;
% NumSamples = 75;

%subtractBaseline is if we want to eliminate the non-zero pre-stimulus
%activity
SubtractBaseline=0;

% Analyze.FaceScratch3.SlidingTimeDPrime(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
%     'Dates',Dates,'Phase',Phase,'Labels',Labels,...
%     'GlobalTrialCriteria',GlobalTrialCriteria,...
%     'Pmethod','fdr','saveFig',...
%     'AllTrialData', AllTrialData,...
%     'ImagTrialData', ImagTrialData,...
%     'TimeStep',TimeStep,...
%     'StepSize',StepSize,...
%     'NumSamples',NumSamples);
% Analyze.FaceScratch3.SlidingTimeDPrimeSplitEffectors(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
%     'Dates',Dates,'Phase',Phase,'Labels',Labels,...
%     'GlobalTrialCriteria',GlobalTrialCriteria,...
%     'Pmethod','fdr','saveFig',...
%     'AllTrialData', AllTrialData,...
%     'ImagTrialData', ImagTrialData,...
%     'TimeStep',TimeStep,...
%     'StepSize',StepSize,...
%     'NumSamples',NumSamples);
Analyze.FaceScratch3.SlidingTimeDPrimeSplitEffectorsStatic(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Pmethod','fdr','saveFig',...
    'AllTrialData', AllTrialData,...
    'ImagTrialData', ImagTrialData,...
    'TimeStep',TimeStep,...
    'StepSize',StepSize,...
    'NumSamples',NumSamples);

%% Plot SingleUnitERAs for specific interesting units
% [Units,~] = Analyze.FaceScratch3.FindInterestingUnits(GlobalTrialCriteria{2});

Units = [1 8 1 1513;...
    1 10 1 1513;...
    1 22 1 1513;...
    1 23 3 1513;...
    1 29 1 1513;...
    1 32 1 1513;...
    1 60 3 1513;...
    1 64 2 1513;...
    1 94 3 1513;...
    1 96 2 1513;...
    
    1 9 1 1534;...
    1 10 3 1534;...
    1 64 1 1534;...
    1 81 2 1534;...
    1 61 1 1541;...
    1 64 1 1541;...
    1 30 2 1555;...
    1 30 3 1555;...
    1 32 1 1555;...
    ];

clr=[];
clr(1,:) = repmat([237 0 38]./255,1,1);
clr(2,:) = repmat([31 120 180]./255,1,1);

LabelsAbbrev = {'FeelCheek','FeelShoulder','ObsCheek','ObsShoulder'};

PlotConds = {{'Condition','NancyCheek'},...
    {'Condition','NancyShoulder'},...
    {'Condition','TysonCheek'},...
    {'Condition','TysonShoulder'}};
SubplotIds = [1 1 2 2];

for i=1:length(Units)
    
    Analyze.FaceScratch3.PlotSingleUnitERAs(AllTrialData,PlotConds,Phases,TimeWindows,OutDir,...
        'SubplotIds',SubplotIds,'Units',Units(i,:),'FigSize',[4 5],...
        'LabelsAbbrev','',...
        'Colormap',clr,...
        'TaskLabel', Tag,'saveFig');
end

%% Specificity analysis
Labels2Incl = {'NancyCheek','NancyShoulder','TysonCheek','TysonShoulder'};
LabelsAbbrev = {'Felt Cheek','Felt Shoulder','Observed Cheek','Observed Shoulder'};

Analyze.FaceScratch3.SpecificityAnalysis(Tag,'Labels',Labels,'Labels2Incl',Labels2Incl,...
    'LabelsAbbrev',LabelsAbbrev,'saveFig','AllTrialData', AllTrialData);
% FaceScratch.SpecificityAnalysis(Tag,'Labels',Labels,'Labels2Incl',Labels2Incl,...
%     'LabelsAbbrev',LabelsAbbrev,'saveFig');

%% R2 by ANOVA type
% FaceScratch.ANOVATuningProperties(Tag,'Labels',Labels,'Labels2Incl','Labels2


%% Plot PCAs
Phases = {'CueTarget','Delay','Go'};
TimeWindows = [-.5 2.5; 0 3; 0 4.5];
PlotConds = {{'Condition','NancyCheek'},...
    {'Condition','TysonCheek'},...
    {'Condition','NancyShoulder'},...
    {'Condition','TysonShoulder'}};
Analyze.FaceScratch3.PlotPCAs(AllTrialData,PlotConds,Phases,TimeWindows,OutDir,'PlotERA','NumPC',4);


%% Plot PCAs (normalized)
Phases = {'CueTarget','Delay','Go'};
TimeWindows = [-.5 2.5; 0 3; 0 4.5];
PlotConds = {{'Condition','NancyCheek'},...
    {'Condition','TysonCheek'},...
    {'Condition','NancyHead'},...
    {'Condition','TysonHead'},...
    {'Condition','NancyShoulder'},...
    {'Condition','TysonShoulder'}};
Analyze.FaceScratch3.PlotPCAs(AllTrialData,PlotConds,Phases,TimeWindows,OutDir,'PlotERA','NumPC',4,'Normalize',1);

%% Plot PCA as 2D movies
Phases = {'CueTarget','Delay','Go'};
TimeWindows = [0 2.5; 0 3; 0 4.5];
Frames2Save = {{'CueTarget','Delay','Go'},[25,80,155]};
% SaveFramesIdx = [{'CueTarget',25},{'Delay',80}{'Go',155}];
Analyze.FaceScratch3.PlotPCAs(AllTrialData,PlotConds,Phases,TimeWindows,OutDir,'Plot2D','PCPairsMaxIdx',6,...
    'Frames2Save',Frames2Save);

%% Plot correlation matrix

DecodeLabels = {'NancyCheek','NancyShoulder','TysonCheek','TysonShoulder'};
LabelsAbbrev = {{'Felt', 'Observed'},{'Cheek','Shoulder'}};

Analyze.FaceScratch3.CorrMatrix(Labels,DecodeLabels,OutDir,Tag,'reorder',[1 3 2 4],...
    'linkagetype','weighted',...
    'LabelsAbbrev',LabelsAbbrev,...
    'Pmethod','fdr',...
    'AllTrialData', AllTrialData,'saveFig',...
    'PhaseStart', PhaseStart, 'PhaseDur', PhaseDur, 'BaselineWindow', BaselineWindow);

%% Decode Analysis (Confusion Matrix)

TimeWindow = [1 4];
DecodeLabels = {'NancyCheek','NancyShoulder','TysonCheek','TysonShoulder'};
LabelsAbbrev = {{'Felt', 'Observed'},{'Cheek','Shoulder'}};

Analyze.FaceScratch3.DecodeAnalysis(AllTrialData,Labels,DecodeLabels,TimeWindow,...
    OutDir,Tag,'reorder',[1 3 2 4],'LabelsAbbrev',LabelsAbbrev,'saveFig');

%% Sliding Time Window for the Confusion Matrix

TimeWindow = [-1:7.5];
TimeStep = 0.15;
StepSize = 0.3;

DecodeLabels = {'NancyCheek','NancyShoulder','TysonCheek','TysonShoulder'};
LabelsAbbrev = {{'Felt', 'Observed'},{'Cheek','Shoulder'}};

Analyze.FaceScratch3.SlidingTimeDecode(AllTrialData,...
    Labels,DecodeLabels,TimeWindow,TimeStep,StepSize,OutDir,Tag,...
    'reorder',[1 3 2 4],'LabelsAbbrev',LabelsAbbrev,'saveFig');


%% Cross decode analysis
TimeWindow = [1 4];
DecodeLabels = {'NancyCheek','NancyShoulder','TysonCheek','TysonShoulder'};
LabelsAbbrev = {{'Felt', 'Observed'},{'Cheek','Shoulder'}};

% Analyze.FaceScratch3.CrossDecodeAnalysis(AllTrialData,Labels,...
%     DecodeLabels,TimeWindow,OutDir,Tag,...
%     'LabelsAbbrev',LabelsAbbrev,...
%     'saveFig');
Analyze.FaceScratch3.CrossDecodeAnalysisImproved(AllTrialData,Labels,...
    DecodeLabels,TimeWindow,OutDir,Tag,...
    'LabelsAbbrev',LabelsAbbrev,...
    'saveFig');

%% Cross decode analysis - New Version

TimeWindow = [1 4];
DecodeLabels = {'NancyCheek','TysonCheek','NancyShoulder','TysonShoulder'};

LabelsAbbrev = {{'Felt', 'Observed'},{'Cheek','Shoulder'}};

SplitDims = {...
    
    'NancyCheek',        1,1;...
    'TysonCheek',        2,1;...
    'NancyShoulder',     1,2;...
    'TysonShoulder',     2,2;...
    
    };

SplitValNames = {{'Nancy','Tyson'},{'Cheek','Shoulder'}}; % Label names, correspondign to splits (1,2, etc)

DimNames = {'Person','BodyPart'};

Analyze.FaceScratch3.CrossDecodeAnalysisTrial(AllTrialData,Labels,DecodeLabels,TimeWindow,OutDir,Tag,...);
    'SplitDims',SplitDims,'DimNames',DimNames,'SplitValNames',SplitValNames,...
    'TaskLabel', Tag,'LabelsAbbrev',LabelsAbbrev,...
    'All','saveFig')

%% Dropping curve analysis
TimeWindow = [1 4];
DecodeLabels = {'NancyCheek','NancyShoulder','TysonCheek','TysonShoulder'};
Analyze.FaceScratch.DroppingCurve(AllTrialData,Labels,DecodeLabels,TimeWindow,OutDir,Tag);

%% Single unit  - Stimulus Profiles

AU.SingleUnitAnalysis(@Analyze.FaceScratch3.Config_LinearAnalysisText,'Tag','StimProfiles',...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,'GlobalTrialCriteria',GlobalTrialCriteria);
%% Stimulus Profiles PCA

Analyze.RFHandEye2.SUP_StimulusProfiles(@Analyze.FaceScratch3.Config_LinearAnalysisText,'Tag', 'StimProfiles',...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Pmethod','fdr',...
    'saveFig',...
    'AllTrialData', AllTrialData);

%% Single unit  - ModelComparison

ModelCompare = {...
    {'NancyCheek',      1,1,1,0,1,1;...
    'NancyShoulder',    1,2,2,0,2,1;...
    'TysonCheek',       1,3,0,1,1,2;...
    'TysonShoulder',    1,4,0,2,2,2}};
ModelCompareNames = {{'Invariant','Idiosyncratic',...
    'BPNancyOnly','BPTysonOnly','BPSpec',...
    'PersonSpec'}};

AU.SingleUnitAnalysis(@Analyze.FaceScratch3.Config_LinearAnalysisText,'Tag','ModelComparison',...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,'GlobalTrialCriteria',GlobalTrialCriteria,...
    'ModelCompare',ModelCompare,'ModelCompareNames',ModelCompareNames);
%% Stimulus Profiles PCA

ModelCompare = {...
    {'NancyCheek',      1,1,1,0,1,1;...
    'NancyShoulder',    1,2,2,0,2,1;...
    'TysonCheek',       1,3,0,1,1,2;...
    'TysonShoulder',    1,4,0,2,2,2}};
ModelCompareNames = {{'Invariant','Idiosyncratic',...
    'BPNancyOnly','BPTysonOnly','BPSpec',...
    'PersonSpec'}};

models2use = logical([1 1 1 1 1 1]);

Analyze.FaceScratch3.SUP_ModelComparison(@Analyze.FaceScratch3.Config_LinearAnalysisText,'Tag', 'ModelComparison',...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'ModelCompare',ModelCompare,'ModelCompareNames',ModelCompareNames,...
    'Pmethod','fdr',...
    'saveFig',...
    'AllTrialData', AllTrialData,...
    'models2use',models2use);


end

function ActionType(varargin)
Dates = {'20171215','20171220','20180105','20180110','20180112','20180115','20180117'};

% 20171117 was gabe touching nancy, carey touching sri
GlobalTrialCriteria = {'TaskLabel','XActionType'};

Phase = 'Go';
Labels = {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek',...
    'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder',...
    'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek',...
    'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder',...
    'Null'};

ModelGroups = {...
    {{{'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek','NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'},...
    {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek','TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'}},...
    {{'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek','TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'},...
    {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder','TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'}},...
    {{'NancyPinchCheek','NancyPinchShoulder','TysonPinchCheek','TysonPinchShoulder'},...
    {'NancyPressCheek','NancyPressShoulder','TysonPressCheek','TysonPressShoulder'},...
    {'NancyRubCheek','NancyRubShoulder','TysonRubCheek','TysonRubShoulder'},...
    {'NancyTapCheek','NancyTapShoulder','TysonTapCheek','TysonTapShoulder'}}},...
    };
ModelNames = {{'Person','BodyPart','Action'}};

Tag = GlobalTrialCriteria{2};
OutDir = fullfile(env.get('result'),'FaceScratch3','SUAnal',[Tag '-Go']);

PhaseStart = 0.5;
PhaseDur = 2.5;
BaselinePhase = {'Phase','Go','Condition','Null'};
BaselineWindow = [.5 2.5];
%%
AllTrialData = Analyze.LoadConvertedData('FaceScratch3',Dates);
AllTrialData = Analyze.SubSelectTrials(AllTrialData,GlobalTrialCriteria{:});

%% demixed PCA

timeWindow=[-6:0.05:4.5];
timeEvents=[-5.5 -3 0]; % cue onset, delay onset, go onset

SmoothingKernel=.85;
Phase='Go';

[~,FR]=Analyze.FaceScratch3.getDataForDPCAAction(AllTrialData,timeWindow,...
    Phase,'meanThresh',1.5,'SmoothingKernel',SmoothingKernel,GlobalTrialCriteria{:});

[FRcond,FRall,trialNum]=Cell2dPCA2Mod(FR);

Out=Analyze.FaceScratch3.dpcaAnalysisAction(FRcond,FRall,trialNum,timeWindow(1:end-1),'timeEvents',timeEvents,'numComponents',12,...
    'decode_Reg','plotPCABasic',1,'plotPCAMarg',0);
%% Single unit analysis: Action Sensation

tic
AU.SingleUnitAnalysis(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,'GlobalTrialCriteria',GlobalTrialCriteria,...
    'ModelGroups',ModelGroups,'ModelNames',ModelNames,...
    'BaselinePhase',{'Phase','Go','Condition','Null'},'BaselineWindow',[.5 2.5],'PhaseStart',0.5,'PhaseDur',2.5);
toc
%% PlotAnalysisThroughTime: BodySpec
AU.PlotAnalysisThroughTime(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'ModelGroups',ModelGroups,'ModelNames',ModelNames,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'BaselinePhase',{'Phase','Go','Condition','Null'},'BaselineWindow',[.5 2.5],'PhaseStart',0.5,'PhaseDur',2.5);

%% Single unit  - Stimulus Profiles

AU.SingleUnitAnalysis(@Analyze.FaceScratch3.Config_LinearAnalysisText,'Tag','StimProfiles',...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,'GlobalTrialCriteria',GlobalTrialCriteria);
%% Stimulus Profiles PCA

Analyze.FaceScratch3.SUP_StimulusProfilesAction(@Analyze.FaceScratch3.Config_LinearAnalysisText,'Tag', 'StimProfiles',...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Pmethod','fdr',...
    'saveFig',...
    'AllTrialData', AllTrialData);

%% Stimulus Profiles Action Modified 

Analyze.FaceScratch3.SUP_StimulusProfilesActionMod(@Analyze.FaceScratch3.Config_LinearAnalysisText,'Tag', 'StimProfiles',...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Pmethod','fdr',...
    'saveFig',...
    'AllTrialData', AllTrialData);
%% Tuning and Discriminability

% models2use = logical([1 1 1 1 1 0 0 1]);

Analyze.FaceScratch3.TuningDPrimeAction(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Pmethod','fdr','saveFig',...
    'AllTrialData',AllTrialData);

%% Tuning and Discriminability - using anova rather than regression


Analyze.FaceScratch3.TuningDPrimeActionVersion2(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Pmethod','fdr','saveFig',...
    'AllTrialData',AllTrialData);

%% Correlation Matrix

Labels2Incl = {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek',...
    'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek',...
    'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder',...
    'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'};
LabelsAbbrev = {'NPiC','NPrC','NRbC','NTpC',...
    'TPiC','TPrC','TRbC','TTpC',...
    'NPiS','NPrS','NRbS','NTpS',...
    'TPiS','TPrS','TRbS','TTpS'
    };

reorder = [1 2 3 4 9 10 11 12 5 6 7 8 13 14 15 16];

LabelsNames = {{'Felt Pinch','Obs Pinch','Felt Press','Obs Press',...
    'Felt Rub','Obs Rub','Felt Tap','Obs Tap',...
    'Felt Pinch','Obs Pinch','Felt Press','Obs Press',...
    'Felt Rub','Obs Rub','Felt Tap','Obs Tap'},{'Cheek', 'Shoulder'}};

Analyze.FaceScratch3.CorrMatrixAction(Labels,Labels2Incl,OutDir,Tag,...
    'reorder',reorder,...
    'linkagetype','weighted',...
    'LabelsAbbrev',LabelsAbbrev,...
    'LabelsNames',LabelsNames,...
    'Pmethod','fdr',...
    'AllTrialData', AllTrialData,'saveFig',...
    'PhaseStart', PhaseStart, 'PhaseDur', PhaseDur, 'BaselineWindow', BaselineWindow);


%% Correlation Matrix - Limited to 8 conditions

Labels2Incl = {'NancyPressCheek','NancyRubCheek',...
    'TysonPressCheek','TysonRubCheek',...
    'NancyPressShoulder','NancyRubShoulder',...
    'TysonPressShoulder','TysonRubShoulder'};
LabelsAbbrev = {'NPrC','NRbC',...
    'TPrC','TRbC',...
    'NPrS','NRbS',...
    'TPrS','TRbS'};

reorder = [1 2 5 6 3 4 7 8];

LabelsNames = {{'NPrC','NRbC',...
    'TPrC','TRbC',...
    'NPrS','NRbS',...
    'TPrS','TRbS'},{'Cheek', 'Shoulder'}};

Analyze.FaceScratch3.CorrMatrixAction(Labels,Labels2Incl,OutDir,Tag,...
    'reorder',reorder,...
    'linkagetype','weighted',...
    'LabelsAbbrev',LabelsAbbrev,...
    'LabelsNames',LabelsNames,...
    'Pmethod','fdr',...
    'AllTrialData', AllTrialData,'saveFig',...
    'PhaseStart', PhaseStart, 'PhaseDur', PhaseDur, 'BaselineWindow', BaselineWindow);

%% Cross-format population correlation

Analyze.FaceScratch3.CrossCorr(OutDir,Tag,...
    'Pmethod','fdr',...
    'AllTrialData', AllTrialData,'saveFig',...
    'PhaseStart', PhaseStart, 'PhaseDur', PhaseDur, 'BaselineWindow', BaselineWindow);

% %% Cross-format population correlation - 8 conditions only
% 
% Analyze.FaceScratch3.CrossCorrLimitedConditions(OutDir,Tag,...
%     'Pmethod','fdr',...
%     'AllTrialData', AllTrialData,'saveFig',...
%     'PhaseStart', PhaseStart, 'PhaseDur', PhaseDur, 'BaselineWindow', BaselineWindow);


%% Decode Analysis (Confusion Matrix)

TimeWindow = [0.5 2.5];
Labels2Incl = {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek',...
    'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek',...
    'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder',...
    'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'};
% LabelsAbbrev = {{'Felt', 'Observed'},{'Cheek','Shoulder'}};

LabelsAbbrev = {{'Felt Pinch','Felt Press','Felt Rub','Felt Tap',...
    'Obs Pinch','Obs Press','Obs Rub','Obs Tap',...
    'Felt Pinch','Felt Press','Felt Rub','Felt Tap',...
    'Obs Pinch','Obs Press','Obs Rub','Obs Tap'},{'Cheek', 'Shoulder'}};

Analyze.FaceScratch3.DecodeAnalysisAction(AllTrialData,Labels,Labels2Incl,TimeWindow,...
    OutDir,Tag,'saveFig','LabelsAbbrev',LabelsAbbrev);

%% sliding time window decode analysis (confusion matrix)

TimeWindow = [-.5 5.5];
TimeStep = 0.10;
StepSize = 1.0;
% TimeStep =  0.85;
% StepSize = 1.5;



Labels2Incl = {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek',...
    'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek',...
    'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder',...
    'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'};
LabelsAbbrev = {{'Felt Pinch','Felt Press','Felt Rub','Felt Tap',...
    'Obs Pinch','Obs Press','Obs Rub','Obs Tap',...
    'Felt Pinch','Felt Press','Felt Rub','Felt Tap',...
    'Obs Pinch','Obs Press','Obs Rub','Obs Tap'},{'Cheek', 'Shoulder'}};

DecodeLabels = Labels2Incl;

Analyze.FaceScratch3.SlidingTimeDecodeAction(AllTrialData,...
    Labels2Incl,DecodeLabels,TimeWindow,TimeStep,StepSize,OutDir,Tag,...
    'LabelsAbbrev',LabelsAbbrev,'saveFig');

%% Plot SingleUnitERAs for specific interesting units

Labels = {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek',...
    'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder',...
    'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek',...
    'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder',...
    'Null'};
DecodeLabels = {'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek',...
    'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder',...
    'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek',...
    'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'};

[Units] = Analyze.FaceScratch3.FindSignificantUnits(GlobalTrialCriteria{2},'Labels',Labels,'Labels2Incl',DecodeLabels, 'All');


% clr=[];
% clr(1,:) = repmat([237 0 38]./255,1,1);
% clr(2,:) = repmat([31 120 180]./255,1,1);
% clr(3,:) = repmat([77 175 74]./255,1,1);
% clr(4,:) = repmat([152 78 163]./255,1,1);

% LabelsAbbrev = {'FeelCheek','FeelShoulder','ObsCheek','ObsShoulder'};

PlotConds = {{'Condition','NancyPinchCheek'},...
    {'Condition','NancyPressCheek'},...
    {'Condition','NancyRubCheek'},...
    {'Condition','NancyTapCheek'},...
    {'Condition','NancyPinchShoulder'},...
    {'Condition','NancyPressShoulder'},...
    {'Condition','NancyRubShoulder'},...
    {'Condition','NancyTapShoulder'},...
    {'Condition','TysonPinchCheek'},...
    {'Condition','TysonPressCheek'},...
    {'Condition','TysonRubCheek'},...
    {'Condition','TysonTapCheek'},...
    {'Condition','TysonPinchShoulder'},...
    {'Condition','TysonPressShoulder'},...
    {'Condition','TysonRubShoulder'},...
    {'Condition','TysonTapShoulder'}};
SubplotIds = [1 1 1 1 3 3 3 3 2 2 2 2 4 4 4 4];

% Units = [1 27 1 1830;
%     1 27 2 1835;...
%     ];

% LabelsAbbrev = {{'Felt Cheek','Felt Shoulder','Obs Cheek','Obs Shoulder'},...
%     {'Pinch','Press','Rub','Tap'}};

Phases = {'CueTarget','Delay','Go'};
TimeWindows = [0 2; 0 1; 0 3];

for i=1:length(Units)
    
    Analyze.FaceScratch3.PlotSingleUnitERAsAction(AllTrialData,PlotConds,Phases,TimeWindows,OutDir,...
        'SubplotIds',SubplotIds,'Units',Units(i,:),'FigSize',[4 8],...
        'TaskLabel', Tag,'saveFig');
end


%% Plot cross decode analysis, looking across BP and person.

% to plot along with shuffled formats, include 'withShuffle' in the
% varargin.

% TimeWindow = [0.5 2.5];
TimeWindow = [0.5 3];
% to account for the differences in times, use offsets
% TimeWindow={[0.5 3.15],[0.05 3.35],[0.75 2.85],[0.75 3.05]};
NumRandIterations = 20;
% NumRandIterations = 5;
NumRandUnits = 50;
% TimeWindow = [-1 3];
 
Conditions = {{'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek'},...
    {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'},...
    {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'},...
    {'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'}};

DimNames = {'PersonBP','Action'};
% SplitValNames = {{'NancyCheek','TysonCheek','NancyShoulder','TysonShoulder'},{'Pinch','Press','Rub','Tap'}}; % Label names, correspondign to splits (1,2, etc)
Analyze.FaceScratch3.CrossDecodeAnalysisActionBaseline(AllTrialData,Labels,TimeWindow,OutDir,Tag,...
    'DimNames',DimNames,'saveFig',...
    'Conditions',Conditions,'NumRandIterations',NumRandIterations,'NumRandUnits',NumRandUnits);

%% Cross-Mahlanobis Distance

% to plot along with shuffled formats, include 'withShuffle' in the
% varargin.

% TimeWindow = [0.5 2.5];
TimeWindow = [0.5 3];
% TimeWindow = [-1 3];
 
Conditions = {{'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek'},...
    {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'},...
    {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'},...
    {'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'}};

DimNames = {'PersonBP','Action'};
% SplitValNames = {{'NancyCheek','TysonCheek','NancyShoulder','TysonShoulder'},{'Pinch','Press','Rub','Tap'}}; % Label names, correspondign to splits (1,2, etc)
Analyze.FaceScratch3.CrossMahalanobisDistance(AllTrialData,Labels,TimeWindow,OutDir,Tag,...
    'DimNames',DimNames,'saveFig',...
    'Conditions',Conditions);

%% Cross format classification Dynamic


TimeWindow = [-0.5:.2:4.5];
% TimeWindow = [-.1:.1:4.5];

Conditions = {{'NancyPinchCheek','NancyPressCheek','NancyRubCheek','NancyTapCheek'},...
    {'NancyPinchShoulder','NancyPressShoulder','NancyRubShoulder','NancyTapShoulder'},...
    {'TysonPinchCheek','TysonPressCheek','TysonRubCheek','TysonTapCheek'},...
    {'TysonPinchShoulder','TysonPressShoulder','TysonRubShoulder','TysonTapShoulder'}};

DimNames = {'PersonBP','Action'};
% SplitValNames = {{'NancyCheek','TysonCheek','NancyShoulder','TysonShoulder'},{'Pinch','Press','Rub','Tap'}}; % Label names, correspondign to splits (1,2, etc)
Analyze.FaceScratch3.CrossDecodeAnalysisActionDynamic(AllTrialData,Labels,TimeWindow,OutDir,Tag,...
    'DimNames',DimNames,'sigUnits','saveFig',...
    'Conditions',Conditions,'MinJerk',1);


%% Single unit  - Model Analysis - With Pairs and Triplets

AU.SingleUnitAnalysis(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Tag','BayesModelTwoThree');

%% Single unit  - Model Analysis - WithinFormatR2

AU.SingleUnitAnalysis(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Tag','WithinFormat');

%% Population Model Analysis = Two formats/Pair

Analyze.FaceScratch3.SUP_BayesModelPairs(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Tag', 'BayesModelTwoThree',...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Pmethod','fdr',...
    'SelectionMode','sigUnits',...
    'saveFig',...
    'AllTrialData', AllTrialData);

%% Single unit  - Model Analysis - All 4 formats together

AU.SingleUnitAnalysis(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Tag','BayesModel');

%% Population Model Analysis = All 4 formats together 

Analyze.FaceScratch3.SUP_BayesModel4Format(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Tag', 'BayesModel',...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Pmethod','fdr',...
    'SelectionMode','sigUnits',...
    'saveFig',...
    'AllTrialData', AllTrialData);

%% Population Model Analysis = All 4 formats together - New grouping instead of by triplets

Analyze.FaceScratch3.SUP_BayesModel4Format4Groups(@Analyze.FaceScratch3.Config_LinearAnalysisText,...
    'Tag', 'BayesModel',...
    'Dates',Dates,'Phase',Phase,'Labels',Labels,...
    'GlobalTrialCriteria',GlobalTrialCriteria,...
    'Pmethod','fdr',...
    'SelectionMode','sigUnits',...
    'saveFig',...
    'AllTrialData', AllTrialData);

%% Subspace Analyses

Analyze.FaceScratch3.SC_FormInvarSubsAnal

%% Partial Correlation version 2
% use the version above, in the section "Cross-format population correlation"
Analyze.FaceScratch3.SC_FormatCorr

%% Dynamic Cross Correlation

Analyze.FaceScratch3.SC_FormatCorrDynamic

%% Dynamic Cross Classification (Decode)

Analyze.FaceScratch3.SC_FormatDecodeDynamic
end


