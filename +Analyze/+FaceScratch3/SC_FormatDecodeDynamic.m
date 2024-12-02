
%% Global  Default Options

GlobalOpts.Subject='p1';%'p1','p3';
GlobalOpts.TaskName='FaceScratch3'; % used for folders and file names
GlobalOpts.TaskFileName='FaceScratch3'; % To use with LoadConvertedData

% Analyses phases and windows for trial response and temporally resovled
% response (TR)
GlobalOpts.Phases = {'Go'};
GlobalOpts.windowStarts.('Go')=[.5];
GlobalOpts.windowDuration.('Go')=2.000;

GlobalOpts.windowStartsTR.('Gos')=[-.5:.1:3];   
GlobalOpts.windowDurationTR.('Go')=.5;

% constraints on units.
GlobalOpts.unitType = 'UnitIDsSorted';
GlobalOpts.FRthreshold=.5;

% addtional GlobalOpts.s
GlobalOpts.overwrite=1;
GlobalOpts.shuffle=0;
GlobalOpts.NShuffles=200;
GlobalOpts.saveResults=1;
GlobalOpts.PlotResults=1;

GlobalOpts.PopulationOverwrite=0;
GlobalOpts.PlotDates='*';
GlobalOpts.PlotSubject=GlobalOpts.Subject;
GlobalOpts.unitWFDPThresh=[];

GlobalOpts.RunInline=1;
GlobalOpts.RunBatch=0;
GlobalOpts.RunPlots = 1;

env.set('subject',GlobalOpts.Subject);

Dates = {'20171215','20171220','20180105','20180110','20180112','20180115','20180117'};
GlobalOpts.Dates=Dates;

GlobalOpts.GlobalTrialCriteria = {'TaskLabel','XActionType'};
GlobalOpts.Tag = GlobalOpts.GlobalTrialCriteria{2};


%%
%% Within Format Analyses (e.g. Fig 4 etc.)
GlobalOpts.AnalysesName='FormatDecodeDynamic';
GlobalOpts.PlotShuffle=0;
Dates=GlobalOpts.Dates;
%%
if GlobalOpts.RunInline
    for i=1:length(Dates)        
         AU.PopulationAnalysisThroughTime(@Analyze.FaceScratch3.Config_FormatDecodeDynamic,GlobalOpts,'Dates',Dates(i),'overwrite',1);
    end  
end

%%
if GlobalOpts.RunPlots
  AU.PlotPopulationAnalThroughTime(@Analyze.FaceScratch3.Config_FormatDecodeDynamic ,GlobalOpts,'PopulationOverwrite',1);    
end
