function opts = Config_LinearAnalysisText(varargin)
%CONFIG_LINEARANALYSISCUSTOM Calls Config_LinearAnalysis with custom labels
%or globaltrial criteria
[varargin,Labels]   = Utilities.ProcVarargin(varargin,'Labels',{});
[varargin,GlobalTrialCriteria]   = Utilities.ProcVarargin(varargin,'GlobalTrialCriteria',{});

[varargin,PThresh]   = Utilities.ProcVarargin(varargin,'PThresh',0.001);

[varargin,opts.ModelGroups]   = Utilities.ProcVarargin(varargin,'ModelGroups',{});
[varargin,opts.ModelNames]   = Utilities.ProcVarargin(varargin,'ModelNames',{});
[varargin,opts.IsOrder]   = Utilities.ProcVarargin(varargin,'OrderMode');

[varargin,ConditionFieldLabel] = Utilities.ProcVarargin(varargin,'ConditionFieldLabel','Condition');
[varargin,opts.ModelCompare]   = Utilities.ProcVarargin(varargin,'ModelCompare',{});
[varargin,opts.ModelCompareNames]   = Utilities.ProcVarargin(varargin,'ModelCompareNames',{});
[varargin,BaselinePhase]   = Utilities.ProcVarargin(varargin,'BaselinePhase','Delay');
[varargin,BaselineWindow]   = Utilities.ProcVarargin(varargin,'BaselineWindow',[1.5 2.5]);

% [varargin,BaselinePhase]   = Utilities.ProcVarargin(varargin,'BaselinePhase','Go');
% [varargin,BaselineWindow]   = Utilities.ProcVarargin(varargin,'BaselineWindow',[-1.25 -.25]);

[varargin,PhaseStart]   = Utilities.ProcVarargin(varargin,'PhaseStart',1);
[varargin,PhaseDur]   = Utilities.ProcVarargin(varargin,'PhaseDur',3);
[varargin,PlotAvg]   = Utilities.ProcVarargin(varargin,'PlotAvg'); % for plotting average response of units with a certain cond pref
[varargin,Units] = Utilities.ProcVarargin(varargin,'Units',[]);
[varargin,PlotLabel] = Utilities.ProcVarargin(varargin,'PlotLabel',[]);
[varargin,AlignArgs] = Utilities.ProcVarargin(varargin,'AlignArgs',[]);
[varargin,Tag] = Utilities.ProcVarargin(varargin,'Tag',{});
[varargin,SubLabelIdx] = Utilities.ProcVarargin(varargin,'SubLabelIdx',[]);
[varargin,Dates]   = Utilities.ProcVarargin(varargin,'Dates','Dates');
[varargin,PlotDates]   = Utilities.ProcVarargin(varargin,'PlotDates','*');



opts.Units = Units;
opts.PlotLabel = PlotLabel;
opts.AlignArgs = AlignArgs;

% ITI=2; Cue=1.5; Delay=1
%% Arguements for wrapping function
% opts.Dates={ '20160316','20160318'  '20160323','20160401'};
opts.Dates=Dates;
[varargin,Phase]   = Utilities.ProcVarargin(varargin,'Phase','Go');
[varargin,opts.UseSparse]   = Utilities.ProcVarargin(varargin,'UseSparse');
% [varargin,WinSize]   = Utilities.ProcVarargin(varargin,'WinSize',1);
% [varargin,Start]   = Utilities.ProcVarargin(varargin,'Start',.5);

opts.Labels = Labels;
opts.ConditionFieldLabel=ConditionFieldLabel;
opts.Subject='p1';

opts.GlobalTrialCriteria = GlobalTrialCriteria;

opts.limitWFSNR=1;

opts.NShuffles=1000;
opts.RemoveTaskVariance=0;

opts.TaskName='FaceScratch3'; % used for folders and file names
opts.TaskFileName='FaceScratch3'; % To use with LoadConvertedData
if ~isempty(Tag)
    AnalysisName = Tag;
else
    AnalysisName = GlobalTrialCriteria{2};
end
if iscell(AnalysisName)
    AnalysisName = strcat(AnalysisName{:});
end
opts.AnalysisName=AnalysisName;  % used for folders and file names
if opts.IsOrder
    opts.ResultsDir=fullfile(env.get('result'),'FaceScratch3','SUAnal',[opts.AnalysisName '-' Phase],'Order');
end

if strcmpi(opts.AnalysisName,'StimProfiles') & strcmpi(GlobalTrialCriteria{2},'XactionType')
    opts.ResultsDir=fullfile(env.get('result'),'FaceScratch3','SUAnal',[opts.AnalysisName 'Action' '-' Phase]);
else
    opts.ResultsDir=fullfile(env.get('result'),'FaceScratch3','SUAnal',[opts.AnalysisName '-' Phase]);
end

opts.overwrite=1;
opts.shuffle=0;

StepSize=.25;
WinSize=2;

opts.Phases ={Phase};  % Phases to Analyze
opts.windowStarts.(Phase)=PhaseStart;

if strcmpi(GlobalTrialCriteria{2},'xactiontype')
    opts.windowDuration.(Phase)=2; %winsize
else
    opts.windowDuration.(Phase)=3; %winsize
end

opts.unitType='UnitIDsSorted'; % UnitIDsSorted or UnitIDs
opts.FRthreshold=.5;

opts.saveResults=1;
opts.PlotResults=1;
%%
% opts.AnalysisFCN={@Analyze.CueGoNo.SUA_BasicAnalysis,'';
%                   @Analyze.CueGoNo.SUA_BasicAnalysis,'Info'};  % Action & Baseline
% opts.AnalysisFCN={ @AU.SUA_BasicAnalysis,''};  % Action & Baseline
% opts.AnalysisFCN={ @AU.SUA_ModelAnalysis,''};  % Action & Baseline

cLabel = GlobalTrialCriteria{2};

if strcmp(AnalysisName,'XTouchRF-Bilat')
    opts.AnalysisFCN={ @Analyze.FaceScratch.SUA_BilatWithinSide, '' };%...
elseif strcmpi(AnalysisName,'StimProfiles') & strcmpi(cLabel,'XactionType')
    opts.AnalysisFCN={ @Analyze.FaceScratch3.SUA_StimProfilesAction, '' };
elseif strcmpi(AnalysisName,'StimProfiles') & ~strcmpi(cLabel,'XactionType')
    opts.AnalysisFCN={ @Analyze.FaceScratch3.SUA_StimProfiles, '' };
elseif strcmpi(AnalysisName,'BayesModel') & strcmpi(cLabel,'XactionType')
    opts.AnalysisFCN={ @Analyze.FaceScratch3.SUA_BayesModel4Format, '' };
elseif strcmpi(AnalysisName,'BayesModelTwoThree') & strcmpi(cLabel,'XactionType')
    opts.AnalysisFCN={ @Analyze.FaceScratch3.SUA_BayesModelTwoThree, '' };
elseif strcmpi(AnalysisName,'WithinFormat') & strcmpi(cLabel,'XactionType')
    opts.AnalysisFCN={ @Analyze.FaceScratch3.SUA_WithinFormatR2, '' };
elseif strcmpi(AnalysisName,'ModelComparison')
    opts.AnalysisFCN={ @Analyze.FaceScratch3.SUA_ModelCompare, '' };
elseif strcmpi(AnalysisName,'XMotor')
    opts.AnalysisFCN={ @AU.SUA_BasicAnalysis, '' };
else
    opts.AnalysisFCN={ @AU.SUA_BasicAnalysis,'';...
        @Analyze.FaceScratch3.SUA_ModelAnalysis,'';...
        @Analyze.FaceScratch3.SUA_ModelCompare,''};  % Action & Baseline
end

opts.ModelCompareIncludeBaseline=0;
% opts.BaselinePhase='CueTarget';
% % opts.BaselineWindow=[0 1];
% opts.BaselineWindow=[1.5 4];

opts.BaselinePhase=BaselinePhase;
% opts.BaselineWindow=[1.5 2.5];
opts.BaselineWindow=BaselineWindow;
%%
opts.ERADir=fullfile(env.get('result'),opts.TaskName,'ERA',[opts.AnalysisName opts.TaskFileName]);
opts.ERAPanelCondNames={'*'};

% What values from the above fields to plot in each panel row.
opts.ERAPanelCondVals={{'*'}};
% opts.ERAPanelCondVals={{'EarlyVisual','LateVisual','Match2Sample'}};
% What conditions to draw in each plot.

% What conditions to draw in each plot.
EffectorPlotConds1 = {{'Condition','1'},{'Condition','2'},{'Condition','3'},{'Condition','4'},{'Condition','5'},{'Condition','6'}};
% EffectorPlotConds2 = {{'CueGoNo','Go'},{'CueGoNo','NoGo'}};
% opts.ERAPlotConds={EffectorPlotConds1,EffectorPlotConds2};
opts.ERAPlotConds={EffectorPlotConds1};



% opts.ERAGlobalTrialCriteria={'>MovementWaitTime',4.5};
opts.ERAPlotFunction={@()plt.fig('units','inches','width',6,'height',5,'font','Arial','fontsize',9)};

opts.ERASmoothingKernel= struct('mode','BoxWin','value',.6);

opts.plotRaster=1;
%     opts.Phases ={'ITI' 'CueTarget' 'Delay'  'Go1'};  % Phases to Analyze

StepSize=.1;
opts.ERAPlotBins.('CueTarget')=[-.5:StepSize:2.5];
opts.ERAPlotBins.('Delay')=[0:StepSize:3];
opts.ERAPlotBins.('Go')=[-.5:StepSize:4.5];

opts.ERASideMargin= 4;
% opts.plotArgs={'NoLine','Colors',[94 158 90; 158 90 90; 5 255 26; 255 5 10]/255};
opts.plotArgs={'NoLine','Colors',lines(6)};
%% Arguements for called function



% opts.Dates={'20150812' '20150814' '20150817' '20150821','20160427','20160504',...
%     '20160511','20160516','20160520','20160523'};
opts.PopulationOverwrite=0;
% opts.PlotDates={ '20170630'};
opts.PlotDates={'*'};
[varargin,opts.PlotDates]   = Utilities.ProcVarargin(varargin,'PlotDates',opts.PlotDates);

opts.PlotSubject='p1';
opts.PlotShuffle=0;
opts.PlotAvg = PlotAvg;
% ConfigureSubject('p1')

% opts.PlotFields={'BA_pEarlyVisualGo','BA_pEarlyVisualNo','BA_pLateVisualGo','BA_pLateVisualNo',...
%  'BA_pDelayGo','BA_pDelayNo','BA_pTrainingGo','BA_pTrainingNo','BA_pReactionTimeGo','BA_pReactionTimeNo'};

opts.PlotFields={};

opts.AddArg={'>UnitmeanRate',0.5};

if strcmp(AnalysisName,'XTouchRF-Bilat')
    Analysis2Perform=[0 0 0 0 1];
    %     Analysis2Perform=[1 1 1 1 1];
    
    opts.PlotFCN = {...
        @Analyze.FaceScratch.Plot_SideSpecificity	,'SideSpecificityDPLimQuality',{'Labels',opts.Labels,'Meas','dp', 'LimQuality'};...
        @Analyze.FaceScratch.Plot_SideSpecificity	,'SideSpecificityDP',{'Labels',opts.Labels,'Meas','dp'};...
        @Analyze.FaceScratch.Plot_SideSpecificity	,'SideSpecificityR2',{'Labels',opts.Labels,'Meas','R2'};...
        @Analyze.FaceScratch.Plot_PairwiseR2        ,'PairwiseR2',{'Labels',opts.Labels};...
        @Analyze.FaceScratch.Plot_WithinAcross      ,'WithinAcross',{'Labels',opts.Labels};...
        };
    opts.PlotFCN=opts.PlotFCN(logical(Analysis2Perform(:)),:);
elseif strcmpi(AnalysisName,'XMotor')
        Analysis2Perform=[1 0 0 0 0];
        clr=lines(29);

    opts.PlotFCN = {...
        @DataPlot.SUP_BasicTuning	,'Tuning',{'Labels',opts.Labels, 'PThresh', PThresh};...
        @DataPlot.SUP_Specificity,'Specif',{'Labels',opts.Labels, 'PThresh', PThresh};...
        @DataPlot.SUP_ActionCorrelations,'Relat',{'Labels',opts.Labels,'clr',clr,'distanceFCN','corr','PThresh', PThresh,'ThreshType','AnySig','PlotAvg',0};...
        };
    opts.PlotFCN=opts.PlotFCN(logical(Analysis2Perform(:)),:);

else
    % opts.AddArg={[]};
    opts.PlotResults=1; % handle plots internally
    Analysis2Perform=[1 1 1];
    Analysis2Perform=[1 0 1];
    Analysis2Perform=[1 1 1 1];
    Analysis2Perform=[1 1 1 0];
    Analysis2Perform=[1 0 1 0];
    Analysis2Perform=[0 0 1 0];
    % Analysis2Perform=[0 0 0 0 0 1];
    Analysis2Perform=[1 1 1 0 0 0 0];
    Analysis2Perform=[0 0 1 0 0 0 0];
    Analysis2Perform=[0 0 1 0 0 0 0 0 0];
    %     Analysis2Perform=[0 0 0 1 0 0 0 0 0];
    % Analysis2Perform=[1 0 0 0 0 0 0];
    % Analysis2Perform=[1 0 0 1 0 0 0];
    % Analysis2Perform=[0 0 0 1 0 0 0];
    % Analysis2Perform=[1 1 1 0 0 0 0];
    Analysis2Perform=[1 1 1 1 1 1 1 1];
    
    Analysis2Perform=[1 1 1 1 0 0 0 1];
    Analysis2Perform=[1 1 1 1 1 1 1];
    
    clr=[lines(5)*1.05; flipud(lines(5))*.85];
    % bc=[1 0 0; 1 1 0 ; 0 1 0;0 1 1 ;0 0 1]*.75;
    % clr=[bc; flipud(bc)*.75];
    
    opts.PlotFCN={...
        @DataPlot.SUP_BasicTuning	,'Tuning',{'Labels',opts.Labels,'PThresh',PThresh};...
        @DataPlot.SUP_Specificity,'Specif',{'Labels',opts.Labels,'PThresh',PThresh};...
        @DataPlot.SUP_ActionCorrelations,'Relat',{'Labels',opts.Labels,'distanceFCN','corr','PThresh',PThresh,'ThreshType','AnySig','SubLabelIdx',SubLabelIdx,'PerDay',0};...
        @DataPlot.SUP_ActionCorrelations,'RelatImag',{'Labels',opts.Labels,'distanceFCN','corr','PThresh',PThresh,'ThreshType','AnySig','SubLabelIdx',SubLabelIdx,'PerDay',1};...
        %        @DataPlot.SUP_ActionCorrelations,'Relat',{'Labels',opts.Labels,'distanceFCN','corr','PThresh',0.001,'ThreshType','AnySig'};...
        %         @DataPlot.SUP_ActionCorrelations,'Relat',{'Labels',opts.Labels,'distanceFCN','euc','PThresh',0.001,'ThreshType','AnySig'};...
        %         @Analyze.FaceScratch.CorrelationsSubtracted,'Relat',{'Labels',opts.Labels,'clr',clr,'distanceFCN','corr','PThresh',0.1,'ThreshType','AnySig'};...
        @Analyze.FaceScratch3.SUP_ModelAnalysis,'Models',{'Labels',opts.Labels,'PThresh',PThresh};...
        @Analyze.FaceScratch3.SUP_ModalityAnalysis,'Models',{'Labels',opts.Labels,'PThresh',PThresh};...
        %         @Analyze.FaceScratch3.SUP_Other,'Models',{'Labels',opts.Labels,'PThresh',PThresh};...
        %                 @Analyze.FaceScratch.SUP_DecodeAnalysis,'Models',{'Labels',opts.Labels,'PThresh',PThresh};...
        @Analyze.FaceScratch3.SUP_ModelCompare,'Models',{'Labels',opts.Labels,'PThresh',PThresh};...
        };
    
    
    opts.PlotFCN=opts.PlotFCN(logical(Analysis2Perform(:)),:);
end
opts.PlotSaveType={'PNG','SVG'};
end

