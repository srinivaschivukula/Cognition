function opts=Config_FormatCorr(GlobalOpts, varargin);

opts=GlobalOpts;
opts.PlotShuffle=0;
[varargin,opts.Dates]   = Utilities.ProcVarargin(varargin,'Dates',{opts.Dates});
[varargin,opts.shuffle]   = Utilities.ProcVarargin(varargin,'Shuffle',opts.shuffle);
[varargin,opts.overwrite]   = Utilities.ProcVarargin(varargin,'overwrite',opts.overwrite);
[varargin,opts.NShuffles]   = Utilities.ProcVarargin(varargin,'NShuffles',opts.NShuffles);
[varargin,opts.unitType]   = Utilities.ProcVarargin(varargin,'unitType',opts.unitType);
[varargin,opts.FRthreshold]   = Utilities.ProcVarargin(varargin,'FRthreshold',opts.FRthreshold);

if opts.shuffle; opts.AnalysisName=['Shuff' opts.AnalysisName];opts.PlotShuffle=1;else;  end


opts.ResultsDir=fullfile(env.get('result'),opts.TaskName,'SUAnal',opts.AnalysesName);


opts.NumShuffles=200;
opts.cvOptions.ValidationType='ClassicCrossValidation';
opts.cvOptions.NReps=30;
opts.cvOptions.NFolds=10;
opts.dec=Predictor.FWClassifier(@Analyze.FaceScratch3.BasicClassifier);
%%

%% Arguements for Analysis function
Analysis2Perform=[1];

opts.AnalysisFCN={@Analyze.FaceScratch3.PopAnal_FormatCorrShuff,''};

opts.AnalysisFCN=opts.AnalysisFCN(logical(Analysis2Perform(:)),:);

%% Arguements for population analysis and plotting

[varargin,opts.PopulationOverwrite]   = Utilities.ProcVarargin(varargin,'PopulationOverwrite',opts.PopulationOverwrite);
[varargin,opts.PlotDates]   = Utilities.ProcVarargin(varargin,'PlotDates',{opts.PlotDates});
[varargin,opts.PlotSubject]   = Utilities.ProcVarargin(varargin,'PlotSubject',{opts.PlotSubject});
[varargin,opts.unitWFDPThresh]   = Utilities.ProcVarargin(varargin,'unitWFDPThresh',[]);

% ConfigureSubject('p1')

opts.PlotSignificance={};
opts.PlotMagnitude={};

Analysis2Perform=[1 0];

    opts.PlotFCN={ ...
        @Analyze.FaceScratch3.PlotPop_FormatCorr ,'CorrThroughTime',{};...
        @Analyze.FaceScratch3.PlotPop_FormatCorrThroughTime ,'CorrThroughTime',{};...
        };
opts.PlotFCN=opts.PlotFCN(logical(Analysis2Perform(:)),:);