function opts=Config_FormInvarSubsAnal(GlobalOpts, varargin);

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
opts.cvOptions.NReps=10;
opts.cvOptions.NFolds=10;
opts.dec=Predictor.FWClassifier(@Analyze.FaceScratch3.BasicClassifier);
%%

%% Arguements for Analysis function
Analysis2Perform=[1];

opts.AnalysisFCN={@Analyze.FaceScratch3.Dec_SubspaceProj,''};

opts.AnalysisFCN=opts.AnalysisFCN(logical(Analysis2Perform(:)),:);

%% Arguements for population analysis and plotting

if 1==1
    
    opts.PopulationOverwrite=1;
    opts.PlotDates='*';
    opts.PlotSubject='p1';
    
    [varargin,opts.PopulationOverwrite]   = Utilities.ProcVarargin(varargin,'PopulationOverwrite',opts.PopulationOverwrite);
    
    [varargin,opts.PlotDates]   = Utilities.ProcVarargin(varargin,'PlotDates',{opts.PlotDates});
    [varargin,opts.PlotSubject]   = Utilities.ProcVarargin(varargin,'PlotSubject',{opts.PlotSubject});
    
    
    opts.PlotShuffle=0;
    opts.BaselineWindow=[-.9 .1];
    % ConfigureSubject('p1')
    
    opts.PlotSignificance={};
    
    opts.PlotMagnitude={};
    
    Analysis2Perform=[1 0 0 0];
    
    opts.PlotFCN={ ...
        @Analyze.FaceScratch3.PlotPop_SubSpaceProj,{'SubSpaceProjPerActB','SubSpaceProjAvgB'},{};...
        @Analyze.FaceScratch3.PlotPop_SubSpaceProjVer2,{'SubSpaceProjPerActTLL','SubSpaceProjTLL'},{};...
        @Analyze.ActionObsInvarianceText.PlotPop_SubSpaceProjVer2,{'SubSpaceProjPerActTL1','SubSpaceProjTL1'},{'Type','TL1'};...
        @Analyze.ActionObsInvarianceText.PlotPop_SubSpaceProjVer2,{'SubSpaceProjPerActFL' ,'SubSpaceProjFL'} ,{'Type','FL'};...
        % @Analyze.ActionObsInvarianceText.PlotPop_CrossAccuracy,'PopPlot',{};...
        %     @Analyze.ActionObsInvarianceText.PlotPop_CrossAccuracyThroughTime,'PopPlotTime',{};...
        };
    
    opts.PlotFCN=opts.PlotFCN(logical(Analysis2Perform(:)),:);
    
end

