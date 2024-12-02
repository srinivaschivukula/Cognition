function [ output_args ] = PlotSingleUnitERAsAction(AllTrialData,PlotConds,Phases,TimeWindows,OutDir,varargin)
%PLOTSINGLEUNITERAS Plots ERAs for specified units and conditions.
%   Units: Matrix of unit IDs (e.g. ASF). All by default
%   PlotConds: Cell array containing subselection criteria use for
%       extracting specific conditions. One cell per set of conditions to
%       plot.
%   Phases: Which phases to plot over.
%   TimeWindows: Timewindows of the phases (start and end times). Each row
%       should match a phase. Times are relative  to phase start (0) and in
%       seconds.
%   SubplotIDs: Vector of index labels (same dim as PlotConds) indicating
%       which subplot to draw the corresponding PlotCond element to. Subplots
%       are arranged as rows. Draws all onto the same plot by default.
%       Indices should start at 1.
%   Colormap: Colormap to draw plots with. Line by default.
%   Label: Label for saving file name. Empty by default.

[varargin,Units] = Utilities.ProcVarargin(varargin,'Units',[]);
[varargin,SubplotIDs] = Utilities.ProcVarargin(varargin,'SubplotIDs',[]);
[varargin,Colormap] = Utilities.ProcVarargin(varargin,'Colormap',[]);
[varargin,Label] = Utilities.ProcVarargin(varargin,'Label','');
[varargin,UnitLabels] = Utilities.ProcVarargin(varargin,'UnitLabels',[]);
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});


[varargin,SmoothingKernel]=Utilities.ProcVarargin(varargin,'SmoothingKernel', struct('mode','MinJerk','value',.5));
[varargin,TimeStep]=Utilities.ProcVarargin(varargin,'TimeStep',.05);
[varargin,SubSelectPhases] = Utilities.ProcVarargin(varargin,'SubSelectPhases',1);
[varargin,FigSize] = Utilities.ProcVarargin(varargin,'FigSize',[]);


% if ~exist(OutDir,'dir'); mkdir(OutDir); end

% All available (sorted) units by default
ASF = Analyze.returnUniqueFieldValues(AllTrialData,'UnitIDsSorted');
Dates = Analyze.returnFieldValues(AllTrialData,'Date');
if isempty(Units)
    Units = ASF;
    if isstr(UnitLabels)
        tmp = UnitLabels;
        UnitLabels = cell(1,size(Units,1));
        UnitLabels(:) = {tmp};
    end
elseif any(Units(:,4)/10000000 > 1) % Can also specify units in format array-chan-unit-date, convert if necessary
    warning('Converting from date to relative day number.');
    relday = Blackrock.Helper.date2unitId(Units(:,4));
    Units(:,4) = relday;
end

% All plots on same subplot by default
if isempty(SubplotIDs)
    SubplotIDs = ones(size(PlotConds));
end

if isempty(Label)
    SavePre = '%s';
else
    SavePre = ['%s-' Label];
end


NumSubPlots = length(unique(SubplotIDs));
% NumPhases = length(Phases);
TotalDuration = sum(diff(TimeWindows,1,2));
%%
StyleArgs = {'FontName','Arial','FontSize',9,'FontWeight','bold'};
if ~isempty(FigSize)
    figH = plt.fig('units','inches','width',FigSize(1),'height',FigSize(2),'font','Arial','fontsize',9);
else
    figH = plt.fig('units','inches','width',.8*TotalDuration,'height',4*NumSubPlots,'font','Arial','fontsize',9);
end

DateStr = Blackrock.Helper.unitId2date(Units);

%% For each unit, plot the ERAs for the specified conditions.
for cUnitIdx = 1:size(Units,1) % each unit
    cUnit = Units(cUnitIdx,:);
    [~,i] = ismember(ASF,cUnit,'rows');
    i = find(i);
    
    cUnitStr = sprintf('%s-%d-%d-%d',DateStr{cUnitIdx},cUnit(1),cUnit(2),cUnit(3));
    %     cUnitStr = sprintf('%s-%d-%d-%d',DateStr,cUnit(1),cUnit(2),cUnit(3));
    fprintf('Plotting ERA for %s (%d/%d)\n',cUnitStr,cUnitIdx,size(Units,1));
    
    % Set up figure, reusing it.
    figure(figH);
    clf; pnl = panel(); pnl.margin=20;
    pnl.pack(NumSubPlots,1);
    plotERAs(AllTrialData,cUnit,Phases,TimeWindows,PlotConds,SubplotIDs,...
        pnl,SmoothingKernel,TimeStep,StyleArgs,'SubSelectPhases',SubSelectPhases,'LabelsAbbrev',LabelsAbbrev,...
        'Colormap',Colormap);
    
    % Save
    uiLabel = [SavePre '-%d'];
    if cUnitIdx <= length(UnitLabels) && ~isempty(UnitLabels{cUnitIdx})
        filename = sprintf([SavePre '-' UnitLabels{cUnitIdx}],cUnitStr);
        uiLabel = sprintf([uiLabel '-' UnitLabels{cUnitIdx}],cUnitStr,cUnit(4));
    else
        filename = sprintf(SavePre,cUnitStr);
        uiLabel = sprintf(uiLabel,cUnitStr,cUnit(4));
    end
    
    
    
    txtHandle=uicontrol('Style', 'text',...
        'String', uiLabel,... %replace something with the text you want
        'Units','centimeters',...
        'Position', [0 0 length(uiLabel)/4 .5]);
    
    OutDir2 = fullfile(OutDir,'ERA');
    if ~exist(OutDir2,'dir'); mkdir(OutDir2); end
    plt.SaveFigure(saveFig, OutDir2,[filename],'PNG','SVGI');
    
end






end

function plotERAs(AllTrialData,unit,Phases,TimeWindows,PlotConds,SubplotIDs,pnl,SmoothingKernel,TimeStep,StyleArgs,varargin)
[varargin,SubSelectPhases] = Utilities.ProcVarargin(varargin,'SubSelectPhases',1);
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});
[varargin,Colormap] = Utilities.ProcVarargin(varargin,'Colormap',[]);

%Plots the ERAs for the specified unit and conditions into the current plot
NumSubPlots = length(unique(SubplotIDs));
UnitDate = Blackrock.Helper.unitId2date(unit);
for cPlotIdx = 1:NumSubPlots % each subplot
    cPlotConds = PlotConds(SubplotIDs == cPlotIdx);
    
    pnl(cPlotIdx,1).select();
    
    TimeCenters = [];
    BinnedData = cell(length(Phases),length(cPlotConds));
    PlotData = cell(1,length(cPlotConds));
    for ii = 1:length(Phases)
        cPhase = Phases{ii};
        cTW = TimeWindows(ii,1):TimeStep:TimeWindows(ii,2);
        
        switch SmoothingKernel.mode
            case 'MinJerk'
                smoothArg=[SmoothingKernel.value TimeStep];
            case 'Box'
                smoothArg= SmoothingKernel.value;
        end
        
        for jj = 1:length(cPlotConds)
            if SubSelectPhases
                %                 varargin = [{'Date',UnitDate,'Phase',cPhase} cPlotConds{jj}];
                varargin = ['Date',UnitDate,'Phase',cPhase,cPlotConds{jj}];
            else
                varargin = ['Date',UnitDate,cPlotConds{jj}];
            end
            cTrialData = Analyze.SubSelectTrials(AllTrialData,varargin{:});
            BinnedData{ii,jj} = Analyze.getNeuralData(cTrialData,unit,cTW,'smooth',smoothArg);
        end
        
        if isempty(TimeCenters)
            TimeCenters = cTW;
        else
            % Remove intersecting value to avoid redundancy
            TimeCenters = [TimeCenters(1:end-1) (TimeCenters(end) + cTW)];
        end
    end
    
    % Pick out center of each time bin (removing last one since not in
    % bounds)
    TimeCenters = TimeCenters(1:end-1) + TimeStep/2;
    
    %% Plot
    % BinnedData{ii,jj} is Trials x TimeBins. Concat timebins
    for jj = 1:length(cPlotConds)
        PlotData{jj} = horzcat(BinnedData{:,jj});
    end
    TotalDuration = sum(diff(TimeWindows,1,2));
    %     timeCenters = linspace(TimeWindows(1,1),TotalDuration,size(PlotData{1},2));
    %     timeCenters = TimeWindows(1,1):TimeStep:(TotalDuration-TimeWindows(1,1));
    Labels = getLegendLabels(cPlotConds,PlotData,'LabelsAbbrev',LabelsAbbrev);
    
    if ~isempty(Colormap)
        Analyze.plotEventRelatedAverage(PlotData,Labels,'TimeVec',TimeCenters,'NoLine',...
           'Legend','Colors',Colormap);
    else
        Analyze.plotEventRelatedAverage(PlotData,Labels,'TimeVec',TimeCenters,'NoLine',...
            'Legend','LegendFontSize',9,'Colors',lines(length(cPlotConds)));
    end
    ax1 = gca;
    set(ax1,'box','off','FontWeight','bold','FontSize',9);
    ax1.XRuler.Axle.LineStyle = 'none';
    set(gca,'TickDir','out');
    set(ax1,'TickLength',[0 0]);
    axis  tight;
    xlabel('Time (s)',StyleArgs{:});
    ylabel('FR (Hz)',StyleArgs{:});
    
end

%% Rescale y axes if multiple subplots
if NumSubPlots > 1
    for cPlotIdx = 1:NumSubPlots
        pnl(cPlotIdx,1).select();
        ylimtmp = ylim();
        if cPlotIdx == 1
            ylimmin = ylimtmp(1);
            ylimmax = ylimtmp(2);
        else
            ylimmin = min(ylimmin,ylimtmp(1));
            ylimmax = max(ylimmax,ylimtmp(2));
        end
    end
    
    for cPlotIdx = 1:NumSubPlots
        pnl(cPlotIdx,1).select();
        ylim([ylimmin ylimmax]);
    end
end

%% Line markers for phase starts/ends and label each phase
for cPlotIdx = 1:NumSubPlots
    pnl(cPlotIdx,1).select();
    
    tEvents = [0 cumsum(TimeWindows(:,2))'];
    plt.vline(tEvents,'k');
    
    % Label each phase
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    for ii = 1:length(Phases)
        text2Use = Phases{ii};
        if strcmpi(text2Use,'cuetarget')
            text2Use='Cue';
        end
        text(tEvents(ii),ytext,text2Use,StyleArgs{:});
    end
end


end

function Labels = getLegendLabels(PlotConds,PlotData,varargin)
[varargin,LabelsAbbrev]=Utilities.ProcVarargin(varargin,'LabelsAbbrev',{});

if isempty(PlotConds{1})
    Labels = {'All'};
    return;
end

if isempty(LabelsAbbrev)
    Labels = {};
    for jj = 1:length(PlotConds)
        tempCond = cell(1,length(PlotConds{jj})/2);
        for c = 1:length(PlotConds{jj})/2
            [~,op] = extractNameFromPlotCond(PlotConds{jj}{2*c-1});
            cond = PlotConds{jj}{2*c};
            if isnumeric(cond)
                cond = num2str(cond);
            end
            tempCond{c} = sprintf('%s%s',op,cond);
        end
        cLabel = sprintf('%s',strcat(tempCond{:}));
        Labels = [Labels cLabel];
    end
else
    Labels =unique(LabelsAbbrev{1},'stable');
end


end

function [Name,op] = extractNameFromPlotCond(plotCond)

if strcmp('>=',plotCond(1:2))
    Name = plotCond(3:end);
    op = '>=';
elseif strcmp('<=',plotCond(1:2))
    Name = plotCond(3:end);
    op = '<=';
elseif strcmp('>',plotCond(1))
    Name = plotCond(2:end);
    op = '>';
elseif strcmp('<',plotCond(1))
    Name = plotCond(2:end);
    op = '<';
else
    Name = plotCond;
    op = '';
end

end