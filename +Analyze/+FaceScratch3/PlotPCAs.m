function [ output_args ] = PlotPCAs(AllTrialData,PlotConds,Phases,TimeWindows,OutDir,varargin)
%PLOTPCAS Plots projection of PCAs over time, separated by conditions.
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
%   TotalVar: Total variance components should add up to.
%   PThresh: Pval threshold to use. p=0.001 by default.
[varargin,Units] = Utilities.ProcVarargin(varargin,'Units',[]);
[varargin,PThresh] = Utilities.ProcVarargin(varargin,'PThresh',0.001);
[varargin,SubplotIDs] = Utilities.ProcVarargin(varargin,'SubplotIDs',[]);
[varargin,Colors] = Utilities.ProcVarargin(varargin,'Colors',[]);
[varargin,Label] = Utilities.ProcVarargin(varargin,'Label','');
[varargin,NumPC] = Utilities.ProcVarargin(varargin,'NumPC',20);
[varargin,Normalize] = Utilities.ProcVarargin(varargin,'Normalize',0);
[varargin,PlotERA] = Utilities.ProcVarargin(varargin,'PlotERA');
[varargin,Plot2D] = Utilities.ProcVarargin(varargin,'Plot2D');
[varargin,PCPairsMaxIdx] = Utilities.ProcVarargin(varargin,'PCPairsMaxIdx',5);
[varargin,Animate] = Utilities.ProcVarargin(varargin,'Animate');
[varargin,Frames2Save] = Utilities.ProcVarargin(varargin,'Frames2Save',{});
[varargin,Groups] = Utilities.ProcVarargin(varargin,'Groups',[]);
[varargin,NewFigForEach] = Utilities.ProcVarargin(varargin,'NewFigForEach');
[varargin,saveFig] = Utilities.ProcVarargin(varargin,'saveFig');


[varargin,SmoothingKernel]=Utilities.ProcVarargin(varargin,'SmoothingKernel', struct('mode','MinJerk','value',.5));
[varargin,TimeStep]=Utilities.ProcVarargin(varargin,'TimeStep',.05);


% if ~exist(OutDir,'dir'); mkdir(OutDir); end

% All available (sorted) units by default
ASF = Analyze.returnUniqueFieldValues(AllTrialData,'UnitIDsSorted');
Dates = Analyze.returnUniqueFieldValues(AllTrialData,'Date');
if isempty(Units)
    Units = ASF;
elseif any(Units(:,4)/10000000 > 1) % Can also specify units in format array-chan-unit-date, convert if necessary
    warning('Converting from date to relative day number.');
    relday = Blackrock.Helper.date2unitId(Units(:,4));
    Units(:,4) = relday;
end

fprintf('Using %d units for PCA\n',size(Units,1));

% All plots on same subplot by default
if isempty(SubplotIDs)
    SubplotIDs = ones(size(PlotConds));
end

if isempty(Label)
    SavePre = '%s';
else
    SavePre = [Label '-%s'];
end

if isempty(Groups)
    Groups = 1:length(PlotConds);
end

NumSubPlots = length(unique(SubplotIDs));
% NumPhases = length(Phases);
TotalDuration = sum(diff(TimeWindows,1,2));

% DateStr = Blackrock.Helper.unitId2date(Units);

switch SmoothingKernel.mode
    case 'MinJerk'
        smoothArg=[SmoothingKernel.value TimeStep];
    case 'Box'
        smoothArg= SmoothingKernel.value;
end

%% Get neural data for units, grouped by condition.
TimeCenters = [];
fprintf('Loading data...\n');
for ii = 1:length(Phases)
    cPhase = Phases{ii};
    cTW = TimeWindows(ii,1):TimeStep:TimeWindows(ii,2);
    
    for dd = 1:length(Dates)
        cDate = Dates{dd};
        dateCode = Blackrock.Helper.date2unitId(cDate);
        cUnits = Units(Units(:,4)==dateCode,:);
    
        for cc = 1:length(PlotConds)
            cPlotCond = PlotConds{cc};
            
            varargin = [{'Date',cDate,'Phase',cPhase} cPlotCond{:}];
%             varargin = [{'Phase',cPhase} cPlotCond{:}];
            cTrialData = Analyze.SubSelectTrials(AllTrialData,varargin{:});
            BinnedData{ii,cc,dd} = Analyze.getNeuralData(cTrialData,cUnits,cTW,'smooth',smoothArg);
        end
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

%% Concatenate data into a matrix that we can do PCA on.

% BinnedData{ii,cc,dd} is Trials x TimeBins x NUnits.
% Concat across days first.
for ii = 1:length(Phases)
    for cc = 1:length(PlotConds)
        BinnedData2{ii,cc} = cat(3,BinnedData{ii,cc,:});
    end
end

% Concat timebins into PlotData{cc} (Trials x TimeBins x NUnits)
for cc = 1:length(PlotConds)
    PlotData{cc} = horzcat(BinnedData2{:,cc});
end

% Average across trials. PlotDataAvg is TimeBins x NUnits
for cc = 1:length(PlotConds)
    PlotDataAvg{cc} = squeeze(mean(PlotData{cc},1));
end

% Now reconfigure into Data matrix: (TimeBins*Conditions) x NUnits
Data = vertcat(PlotDataAvg{:}); 

% Normalize FR of each unit independent by the maxmin (averaging within
% condition) (i.e. so FR ranges from 0 to 1).
if Normalize
    datamin = min(Data);
    datamaxmin = max(Data)-datamin;
    Data = (Data-repmat(datamin,size(Data,1),1))./repmat(datamaxmin,size(Data,1),1);
    minBig(1,:,:) = datamin;
    maxminBig(1,:,:) = datamaxmin;
    minBig = repmat(minBig,size(PlotData{1},1),size(PlotData{1},2),1);
    maxminBig = repmat(maxminBig,size(PlotData{1},1),size(PlotData{1},2),1);
    for cc = 1:length(PlotConds)
        PlotData{cc} = (PlotData{cc}-minBig)./maxminBig;
    end
end

% Remove bad unit data in case (nans)
GoodIdx = ~any(isnan(Data),1);
if any(~GoodIdx)
    warning('Some bad units found (nans)...removing...');
    Data = Data(:,GoodIdx);
end

%% PCA
% coeff is components (ncomponents x ncomponents = nunits x nunits)
% score is coefficients/weights/projection (TimeBins*Conditions x nunits)
% score*coeff'
Data2 = bsxfun(@minus,Data,mean(Data));
[coeff,scores,~,~,explained] = pca(Data2);



%% Project neural data onto PCA space 

ProjData = cell(1,NumPC);
for pcIdx = 1:NumPC
    cPC = coeff(:,pcIdx);
    for cc = 1:length(PlotConds)
        for i = 1:size(PlotData{cc},1)
            cNeural = squeeze(PlotData{cc}(i,:,GoodIdx));
            ProjData{pcIdx}{cc}(i,:) = (cNeural*cPC)';
        end
    end
end

StyleArgs = {'FontName','Arial','FontSize',14};
%% Plot each component like with ERAs
if PlotERA
    if ~NewFigForEach
        figH = plt.fig('units','inches','width',1.2*TotalDuration,'height',4*NumSubPlots,'font','Arial','fontsize',9);
    end
    
    for cPCIdx = 1:NumPC % each PC
        if NewFigForEach
            figH = plt.fig('units','inches','width',1.2*TotalDuration,'height',4*NumSubPlots,'font','Arial','fontsize',9);
        end
        
        cPCData = ProjData{cPCIdx};
        
        cPCStr = sprintf('pc-%d',cPCIdx);
        fprintf('Plotting ERA for %s (%d/%d)\n',cPCStr,cPCIdx,NumPC);
        
        % Set up figure, reusing it.
        figure(figH);
        clf; pnl = panel(); pnl.margin=20;
        pnl.pack(NumSubPlots,1);
        plotERAs(cPCData,TimeCenters,Phases,TimeWindows,PlotConds,SubplotIDs,pnl,StyleArgs,'Colors',Colors);
        
        % Save
        filename = sprintf(SavePre,cPCStr);
        if Normalize
            filename = [filename '-norm'];
        end
        uiLabel = sprintf([SavePre ': var expl = %.2f%%'],cPCStr,explained(cPCIdx));
        
        txtHandle=uicontrol('Style', 'text',...
            'String', uiLabel,... %replace something with the text you want
            'Units','centimeters',...
            'Position', [0 0 length(uiLabel)/4 .5]);
        
        plt.SaveFigure(saveFig,fullfile(OutDir,'PCA'),filename,'PNG','SVGI');
    end
end

%% Plot pairs of PCAs in animated movie.
% Also save frame from middle of each phase.
if Plot2D
    pairs2plot = combnk(1:PCPairsMaxIdx,2);
    
    figH = plt.fig('units','inches','width',8,'height',8,'font','Arial','fontsize',9);
    
    for p = 1:size(pairs2plot,1)
        %%
        figure(figH);
        clf; pnl = panel(); pnl.margin=20;
        pnl.pack(1,1);
        
        pc1Idx = pairs2plot(p,1);
        pc2Idx = pairs2plot(p,2);
        
        pc1 = ProjData{pc1Idx};
        pc2 = ProjData{pc2Idx};
        
        pc1Mean = cellfun(@mean,pc1,'UniformOutput',false);
        pc2Mean = cellfun(@mean,pc2,'UniformOutput',false);
        
%         pc1Traj = vertcat(pc1Mean{:})';
%         pc2Traj = vertcat(pc2Mean{:})';
        
        xstr = sprintf('PC %d',pc1Idx);
        ystr = sprintf('PC %d',pc2Idx);
        tstr = sprintf('Trajectories for PCs %d and %d',pc1Idx,pc2Idx);
        
        SaveFilePref = sprintf('Traj2D-pc%d+%d',pc1Idx,pc2Idx);
        MovieFrames = plotMovie(pc1Mean,pc2Mean,xstr,ystr,tstr,...
            Phases,TimeWindows,TimeCenters,'StyleArgs',StyleArgs,...
            'Animate',Animate,'Frames2Save',Frames2Save,...
            'SaveDir',fullfile(OutDir,'PCA'),'SaveFilePref',SaveFilePref);
        
        %% Save movie
        vidName = fullfile(OutDir,'PCA',[SaveFilePref]);
        vidObj = VideoWriter(vidName,'MPEG-4');
        open(vidObj);
        writeVideo(vidObj,MovieFrames);
        close(vidObj);
    end
end
end

function MovieFrames = plotMovie(x,y,xstr,ystr,tstr,Phases,TimeWindows,TimeCenters,varargin)
% x and y each contain one cell per condition
[varargin,StyleArgs] = Utilities.ProcVarargin(varargin,'StyleArgs',{});
% [varargin,Animate] = Utilities.ProcVarargin(varargin,'Animate',0);
[varargin,Frames2Save] = Utilities.ProcVarargin(varargin,'Frames2Save',[]);
[varargin,SaveDir] = Utilities.ProcVarargin(varargin,'SaveDir','.');
[varargin,SaveFilePref] = Utilities.ProcVarargin(varargin,'SaveFilePref','Traj2D');


Phases2 = {'', Phases{:}};
%%
% figure(figH);
% clf; pnl = panel(); pnl.margin=20;
% pnl.pack(1,1);

tEvents = [0 cumsum(TimeWindows(:,2))'];
clr = lines(length(x));
clr = jet(length(x));
% clr = mat2cell(clr,[1 1 1 1],3);
% First plot just the trajectories
for i = 1:length(x)
    hTraj=plot(x{i},y{i},'Color',clr(i,:));
    hold on;
end
% axis square;
xlabel(xstr,StyleArgs{:});
ylabel(ystr,StyleArgs{:});
title({tstr,''},StyleArgs{:});
set(gca,'XTick',[],'YTick',[]);
xruler = get(gca,'XRuler');

% Hides border boxes
set(xruler.Axle,'Visible','off');  
yruler = get(gca,'YRuler');
set(yruler.Axle,'Visible','off');  
box off;

% % Changes backgroudn color
% bg = [1 1 1]*.95;
% set(gcf,'Color',bg);
% set(gca,'Color',bg);

% Save an image of just the trajectories, without any markers
plt.SaveFigure(1,SaveDir,[SaveFilePref],'PNG','SVGI');
%%

% [h,icons] = legend(PlotLabels,'Location','SouthEast');
% set(h,StyleArgs{:});
% for i = 1:length(icons)
%     if ~isempty(strfind(class(icons(i)),'Text'))
%         set(icons(i),StyleArgs{:});
%     elseif ~isempty(strfind(class(icons(i)),'Line'))
%         set(icons(i),'LineWidth',10);
%     end
% end

% Now animate
marks = [];
for tIdx = 1:length(TimeCenters)
    t = TimeCenters(tIdx);
    cPhaseIdx = find(t < tEvents,1); % first idx is idx of current phase
    cPhase = Phases2{cPhaseIdx};
    title({tstr,sprintf('%s: %.2fs',cPhase,t)},StyleArgs{:});
    
    % Draw a marker for each trace at the current time step
%     tmp = plot(x(tIdx,:),y(tIdx,:),'.','MarkerSize',10);
    for i = 1:length(x)
        if tIdx == 1
            marks(i) = plot(x{i}(tIdx),y{i}(tIdx),'.','MarkerSize',40,'Color',clr(i,:));
        else
            set(marks(i),'XData',x{i}(tIdx));
            set(marks(i),'YData',y{i}(tIdx));
        end
    end
    MovieFrames(tIdx) = getframe(gcf);
    
    if ~isempty(Frames2Save) && any(tIdx == Frames2Save{2})
        matchIdx = find(tIdx == Frames2Save{2});
        plt.SaveFigure(1,SaveDir,[SaveFilePref '-' Frames2Save{1}{matchIdx}],'PNG','SVGI');
    end
end


end

function plotERAs(PlotData,TimeCenters,Phases,TimeWindows,PlotConds,SubplotIDs,pnl,StyleArgs,varargin)
[varargin,Colors] = Utilities.ProcVarargin(varargin,'Colors',[]);

% Plots the ERAs for the specified unit and conditions into the current plot
NumSubPlots = length(unique(SubplotIDs));
for cPlotIdx = 1:NumSubPlots % each subplot
    cPlotConds = PlotConds(SubplotIDs == cPlotIdx);
    cPlotData = PlotData(SubplotIDs == cPlotIdx);
    pnl(cPlotIdx,1).select();

    %% Plot
    % cPlotData{ii,jj} is Trials x TimeBins.
    Labels = getLegendLabels(cPlotConds,cPlotData);
    if ~isempty(Colors)
        Analyze.plotEventRelatedAverage(cPlotData,Labels,'TimeVec',TimeCenters,'NoLine',...
            'Legend','LegendLoc','NorthWest','LegendFontSize',12,'Colors',Colors);
    else
        Analyze.plotEventRelatedAverage(cPlotData,Labels,'TimeVec',TimeCenters,'NoLine',...
            'Legend','LegendLoc','NorthWest','LegendFontSize',12,'useBootStrap');
    end
%     Analyze.plotEventRelatedAverage(cPlotData,Labels,'TimeVec',TimeCenters,'NoLine',...
%         'Legend','LegendLoc','NorthWest','LegendFontSize',12,'Colors',jet(length(cPlotConds)),...
%         'std');
    axis tight;
    xlabel('Time (s)',StyleArgs{:});
    ylabel('Score',StyleArgs{:});
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
        text(tEvents(ii),ytext,Phases{ii},StyleArgs{:});
    end
end

end

function Labels = getLegendLabels(PlotConds,PlotData)
if isempty(PlotConds{1})
    Labels = {'All'};
    return;
end

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
    cLabel = sprintf('%s: %d',strcat(tempCond{:}),size(PlotData{jj},1));
    Labels = [Labels cLabel];
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