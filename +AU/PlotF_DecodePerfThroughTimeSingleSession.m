%% TODO: turn this into a function, rather than script

%%
%Plot the basic elements of a linear model that accounts for neural firing
%as a function of the subject's choice and the cued action.  From this
%dataset it appears that there is no significant tuning prior to the action
%cue suggesting some type of heirarchical neural architecture where the
%high-level decision was performed elsewhere.
args=opts.PlotFCN{plotFCNidx,3};
[args,holdPanel]=ProcVarargin(args,'holdPanel');

if holdPanel
    global pnl
    if isempty(pnl)
        plt.fig('units','inches','width',5,'height',4,'font','Arial','fontsize',12);
        pnl = panel();  pnl.pack('h', hPack); pnl.margin=10;
    end
else
    clear global pnl
        plt.fig('units','inches','width',5,'height',4,'font','Arial','fontsize',12);
    pnl = panel();  pnl.pack('h', hPack); pnl.margin=10;
end


[args,ACCField]=ProcVarargin(args,'ACC',[]);
[args,ShuffleRankField]=ProcVarargin(args,'ShuffleRank',[]);
[args,SigThresh]=ProcVarargin(args,'SigThresh',0.05);
[args,nullDistField]=ProcVarargin(args,'nullDist',[]);
[args,nullDistArgs]=ProcVarargin(args,'nullDistArgs',{'NoEdges','NoLine'});
[args,CIfield]=ProcVarargin(args,'CI',[]);
[args,titleTXT]=ProcVarargin(args,'title','');
[args,ylimVals]=ProcVarargin(args,'ylim',[]);
[args,xlimVals]=ProcVarargin(args,'xlim',[]);
[args,winType]=ProcVarargin(args,'winType','WindowStarts');
[args,plotAverage]=ProcVarargin(args,'plotAverage',false);

% nullDistField=[];

warning off
% Phases=cellfun(@(x)x.Phase,PlotData);
for phaseIDX=1:length(Phases)
    try
    ax(phaseIDX)=pnl(phaseIDX).select(); axes(ax(phaseIDX));hold on
    cPlotData=PlotData{phaseIDX};
    ACC=Analyze.returnFieldValues(cPlotData,ACCField);
    %     ACC=smooth(ACC,2);
    
    if ~isempty(nullDistField)
        nullDist=Analyze.returnFieldValues(cPlotData,nullDistField);
        nullDist=cellfun(@transpose,nullDist,'UniformOutput' , false);
        nullDist=[nullDist{:}];
        mu=mean(nullDist,1);
        nullCI=prctile(nullDist,[5 90]);
        if ischar(nullDistArgs); nullDistArgs={nullDistArgs}; end
        plt.shadedErrorBar(Time{phaseIDX}.(winType),mu,nullCI,'CI',nullDistArgs{:})
    end
%     if ~isempty(CIfield)
%         CI=Analyze.returnFieldValues(cPlotData,CIfield);
%         CI=cellfun(@transpose,CI,'UniformOutput' , false);
%         CI=[CI{:}];
% %         errorbar(Time{phaseIDX}.(winType),ACC,ACC-CI(1,:),CI(1,:)-ACC,'k.')
%                 plt.shadedErrorBar(Time{phaseIDX}.(winType),ACC,CI,'CI',nullDistArgs,'transparent',.5)
%     end
    
    
    % Arrays with elements => array rows
    if iscell(ACC)
        % Each row is a separate element
        ACC = cell2mat(reshape(ACC, numel(ACC), 1));
    end
    
    % TODO: remove this NaN-filling
    % heuristic is that fully 0-accuracy rows must have been the result of
    % matrix zero-filling from trial dates without any of that type of
    % accuracy
    ACC_zero_rows = ~any(ACC,2);
    ACC(ACC_zero_rows,:) = NaN;
    
    % Plot the average
    if plotAverage
        Analyze.plotEventRelatedAverage({ACC.'},[ACCField],'useBootStrap','TimeVec',Time{phaseIDX}.(winType))
    else
        plot(Time{phaseIDX}.(winType),ACC,'o-','markersize',7,'markerfacecolor',[1 1 1],'linewidth',2)
    end
    
    % Draw the shuffle-test significance threshold
    SigThresh=.1;
    if ~isempty(ShuffleRankField)
        ShuffleRank=Analyze.returnFieldValues(cPlotData,ShuffleRankField);
        isSig=ShuffleRank<SigThresh;
        plot(Time{phaseIDX}.(winType),sum(isSig,1)/size(isSig,1),'-o','markersize',6,'markerfacecolor',[1 0 0])
    end    
    
    legend({'Accuracy'},'Location','best')
    axis tight
    ylim([0 100]);
    xlabel('Time (s)')
    ylabel('Decode Accuracy (%)')
    titleTXT=strrep(titleTXT,'_','-');
    title([titleTXT ' ' Phases{phaseIDX}])
    plt.hline(100/6,{'--','color',[.5 .5 .5]});
    catch err
        warning(err.identifier);
    end
end
yl=[ax.YLim];
linkaxes(ax,'y')
if isempty(ylimVals)
    ylim([min(yl) max(yl)])
else
    ylim(ylimVals)
end

if ~isempty(ylimVals)
    xlim(xlimVals)
end
warning on

%%
