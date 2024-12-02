

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

% nullDistField=[];

warning off
% Phases=cellfun(@(x)x.Phase,PlotData);
for phaseIDX=1:length(Phases)
    
    ax(phaseIDX)=pnl(phaseIDX).select(); axes(ax(phaseIDX));hold on
    cPlotData=PlotData{phaseIDX};
    ACC=Analyze.returnFieldValues(cPlotData,ACCField);
    %     ACC=smooth(ACC,2);
    
    try
    if ~isempty(nullDistField)
        nullDist=Analyze.returnFieldValues(cPlotData,nullDistField);
        nullDist=cellfun(@transpose,nullDist,'UniformOutput' , false);
        nullDist=[nullDist{:}];
        mu=mean(nullDist,1);
        nullCI=prctile(nullDist,[5 90]);
        if ischar(nullDistArgs); nullDistArgs={nullDistArgs}; end
        plt.shadedErrorBar(Time{phaseIDX}.(winType),mu,nullCI,'CI',nullDistArgs{:})
    end
    catch
        lasterr
       warning('problem with null dist') 
    end
    
    try
    if ~isempty(CIfield)
        CI=Analyze.returnFieldValues(cPlotData,CIfield);
        CI=cellfun(@transpose,CI,'UniformOutput' , false);
        CI=[CI{:}];
%         errorbar(Time{phaseIDX}.(winType),ACC,ACC-CI(1,:),CI(1,:)-ACC,'k.')
                plt.shadedErrorBar(Time{phaseIDX}.(winType)',ACC,CI,'CI',nullDistArgs,'transparent',.5)
    end
     catch
        lasterr
       warning('problem with null dist') 
    end
    
    %       ACC=smooth(ACC,3);
    plot(Time{phaseIDX}.(winType),ACC,'o-','markersize',7,'markerfacecolor',[1 1 1],'linewidth',2)
    SigThresh=.1;
    if ~isempty(ShuffleRankField)
        ShuffleRank=Analyze.returnFieldValues(cPlotData,ShuffleRankField);
        isSig=ShuffleRank<SigThresh;
        plot(Time{phaseIDX}.(winType)(isSig),ACC(isSig),'o','markersize',6,'markerfacecolor',[1 0 0])
    end
    %       plot(Time{phaseIDX}.(winType),Time{phaseIDX}.WindowMid*0+50,'--','color',[.7 .7 .7],'linewidth',2)
    %       plot(Time{phaseIDX}.(winType),mean(nullCI),'--','color',[.7 .7 .7],'linewidth',2)
    
    
    
    if phaseIDX==1
        tmp=xlim;
        xlim([tmp(1) 4])
        legend({'Null','Acu','Sig'},'Location','best')
    end
    axis tight
    xlabel('Time (S)')
    ylabel(sprintf('Decode Accuracy'))
    titleTXT=strrep(titleTXT,'_','-');
    title([Phases{phaseIDX} ' ' titleTXT])
    plt.hline(33.3,{'--','color',[.5 .5 .5]})
   
end
%%
yl=[ax.YLim];
linkaxes(ax,'y')
if isempty(ylimVals)
    ylim([min(yl) max(yl)])
else
    ylim(ylimVals)
end

if ~isempty(xlimVals)
    xlim(xlimVals)
end
warning on

%%
