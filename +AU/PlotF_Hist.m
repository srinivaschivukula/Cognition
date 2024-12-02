

%%
%Plot the basic elements of a linear model that accounts for neural firing
%as a function of the subject's choice and the cued action.  From this
%dataset it appears that there is no significant tuning prior to the action
%cue suggesting some type of heirarchical neural architecture where the
%high-level decision was performed elsewhere.

plt.fig('units','inches','width',4,'height',3,'font','Arial','fontsize',12);

p = panel(); p.pack('h', hPack); p.margin=10;

args=opts.PlotFCN{plotFCNidx,3};
[args,PlotDataField]=ProcVarargin(args,'PlotData',[]);
[args,titleTXT]=ProcVarargin(args,'title','');

 


warning off
% Phases=cellfun(@(x)x.Phase,PlotData);
for phaseIDX=1:length(Phases)
    ax(phaseIDX)=p(phaseIDX).select(); hold on
    cPlotData=PlotData{phaseIDX};
    ACC=Analyze.returnFieldValues(cPlotData,PlotDataField);
    
    histogram(ACC)
    
    
    if phaseIDX==1
        legend(PlotDataField,'Location','best')
    end
    axis tight
    xlabel('Values')
    ylabel('Frequency')
    titleTXT=strrep(titleTXT,'_','-');
    title(titleTXT)
end
yl=[ax.YLim];
linkaxes(ax,'y')
if isempty(ylimVals)
    ylim([min(yl) max(yl)])
else
    ylim(ylimVals)
end
warning on

%%
