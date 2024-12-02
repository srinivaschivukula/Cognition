

%%
% Plot significance of the population through time with the specified
% threshold value (0.05 by default.) multiple plots allowed.



args=opts.PlotFCN{plotFCNidx,3};
[args,PlotValueField]=ProcVarargin(args,'PlotValue',[]);
[args,PValField]=ProcVarargin(args,'PVal',[]);
[args,SigThresh]=ProcVarargin(args,'SigThresh',0.05);

[args,ylimVals]=ProcVarargin(args,'ylim',[]);
[args,ecdfFlag]=ProcVarargin(args,'ecdf');
[args,LimData]=ProcVarargin(args,'LimData',[]);
[args,LimThresh]=ProcVarargin(args,'LimThresh',0.05);
[args,WinType]=ProcVarargin(args,'WinType','WindowStarts');
[args,barPlot]=ProcVarargin(args,'barPlot');

% if isempty(PValField), error('the arguement PVal must be specified'); end
if ischar(PValField), PValField={PValField}; end

[args,clr]=ProcVarargin(args,'clr',[]);

barPlot=0;
%%
plt.fig('units','inches','width',6,'height',4,'font','Arial','fontsize',12);
p = panel(); p.pack('h', hPack); p.margin=10;
warning off
if ~isempty(LimData)
    AnyTuning=Analyze.returnFieldValues(PlotData{1},LimData);
end

dx=mean(diff(Time{1}.(WinType)));
spacing=dx/length(PlotValueField)/2;
shifts=(0:(length(PlotValueField)-1))*spacing;
shifts=shifts-mean(shifts);
%%
for phaseIDX=1:length(Phases)
    ax(phaseIDX)=p(phaseIDX).select(); hold on
    cPlotData=PlotData{phaseIDX};
    
    %%
    cla; hold on; clr=lines(20);
    for i=1:length(PlotValueField)
        PlotValue=Analyze.returnFieldValues(cPlotData,PlotValueField{i});
        
        [args,titleTXT]=ProcVarargin(args,'title',[Phases{phaseIDX} '-' PlotValueField{i}]);
        if ~isempty(LimData)
            PlotValue=PlotValue(AnyTuning<LimThresh,:);
        end
        
        if ~isempty(PValField)
            PVal=Analyze.returnFieldValues(cPlotData,PValField{i});
            PValthresh=PVal<SigThresh;
            PlotValue(~PValthresh)=nan;
        end
        
        mu=nanmean(PlotValue,1);
        ci=bootci(2000,@nanmean,PlotValue);
        
        if barPlot
            errorbar(Time{phaseIDX}.(WinType)+shifts(i),mu,ci(1,:)-mu,mu-ci(2,:),'.','linewidth',2);
%              plt.errorbarxy(Time{phaseIDX}.(WinType)+shifts(i),mu,[],abs(bsxfun(@minus,mu,ci)),'linewidth',2,'color',clr(i,:));
        else
            plt.shadedErrorBar(Time{phaseIDX}.(WinType)+shifts(i),mu,ci,'lineProps',{'color',clr(i,:)},'CI');
        end
        
    end
    %%
    
    if phaseIDX==1
        legend(PlotValueField,'Location','best')
    end
    
    axis tight
    xlabel('Time (S)')
    ylabel(sprintf('% Signif'))
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
