

%%
% Plot significance of the population through time with the specified
% threshold value (0.05 by default.) multiple plots allowed.



args=opts.PlotFCN{plotFCNidx,3};
[args,PValField]=ProcVarargin(args,'PVal',[]);
[args,SigThresh]=ProcVarargin(args,'SigThresh',0.05);
%%
[args,ylimVals]=ProcVarargin(args,'ylim',[]);
[args,Array]=ProcVarargin(args,'Array',[]);
[args,ecdfFlag]=ProcVarargin(args,'ecdf');
[args,LimData]=ProcVarargin(args,'LimData',[]);
[args,LimThresh]=ProcVarargin(args,'LimThresh',0.05);
[args,WinType]=ProcVarargin(args,'WinType','WindowStarts');
[args,barPlot]=ProcVarargin(args,'barPlot');

if isempty(PValField), error('the arguement PVal must be specified'); end
if ischar(PValField), PValField={PValField}; end

[args,clr]=ProcVarargin(args,'clr',[]);


if ischar(PValField)
    PValField={PValField};
end
%%

ecdfFlag=1;
% barPlot=1;
[args,pnl]=ProcVarargin(args,'pnl',[]);
if isempty(pnl)
    plt.fig('units','inches','width',12,'height',4,'font','Arial','fontsize',12);
    p = panel();
else
    [args,pnlLoc]=ProcVarargin(args,'pnlLoc',[]);
    if ~isempty(pnlLoc)
        if length(pnlLoc)==1
            p=pnl(pnlLoc);
        else
            p=pnl(pnlLoc(1),pnlLoc(2));
        end
    else
        p=pnl; % just get directly
    end
end

p.pack('h', hPack); p.marginright=5; p.marginleft=5;
p.fontsize=12;p.fontname='arial';
warning off
if ~isempty(LimData)
    AnyTuning=Analyze.returnFieldValues(PlotData{1},LimData);
end
if isfield(opts, 'AddArg')
    PlotData=Analyze.SubSelectTrials(PlotData,opts.AddArg{:});
end
for phaseIDX=1:length(Phases)
    ax(phaseIDX)=p(phaseIDX).select(); hold on
    cPlotData=PlotData{phaseIDX};
    
    
    for i=1:length(PValField)
        try
            PVal=Analyze.returnFieldValues(cPlotData,PValField{i});
%             [H,P,adj_p]=Utilities.MultipleComparisonsCorrection(PVal,'method','fdr')
%             PVal=adj_p;
            if ~isempty(Array)
            ArrayValues=cell2mat(Analyze.returnFieldValues(cPlotData,'Array'));
            PVal=PVal(Array==ArrayValues,:);
            end
        catch
            return
        end
        
        if all(isnan(PVal(:))), return; end
        [args,titleTXT]=ProcVarargin(args,'title',[Phases{phaseIDX} '-' PValField{i}]);
        if ~isempty(LimData)
            PVal=PVal(AnyTuning<LimThresh,:);
        end
        
        if ecdfFlag
            clear mu ci
            for j=1:size(PVal,2)
                [f,x,flo,fup] = ecdf(PVal(:,j));
                
                if ischar(SigThresh)
                    [~,effAlphaVal,adj_p]= Utilities.MultipleComparisonsCorrection(PVal(:,j),'method',SigThresh);
                    mu(j)=interp1(x(2:end),f(2:end),effAlphaVal);
                ci(2,j)=interp1(x(2:end),flo(2:end),effAlphaVal);
                ci(1,j)=interp1(x(2:end),fup(2:end),effAlphaVal);
                else
                
                mu(j)=interp1(x(2:end),f(2:end),SigThresh);
                ci(2,j)=interp1(x(2:end),flo(2:end),SigThresh);
                ci(1,j)=interp1(x(2:end),fup(2:end),SigThresh);
                end
            end
            %                    mu=smooth(mu,3);
            %                 ci(2,:)=smooth(ci(2,:),3);
            %                 ci(1,:)=smooth(ci(1,:),3);
            if barPlot
                
                %                  h=bar(Time{phaseIDX}.(WinType),100*mu,'BaseValue',SigThresh,'FaceColor',[.5 .5 .5]);
                %                  errorbar(Time{phaseIDX}.(WinType),100*mu,100*(ci(1,:)-mu(:)'),100*(mu(:)'-ci(2,:)),'.r','linewidth',1);
                %       errorbarxy(Time{phaseIDX}.(WinType), ...
                %                  100*mu,[],[100*(ci(1,:)-mu(:)');100*(mu(:)'-ci(2,:))],'r','linewidth',2);
                %%
                
                % dx=mean(diff(Time{phaseIDX}.(WinType)))/2;
                % h=stairs([Time{phaseIDX}.(WinType); Time{phaseIDX}.(WinType)(end)+dx*2]-dx,100*[mu mu(end)],'Color',[.5 .5 .5],'linewidth',2);
                %          plt.errorbarxy(Time{phaseIDX}.(WinType), ...
                %                  100*mu,[],[100*(ci(1,:)-mu(:)');100*(mu(:)'-ci(2,:))],'Color',[.5 .5 .5],'linewidth',2);
                %%
                plt.barpatch(Time{phaseIDX}.(WinType),{100*mu,100*ci},'patchbar',1,'barwidth',.075)
            else
                plt.shadedErrorBar(Time{phaseIDX}.(WinType),100*mu,100*ci,'CI');
            end
        else
        
            
            if ischar(SigThresh)
                clear PPerc
                  H= Utilities.MultipleComparisonsCorrection(PVal,'method',SigThresh);
                  mu=sum(H,1)/size(H,1);
                  ci=bootci(2000,@(x)sum(x,1)/size(x,1),H)
                  plt.shadedErrorBar(Time{phaseIDX}.(WinType),100*mu,100*ci,'CI');

                  
                for j=1:size(PVal,2)
                    H= Utilities.MultipleComparisonsCorrection(PVal(:,j),'method',SigThresh);
                    PPerc(j)=nnz(H)/length(H);
                end
                
            else
                PThresh=double(PVal<SigThresh)';
                PPerc=sum(PThresh,2)/size(PThresh,2);
            end
            
            
            if isempty(clr)
                area(Time{phaseIDX}.(WinType),PPerc*100)
            else
                area(Time{phaseIDX}.(WinType),PPerc*100,'FaceColor',clr{i})
            end
        end
    end
    
    
    if phaseIDX==1
        legend(PValField,'Location','best')
    end
    
    axis tight
    xlabel('Time (S)')
    ylabel(sprintf('% Signif'))
    titleTXT=strrep(titleTXT,'_','-');
     if ~isempty(Array)
             titleTXT=sprintf('%s: Array %d',titleTXT,Array);
     else
             titleTXT=sprintf('%s: Array All',titleTXT);

     end
     
    title([titleTXT sprintf(' (N=%d)',size(PVal,1))])
    if ~isempty(ylimVals)
        ylim(ylimVals)
    end
    if phaseIDX~=1
        set(gca,'yTick',[])
    end
    if phaseIDX==1
        plt.hline(SigThresh*100,{':','color',[0 0 0],'linewidth',2},'Chance',0);
    else
        plt.hline(SigThresh*100,{':','color',[0 0 0],'linewidth',2});
    end
end
xlabel('Time (S)')
ylabel('%')
yl=[ax.YLim];
linkaxes(ax,'y')
if isempty(ylimVals)
    ylim([0 max(yl)])
else
    ylim(ylimVals)
end
warning on

%%
