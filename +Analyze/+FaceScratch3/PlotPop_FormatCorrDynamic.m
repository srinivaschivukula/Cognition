

%%

[varargin,rescale]   = Utilities.ProcVarargin(pltFCNArgs,'rescale',0);
[varargin,Sessions]   = Utilities.ProcVarargin(pltFCNArgs,'Sessions',[]);

% typeID={'DecodeTest_corrNM'};
%
% cvAccuracy=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,typeID{decodeType}));
%
% FigHandels(1)=plt.fig('units','inches','width',8,'height',8,'font','Arial','fontsize',12);
%
%  imagesc(mean(cvAccuracy,3))
%
% plt.vline([.5:5:20.5],{'k','linewidth',2}); plt.hline([.5:5:20.5],{'k','linewidth',2});
% colormap(jet)s
% colorbar

tickLoc=6; StartInds=1;
%%

F=Analyze.returnFieldValues(PlotData,'DecodeTest_FormatCorrsPC');
FOut=[];
for sessionIDX=1:length(F)
    for format1=1:4
        for format2=1:4
            if sessionIDX==1
                FOut{format1,format2}=[];
            end
            FOut{format1,format2}=cat(3,FOut{format1,format2},F{sessionIDX}{format1,format2}(StartInds:end,StartInds:end,:));
        end
    end
end

timeWindow=opts.timeWindow(StartInds:end);
%%
for format1=1:4
    for format2=1:4
        if isempty(Sessions)
            FOutMu{format1,format2}=mean(FOut{format1,format2},3);
        else
            FOutMu{format1,format2}=mean(FOut{format1,format2}(:,:,Sessions),3);
        end
    end
end


% for format1=1:4
%     for format2=1:4
%       FOutMu{format1,format2}=mean(FOut{format1,format2}(:,:,8:13),3)-mean(FOut{format1,format2}(:,:,1:7),3);
%     end
% end
%%
figure; hold on

v=diag(FOutMu{1,1}); v=v-v(1);v=v/max(v);
plot(v,'g')
v=diag(FOutMu{2,2}); v=v-v(1);v=v/max(v);
plot(v,'k')
v=diag(FOutMu{3,3});v=v-v(1); v=v/max(v);
plot(v,'b')
v=diag(FOutMu{4,4});v=v-v(1); v=v/max(v);
plot(v,'r')
%%
ViewPoints={'Fc','Fs','Oc','Os'};
Conditions=ViewPoints;
plt.fig('units','inches','width',10,'height',8,'font','Arial','fontsize',12);
pnl = panel();  pnl.margin=20; pnl.pack(4,4);
StyleArgs = {'FontName','Arial','FontSize',10,'FontWeight','bold'};

allVals=cat(1,FOutMu{:}); mx=max(allVals(:)); mn=min(allVals(:));
pk95=prctile(allVals(:),[5 95]);



for format1=1:4
    for format2=1:format1
        
        pnl(format1,format2).select();
        if rescale
            imagesc(FOutMu{format1,format2},[0 pk95(2)]);
            %             imagesc(FOutMu{format1,format2},[pk95(1) pk95(2)]);
        else
            allVals=FOutMu{format1,format2}; mx=max(allVals(:)); mn=min(allVals(:));
            pk95=prctile(allVals(:),[5 95]);
            imagesc(FOutMu{format1,format2},[0 pk95(2)]);
        end
        %         colormap(flipud(parula)); colorbar
        colormap(parula); h = colorbar;
        
%         if format2==1
            yticks([1 ceil(length(timeWindow)/2) length(timeWindow)])
            set(gca, 'YTickLabel',timeWindow([1 ceil(length(timeWindow)/2) length(timeWindow)]))
%         end
%         if format1==4
            xticks([1 ceil(length(timeWindow)/2) length(timeWindow)])
            set(gca, 'XTickLabel',timeWindow([1 ceil(length(timeWindow)/2) length(timeWindow)]))
%         end
        
%         if format2~=1
            set(gca,'YTickLabel',[])
%         end
%         if format1~=4
            set(gca,'XTickLabel',[])
%         end
        
        xlabel(Conditions(format1))
        ylabel(Conditions(format2))
        
        plt.vline(tickLoc,{'-','linewidth',1.25,'Color',[240 240 240]./255})
        plt.hline(tickLoc,{'-','linewidth',1.25,'Color',[240 240 240]./255})
%         if format1==1 && format2==1
            ylims = ylim();xlims = xlim();
            ytext = ylims(2) + .05*diff(ylims);
            text(tickLoc,ytext+1,'Go',StyleArgs{:},'fontsize',11,'color',[226 28 26]./255);
            text(xlims(2)+.05*diff(ylims),tickLoc,'Go',StyleArgs{:},'fontsize',10,'color',[226 28 26]./255);
%             h.Visible='off';
%         end
%         if format2~=format1
%             h.Visible='off';
%         end
        
        axis tight
        ax1 = gca;
        set(ax1,'box','off',StyleArgs{:});
        ax1.XRuler.Axle.LineStyle = 'none';
        
        
    end
end

%%

