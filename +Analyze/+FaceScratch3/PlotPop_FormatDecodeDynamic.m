

%%

[varargin,rescale]   = Utilities.ProcVarargin(pltFCNArgs,'rescale',0);

% typeID={'DecodeTest_corrNM'};
%
% cvAccuracy=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,typeID{decodeType}));
%
% FigHandels(1)=plt.fig('units','inches','width',8,'height',8,'font','Arial','fontsize',12);
%
%  imagesc(mean(cvAccuracy,3))
%
% plt.vline([.5:5:20.5],{'k','linewidth',2}); plt.hline([.5:5:20.5],{'k','linewidth',2});
% colormap(jet)
% colorbar
%%

F=Analyze.returnFieldValues(PlotData,'DecodeTest_cvAccuracy');
FOut=[];
for sessionIDX=1:length(F);
    for format1=1:4
        for format2=1:4
            if sessionIDX==1
                FOut{format1,format2}=[];
            end
            FOut{format1,format2}=cat(3,FOut{format1,format2},F{sessionIDX}{format1,format2});
        end
    end
end
%%
for format1=1:4
    for format2=1:4
        
        FOutMu{format1,format2}=mean(FOut{format1,format2},3);
    end
end


%%
Conditions={'Fc','Fs','Oc','Os'};
plt.fig('units','inches','width',8,'height',6,'font','Arial','fontsize',12);
pnl = panel();  pnl.margin=10; pnl.pack(4,4);

allVals=cat(1,FOutMu{:}); mx=max(allVals(:)); mn=min(allVals(:));
pk95=prctile(allVals(:),[5 95]);



for format1=1:4
    for format2=1:4
        
        pnl(format1,format2).select();
        if rescale
            imagesc(FOutMu{format1,format2},[20 pk95(2)]);
        else
            allVals=FOutMu{format1,format2}; mx=max(allVals(:)); mn=min(allVals(:));
            pk95=prctile(allVals(:),[30 95]);
            imagesc(FOutMu{format1,format2},[20 pk95(2)]);
        end
        colormap(jet)
%         colorbar
        axis tight
        axis image
        set(gca, 'XTick',1:5:length(opts.timeWindow))
        set(gca, 'YTick',1:5:length(opts.timeWindow))
        
%         set(gca, 'XTickLabel',opts.timeWindow(1:5:end)-opts.winSize/2)
%         set(gca, 'YTickLabel',opts.timeWindow(1:5:end)-opts.winSize/2)
        
          set(gca, 'XTickLabel',opts.timeWindow(1:5:end))
        set(gca, 'YTickLabel',opts.timeWindow(1:5:end))
        
        xlabel(sprintf('From: %s ',Conditions{format1}))
        ylabel(sprintf('To: %s ',Conditions{format2}))
        
        plt.vline(10,{'y','linewidth',1})
        plt.hline(10,{'y','linewidth',1})
    end
end
