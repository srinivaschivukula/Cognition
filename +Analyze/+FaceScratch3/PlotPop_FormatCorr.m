

%%
decodeType=1;


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

% a=squeeze(cvAccuracy(1,2,:));
% FigHandels(1)=
plt.fig('units','inches','width',5,'height',6,'font','Arial','fontsize',12);
pnl = panel();  pnl.margin=15; pnl.pack('v',2);


typeID={'corr','corrPart'};

for decodeType=1:length(typeID)
    pnl(decodeType).select()
    cvAccuracy=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,typeID{decodeType}));
    
    cvAccuracyShuff=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'corrShuff'));
    
    
    %     cvAccuracy=cell2mat(cvAccuracy);
    
    T={'Fc-Fs','Fc-Oc', 'Fc-Os','Fs-Oc', 'Fs-Os','Oc-Os'};
    IDX=[2 3 4 7 8 12];
    

    clear Vals
    for i=1:length(IDX);
        [I,J] = ind2sub([4 4],IDX(i));
        Vals(:,i)= squeeze(cvAccuracy(I,J,:));
        ValsShuff(:,i)= squeeze(cvAccuracyShuff(I,J,:));
    end

    hold on
    
    xlim([0 length(IDX)+1])
    %     plt.hline([20],{'--','color',[.5 .5 .5]},'Chance')
    % boxplot(Vals,'Labels',T,'plotstyle','traditional','Symbol','')
    
    
    % try
    PR=bootci(5000,{@mean,Vals},'type','per');
    %     PR=bootci(5000,@median,Vals);
    % catch
    %    PR= prctile(Vals,[5 95]);
    % end
    
    t=[.5 .5 .5 1 .5 1];
    t=cumsum(t);
    for i=1:size(Vals,2),
        [~,p(i)]=ttest(Vals(:,i),ValsShuff(:,i),'tail' , 'right');
        %     ptxt={}
                    text(t(i),0,sprintf('%0.3f',p(i)),'rotation',45)
    end
    %
    mu=mean(Vals); muOLD{decodeType}=mu;
    % errorbar(1:length(IDX),mu,mu-PR(1,:),mu-PR(2,:),'ok','linewidth',2);
    plot(repmat(linspace(-.1,.1,size(Vals,1))',1,size(Vals,2))+repmat(t,size(Vals,1),1),Vals,'r.','markersize',10,'color',[0 0 1]);
    
    %     plot(.1*randn(size(Vals,1),size(Vals,2))+repmat(1:length(IDX),size(Vals,1),1),Vals,'r.','markersize',10,'color',[0 0 1])
    plt.barpatch(t,{mu,PR},'barcolor',[0 0 1],'patchbar',1,'barwidth',.3,'barname',T(1:length(mu)),'fontsize',20,'xl','Format','yl','% Tuned (fdr)','t','Action tuning across formats')
    %
    if decodeType==1
       mu=mean(ValsShuff);
    % errorbar(1:length(IDX),mu,mu-PR(1,:),mu-PR(2,:),'ok','linewidth',2);
    plot(repmat(linspace(-.1,.1,size(ValsShuff,1))',1,size(ValsShuff,2))+repmat(t,size(ValsShuff,1),1),ValsShuff,'r.','markersize',10,'color',[1 0 0]);
    PR=bootci(5000,{@mean,ValsShuff},'type','per');
    %     plot(.1*randn(size(Vals,1),size(Vals,2))+repmat(1:length(IDX),size(Vals,1),1),Vals,'r.','markersize',10,'color',[0 0 1])
    plt.barpatch(t,{mu,PR},'barcolor',[1 0 0],'patchbar',1,'barwidth',.3,'barname',T(1:length(mu)),'fontsize',20,'xl','Format','yl','% Tuned (fdr)','t','Action tuning across formats') 
    end
    
    if decodeType==2
        plot(t,muOLD{1},'rd','MarkerFaceColor',[1 0 0]);
    end
    %     plot([8 8 8 8], mu(end-3:end),'k>','markerfacecolor',[.5 .5 .5]);
    %     text([8 8 8 8]+.1, mu(end-3:end),T(end-3:end),'fontweight','bold','fontsize',12);
    plt.hline(0,{'k--'},'Chance',[1 0])
    set(gca,'XTick',t)
    % set(gca,'XTickLabel',T)
    set(gca,'XTickLabelRotation',20)
    title('PopCorr')
    ylabel('Corr')
    %     xlabel('Raw Firing')
    ylim([-.2 1])
    set(gca,'YTick',0:.2:.8)
end

%%