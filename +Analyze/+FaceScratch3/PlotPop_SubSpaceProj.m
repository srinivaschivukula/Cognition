
%      cmps={[1 2],[1 3],[2 3],[ 1 2 3],1,2,3};
%%
FigHandels=[];
MeasName='gdP';
% MeasName='gR2';

Labels={'Pinch'    'Press'    'Rub'    'Tap'};
condName='FcFsOc';
L={'FcFs','FcOc','FsOc','FcFsOc','Fc','Fs','Oc'};

% condName='FcFsOs';
% L={'FcFs','FcOs','FsOs','FcFsOs','Fc','Fs','Os'};
 
% condName='FcOcOs';
% L={'FcOc','FcOs','OcOs','FcOcOs','Fc','Oc','Os'};

% condName='FsOcOs';
% L={'FsOc','FsOs','OcOs','FsOcOs','Fs','Oc','Os'};

act2Avg=[1:4];

Meas=Analyze.returnFieldValues(PlotData,[condName '_' MeasName]);
clear MatchedVals MisMatchedVals
for idx=1:7
    for dateIDX=1:length(Meas)
        
        
        d=eye(4); d(find(d))=nan;
        cMeas=squeeze(Meas{dateIDX}(idx,:,:));
        
        
        
        if idx<5
            A=diag(cMeas);
            B=cMeas; B=B+d;
            B=(B+B')/2;
            B=nanmean(B); B=B(:);
            %     B./A
            MatchedVals(:,idx,dateIDX)=A;
            MisMatchedVals(:,idx,dateIDX)=B;
        else
            A=cMeas(:,1);
            
            MatchedVals(:,idx,dateIDX)=A;
        end
    end
end
%
FigHandels(1)=plt.fig('units','inches','width',16,'height',6,'font','Arial','fontsize',12);
pnl = panel();  pnl.margin=15; pnl.pack('h',7);
pnl(1).select()

for idx=1:7
    pnl(idx).select(); hold on
    
    if idx>=5
        M=squeeze(MatchedVals(:,idx,:));
        yy=.15*randn(size(M,1),size(M,2))+repmat((1:4)',1,size(M,2));
        plot(yy,M,'.','markersize',10,'color',[.5 .5 .5]);
        xlim([.5 5.5])
        
        Vals=M';
        try
            PR=bootci(5000,@mean,Vals);
        catch
            PR= prctile(Vals,[5 95]);
        end
        mu=mean(Vals);
                errorbar(1:4,mu,mu-PR(1,:),mu-PR(2,:),'.k','linewidth',2);
        
        
        if strcmp(MeasName,'AUC')
            ylim([.7 1])
        elseif strcmp(MeasName,'gdP')
            ylim([0 7])
        else
            ylim([-.1 .7])
        end
    else
        M=squeeze(MatchedVals(:,idx,:));
        yy=.15*randn(size(M,1),size(M,2))+repmat((1:4)',1,size(M,2));
        plot(yy,M,'.','markersize',10,'color',[0 1 0]);
        xlim([.5 5.5])
        
        Vals=M';
        
        
        
        M=squeeze(MisMatchedVals(:,idx,:));
        yy=.15*randn(size(M,1),size(M,2))+repmat((1:4)',1,size(M,2));
        plot(yy,M,'.','markersize',10,'color',[1 0 0]);
        xlim([.5 5.5])
        
        try
            PR=bootci(5000,{@mean, Vals}, 'type','per');
        catch
            PR= prctile(Vals,[5 95]);
        end
        mu=mean(Vals);
                errorbar((1:4)-0.05,mu,mu-PR(1,:),mu-PR(2,:),'.g','linewidth',2);
        
        Vals=M';
        try
            PR=bootci(5000,{@mean, Vals}, 'type','per');
        catch
            PR= prctile(Vals,[5 95]);
        end
        mu=mean(Vals);
                errorbar((1:4)+0.05,mu,mu-PR(1,:),mu-PR(2,:),'.r','linewidth',2);
        
        if strcmp(MeasName,'AUC')
            ylim([.7 1])
        elseif strcmp(MeasName,'gdP')
            ylim([0 7])
        else
            ylim([-.1 .7])
        end
    end
    set(gca,'XTickLabel',Labels); xticks([1:4])
    set(gca,'XTickLabelRotation',30)
    title(sprintf('AUC : %s subspace',L{idx}))
end


%
%%
FigHandels(2)=plt.fig('units','inches','width',6,'height',5,'font','Arial','fontsize',12);


MatchedValsA=permute(MatchedVals,[1 3 2]);MatchedValsA=MatchedValsA(act2Avg,:,:);
MatchedValsA=reshape(MatchedValsA,size(MatchedValsA,1)*size(MatchedValsA,2),size(MatchedValsA,3));

MisMatchedValsA=permute(MisMatchedVals,[1 3 2]);MisMatchedValsA=MisMatchedValsA(act2Avg,:,:);
MisMatchedValsA=reshape(MisMatchedValsA,size(MisMatchedValsA,1)*size(MisMatchedValsA,2),size(MisMatchedValsA,3));


% L={'TL0','TF','FL0','TFL0','T','F','L0','','','',''};
%
% boxplot([MatchedValsA MisMatchedValsA],'notch','on','Labels',L)
%



MM=[];
for i=1:4
    MM=[MM MatchedValsA(:,i)  MisMatchedValsA(:,i)];
end
MM=[MM MatchedValsA(:,5:7)];


L={'TL0','','TF','','FL0','','TFL0','','T','F','L0'};
MM=MM(:,[7:8 1:6 ]);
L=L([7:8 1:6]);

ci=bootci(2000,{@mean,MM},'type','per');
mu=mean(MM);
plt.barpatch([1 1 3 3 5 5 7 7]+[-.25 .25 -.25 .25 -.25 .25 -.25 .25],{mu,ci},'barcolor',[0 0 0],'patchbar',1,'barwidth',.3,'barname',L,'fontsize',20,'xl','Format','yl','% Tuned (fdr)','t','Action tuning across formats')


% boxplot(MM,'notch','on','Labels',L,'PlotStyle','compact')

if strcmp(MeasName,'AUC')
    ylim([.7 1.01])
else
    ylim([-.05 .41])
end
plt.hline(0,{'--','color',[.5 .5 .5]})

title('')

%%
%
% plt.fig('units','inches','width',6,'height',3,'font','Arial','fontsize',12);
% 
% 
% MatchedValsA=permute(MatchedVals,[1 3 2]);MatchedValsA=squeeze(mean(MatchedValsA([1 3 5],:,:),1));
% MisMatchedValsA=permute(MisMatchedVals,[1 3 2]);MisMatchedValsA=squeeze(mean(MisMatchedValsA([1 3 5],:,:),1));
% 
% hold on
% plot(MatchedValsA(:,2),'g')
% plot(MisMatchedValsA(:,2),'r')

