

%      cmps={[1 2],[1 3],[2 3],[ 1 2 3],1,2,3};
%%
% [~,Type]   = Utilities.ProcVarargin(pltFCNArgs,'Type','T');


% ConditinList={{'T','F','L0'},...
%     {'F','L0','L1'}}

FigHandels=[];

% if Type=='T'
%     condName='TFL0';
% L={'TFL0','','FL0','','TF','','TL0',''};
% else
% condName='FL0L1';
% L={'FL0L1','','L0L1','','FL0','','FL1',''};
% end

condName='FcFsOc';
L={'FcFsOc','','FsOc','','FcFs','','FcOc',''};

MeasName='gdP';
MeasName='PairComp';
Labels={'Pinch'    'Press'    'Rub'    'Tap'};
cmps={[1 2],[1 3],[2 3],[ 1 2 3],1,2,3};

%  L={'TF','','TL0','','FL0','','TFL0',''};
%  L={'FL0','','FL1','','L0L1','','FL0L1',''};

% cmps={[1 2 0],[1 3 0],[2 3 0],[ 1 2 3],[1 0 0],[2 0 0],[3 0 0]};



cmps={[ 1 2 3],[2 3 0],[1 2 0],[1 3 0],[1 0 0],[2 0 0],[3 0 0]};



Type =  'WithinMin';;%,WithinMean 'WithinMin','AcrossMean','AcrossMin'
comboType='mean';
act2Avg=[1:4];
% if strcmp(opts.PlotSubject,'P3')
% act2Avg=[1 3 4];
% end
RemoveDailyVariance=1;

PairComp=Analyze.returnFieldValues(PlotData,[condName '_PairComp' ]);
PairCompP=Analyze.returnFieldValues(PlotData,[condName '_PairCompP' ]);
RefCond=Analyze.returnFieldValues(PlotData,[condName '_RefCond']);

%
% foo=[8:11]
% PairComp=PairComp(foo);
% PairCompP=PairCompP(foo);
% RefCond=RefCond(foo);

% id2rm=[2 4]
% for i=1:length(PairComp)
%     for j=[(id2rm) id2rm+5 id2rm+10];
%    PairComp{i}(:,:,j)=nan;
%     end
% end



FormMap{1}=[1:4];
FormMap{2}=[5:8];
FormMap{3}=[9:12];
%
for dateIDX=1:length(PairComp)
    cPairComp=PairComp{dateIDX};
    cPairCompP=PairCompP{dateIDX};
    cRefCond=RefCond{dateIDX};
    
    for CMPidx=1:4
        cCMP=cmps{CMPidx};
        ids=cRefCond(:,1,1)==cCMP(1) & cRefCond(:,2,1)==cCMP(2) & cRefCond(:,3,1)==cCMP(3);
        
        redPairCmp=cPairComp(ids,:,:);
        redPairCompP=cPairCompP(ids,:,:);
        redRefCond=cRefCond(ids,:,:);
        clear Mtmp nMtmp
        for actIDX=1:4
            
            matchIDX= squeeze(redRefCond(:,1,2))==actIDX &  squeeze(redRefCond(:,2,2)) ==actIDX;
            mismatchIDX= squeeze(redRefCond(:,1,2))==actIDX &  squeeze(redRefCond(:,2,2)) ~=actIDX;
            %%
            switch CMPidx
                case {2 3 4}
                    
                    F1=cmps{CMPidx}(1);
                    F2=cmps{CMPidx}(2);
                    F1pc= squeeze(redPairCmp(matchIDX,1,:));
                    F2pc= squeeze(redPairCmp(matchIDX,2,:));
                    
                    Mtmp(actIDX,1)=ComputeMeasureMod(Type,F1pc,F1);
                    Mtmp(actIDX,2)=ComputeMeasureMod(Type,F2pc,F2);
                    
                    nF1pc= squeeze(redPairCmp(mismatchIDX,1,:))';
                    nF2pc= squeeze(redPairCmp(mismatchIDX,2,:))';
                    
                    nMtmp(actIDX,1)=ComputeMeasureMod(Type,nF1pc,F1);
                    nMtmp(actIDX,2)=ComputeMeasureMod(Type,nF2pc,F2);
                    
                case 1
                    
                    F1=cmps{CMPidx}(1);
                    F2=cmps{CMPidx}(2);
                    F3=cmps{CMPidx}(3);
                    F1pc= squeeze(redPairCmp(matchIDX,1,:));
                    F2pc= squeeze(redPairCmp(matchIDX,2,:));
                    F3pc= squeeze(redPairCmp(matchIDX,3,:));
                    
                    Mtmp(actIDX,1)=ComputeMeasureMod(Type,F1pc,F1);
                    Mtmp(actIDX,2)=ComputeMeasureMod(Type,F2pc,F2);
                    Mtmp(actIDX,3)=ComputeMeasureMod(Type,F3pc,F3);
                    
                    nF1pc= squeeze(redPairCmp(mismatchIDX,1,:))';
                    nF2pc= squeeze(redPairCmp(mismatchIDX,2,:))';
                    nF3pc= squeeze(redPairCmp(mismatchIDX,3,:))';
                    
                    nMtmp(actIDX,1)=ComputeMeasureMod(Type,nF1pc,F1);
                    nMtmp(actIDX,2)=ComputeMeasureMod(Type,nF2pc,F2);
                    nMtmp(actIDX,3)=ComputeMeasureMod(Type,nF3pc,F3);
                otherwise
                    disp('No Resp')
            end
        end
        
        
        switch comboType
            case 'mean'
                MM=mean(Mtmp,2)';
                nMM=mean(nMtmp,2)';
            case 'min'
                MM=min(Mtmp,[],2)';
                nMM=min(nMtmp,[],2)';
        end
        M1tot{CMPidx}(dateIDX,:)=MM;
        nM1tot{CMPidx}(dateIDX,:)=nMM;
    end
end

MatchedVals=M1tot;
MisMatchedVals=nM1tot;

% %%
%
% for kkk=1:4
%     figure;
% for jjj=1:5
%  subplot(5,1,jjj)
%     plot(MatchedVals{kkk}(:,jjj)-MisMatchedVals{kkk}(:,jjj))
%     axis tight
%     ylim([-1 1])
%     plt.hline(0)
% end
% end
% %
%%
%%
%     %%
%     clear MatchedVals MisMatchedVals
%     for CMPidx=1:7;
%         for dateIDX=1:length(Meas)
%
%
%             d=eye(5); d(find(d))=nan;
%             cMeas=squeeze(Meas{dateIDX}(CMPidx,:,:));
%
%
%
%             if CMPidx<5
%                 A=diag(cMeas);
%                 B=cMeas; B=B+d;
%                 B=(B+B')/2;
%                 B=nanmean(B); B=B(:);
%                 %     B./A
%                 MatchedVals(:,CMPidx,dateIDX)=A;
%                 MisMatchedVals(:,CMPidx,dateIDX)=B;
%             else
%                 A=cMeas(:,1);
%
%                 MatchedVals(:,CMPidx,dateIDX)=A;
%             end
%         end
%     end
%     %
FigHandels(1)=plt.fig('units','inches','width',16,'height',6,'font','Arial','fontsize',12);
pnl = panel();  pnl.margin=15; pnl.pack('h',7);
pnl(1).select()

for CMPidx=1:4
    pnl(CMPidx).select(); hold on
    
    if CMPidx>=5
        M=MatchedVals{CMPidx}';
        yy=.15*randn(size(M,1),size(M,2))+repmat((1:5)',1,size(M,2));
        plot(yy,M,'.','markersize',10,'color',[.5 .5 .5]);
        xlim([.5 5.5])
        
        Vals=M';
        try
            PR=bootci(5000,{@mean,Vals},'type','per');
        catch
            PR= prctile(Vals,[5 95]);
        end
        mu=mean(Vals);
        errorbar(1:4,mu,mu-PR(1,:),mu-PR(2,:),'.k','linewidth',2);
        
        
        if strcmp(MeasName,'AUC')
            ylim([.7 1])
        else
            ylim([-.1 .65])
        end
    else
        M=MatchedVals{CMPidx}';
        
        yy=.15*randn(size(M,1),size(M,2))+repmat((1:4)',1,size(M,2));
        plot(yy,M,'.','markersize',10,'color',[0 1 0]);
        xlim([.5 5.5])
        
        Vals=M';
        
        
        
        M=MisMatchedVals{CMPidx}';
        yy=.15*randn(size(M,1),size(M,2))+repmat((1:4)',1,size(M,2));
        plot(yy,M,'.','markersize',10,'color',[1 0 0]);
        xlim([.5 5.5])
        
        try
            PR=bootci(5000,{@mean,Vals},'type','per');
        catch
            PR= prctile(Vals,[5 95]);
        end
        mu=mean(Vals);
        errorbar((1:4)-0.05,mu,mu-PR(1,:),mu-PR(2,:),'.g','linewidth',2);
        
        Vals=M';
        try
            PR=bootci(5000,{@mean,Vals},'type','per');
        catch
            PR= prctile(Vals,[5 95]);
        end
        mu=mean(Vals);
        errorbar((1:4)+0.05,mu,mu-PR(1,:),mu-PR(2,:),'.r','linewidth',2);
        
        %             if strcmp(MeasName,'AUC')
        %                 ylim([.7 1])
        %             elseif strcmp(MeasName,'PairComp')
        %                 ylim([.7 1])
        %             else
        %                 ylim([-.1 .7])
        %             end
    end
    xticks([1:4])
    set(gca,'XTickLabel',Labels)
    set(gca,'XTickLabelRotation',30)
    title(sprintf('AUC : %s subspace',L{CMPidx}))
end



%     %%
FigHandels(2)=plt.fig('units','inches','width',6,'height',5,'font','Arial','fontsize',12);

MatchedValsA=cellfun(@(x)x(:,act2Avg),MatchedVals,'UniformOutput' , false);
MatchedValsA=cellfun(@(x)x(:),MatchedValsA,'UniformOutput' , false);
MatchedValsA=cat(2,MatchedValsA{:});


MisMatchedValsA=cellfun(@(x)x(:,act2Avg),MisMatchedVals,'UniformOutput' , false);
MisMatchedValsA=cellfun(@(x)x(:),MisMatchedValsA,'UniformOutput' , false);
MisMatchedValsA=cat(2,MisMatchedValsA{:});



MM=[];
for i=1:4
    MM=[MM MatchedValsA(:,i)  MisMatchedValsA(:,i)];
end
%     MM=[MM MatchedValsA(:,5:7)];

ci=bootci(2000,{@mean,MM},'type','per');
mu=mean(MM);
plt.barpatch([1 1 3 3 5 5 7 7]+[-.25 .25 -.25 .25 -.25 .25 -.25 .25],{mu,ci},'barcolor',[0 0 0],'patchbar',1,'barwidth',.3,'barname',L,'fontsize',20,'xl','Format','yl','% Tuned (fdr)','t','Action tuning across formats')
%%
%     L={'TL0','','TF','','FL0','','TFL0','','T','F','L0'};
%     MM=MM(:,[7:8 1:6 ]);
%     L=L([7:8 1:6]);
%
%     ci=bootci(2000,@mean,MM);
%     mu=mean(MM);
%     plt.barpatch([1 1 3 3 5 5 7 7]+[-.25 .25 -.25 .25 -.25 .25 -.25 .25],{mu,ci},'barcolor',[0 0 0],'patchbar',1,'barwidth',.3,'barname',L,'fontsize',20,'xl','Format','yl','% Tuned (fdr)','t','Action tuning across formats')
%
%
%     % boxplot(MM,'notch','on','Labels',L,'PlotStyle','compact')
%
%     if strcmp(MeasName,'AUC')
%         ylim([.7 1.01])
%     else
%         ylim([-.05 .41])
%     end
%     plt.hline(0,{'--','color',[.5 .5 .5]})
%
%     title('')
%
% %     %%
% %     %
% %     % plt.fig('units','inches','width',6,'height',3,'font','Arial','fontsize',12);
% %     %
% %     %
% %     % MatchedValsA=permute(MatchedVals,[1 3 2]);MatchedValsA=squeeze(mean(MatchedValsA([1 3 5],:,:),1));
% %     % MisMatchedValsA=permute(MisMatchedVals,[1 3 2]);MisMatchedValsA=squeeze(mean(MisMatchedValsA([1 3 5],:,:),1));
% %     %
% %     % hold on
% %     % plot(MatchedValsA(:,2),'g')
% %     % plot(MisMatchedValsA(:,2),'r')
