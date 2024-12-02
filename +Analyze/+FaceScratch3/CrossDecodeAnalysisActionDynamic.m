function CrossDecodeAnalysisActionDynamic(AllTrialData,Labels,TimeWindow,OutDir,Tag,varargin)
% Train on dim 1 to classify dim 2, test on how well it generalizes.
% e.g. train on nancy to diff bp, test how well bp diff'd in Tyson.

[~,PThresh]  = Utilities.ProcVarargin(varargin,'PThresh','fdr');
% [~,NIter] = Utilities.ProcVarargin(varargin,'NReps',100);
[~,Phase] = Utilities.ProcVarargin(varargin,'Phase','Go');
% [~,BasePhase] = Utilities.ProcVarargin(varargin,'BasePhase','Delay');
% [~,BaseWindow] = Utilities.ProcVarargin(varargin,'BaseWindow',[1.5 2.5]);
% [~,Dims] = Utilities.ProcVarargin(varargin,'SplitDims',{});
[~,DimNames] = Utilities.ProcVarargin(varargin,'DimNames',{});
% [~,SplitValNames] = Utilities.ProcVarargin(varargin,'SplitValNames',{});
[~,PBPSpec] = Utilities.ProcVarargin(varargin,'PBPSpec');
[~,sigUnits] = Utilities.ProcVarargin(varargin,'sigUnits');
[~,ActSpec] = Utilities.ProcVarargin(varargin,'ActSpec');
% [~,ConditionFieldLabel] = Utilities.ProcVarargin(varargin,'ConditionFieldLabel','Condition');
% [~,ObsPerCond] = Utilities.ProcVarargin(varargin,'ObsPerCond',10);
[~,FileTag] = Utilities.ProcVarargin(varargin,'FileTag','');
[~,saveFig] = Utilities.ProcVarargin(varargin,'saveFig');
[~,withShuffle] = Utilities.ProcVarargin(varargin,'withShuffle');
[~,Conditions] = Utilities.ProcVarargin(varargin,'Conditions',{});
[~,MinJerk] = Utilities.ProcVarargin(varargin,'MinJerk',1);



StyleArgs = {'FontSize',9,'FontName','Arial','FontWeight','bold'};

DecodeLabels = [Conditions{:}];

%% Only include units with tunign to at least one of the relevant variables
Dates = Analyze.returnUniqueFieldValues(AllTrialData,'Date');
if PBPSpec
    ASF = Analyze.FaceScratch3.FindSignificantUnitsAction(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,'Dates',Dates,'UnitGroup','PBPSpec');
elseif ActSpec
    ASF = Analyze.FaceScratch3.FindSignificantUnitsAction(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,'Dates',Dates,'UnitGroup','ActSpec');
elseif sigUnits
    ASF = Analyze.FaceScratch3.FindSignificantUnitsAction(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,'Dates',Dates,'UnitGroup','sigUnits');
else
    ASF = Analyze.FaceScratch3.FindSignificantUnitsAction(Tag,'PThresh',PThresh,...
        'Labels',Labels,'Labels2Incl',DecodeLabels,'Dates',Dates,'UnitGroup','All');
end

%%

opts.nReps=50;
opts.winSize=1.3;
opts.timeWindow=TimeWindow;

opts.NumShuffles=200;
opts.cvOptions.ValidationType='ClassicCrossValidation';
opts.cvOptions.NReps=30;
opts.cvOptions.NFolds=10;
opts.dec=Predictor.FWClassifier(@Analyze.FaceScratch3.BasicClassifier);

%%

useSplits=1;
nReps=opts.nReps;
winSize=opts.winSize;
timeWindow=opts.timeWindow;

kern=Kinematics.MinJerkKernel(winSize+.1, .1);
MinJerk=MinJerk;


%% first do within and across format classification of Action Selectivity

for timeIDX=1:length(timeWindow)
    clear FR
    
    for condIDX=1:length(Conditions)
        cCond=Conditions{condIDX};
        cTrialData=Analyze.SubSelectTrials(AllTrialData,'Phase',Phase,'Condition',cCond);
        
        L{condIDX}=Analyze.returnFieldValues(cTrialData,'Action');
        
        if MinJerk
            tmp=squeeze(Analyze.getNeuralData(cTrialData, ASF,[(timeWindow(timeIDX)-winSize/2-0.05):.1:(timeWindow(timeIDX)+winSize/2+.05)]));
            k2=repmat(kern',[size(tmp,1),1,size(tmp,3)]);
            tmp=tmp.*k2; tmp=squeeze(sum(tmp,2));
            FR{condIDX}=tmp;
        else
            FR{condIDX}=squeeze(Analyze.getNeuralData(cTrialData, ASF,[timeWindow(timeIDX)-winSize/2 timeWindow(timeIDX)+winSize/2]));
        end
        
        [FR{condIDX},L{condIDX}]=Analyze.SampleFeaturePopulation(FR{condIDX},L{condIDX},'Type','Basic');
        
        if 1==0
            [FRCell{condIDX,timeIDX},LabelCell{condIDX,timeIDX}]=DataConvert.list2cell(FR{condIDX},L{condIDX});
            FRCell{condIDX,timeIDX}= cellfun(@(x)x(1:10,:),FRCell{condIDX,timeIDX},'UniformOutput',false);
            [FRTrue{condIDX,timeIDX},LabelTrue{condIDX,timeIDX}]=DataConvert.cell2list(FRCell{condIDX,timeIDX});
        else
            [FRCell{condIDX,timeIDX},LabelCell{condIDX,timeIDX}]=DataConvert.list2cell(FR{condIDX},L{condIDX});
            mu=cellfun(@mean,FRCell{condIDX,timeIDX},'UniformOutput',false);
            mu=cat(1,mu{:});
            mn=min(mu,[],1);
            mx=max(mu,[],1);
            FR{condIDX}=(FR{condIDX}-mn)./(mx-mn);
            [FRCell{condIDX,timeIDX},LabelCell{condIDX,timeIDX}]=DataConvert.list2cell(FR{condIDX},L{condIDX});
            FRCell{condIDX,timeIDX}= cellfun(@(x)x(1:10,:),FRCell{condIDX,timeIDX},'UniformOutput',false);
            [FRTrue{condIDX,timeIDX},LabelTrue{condIDX,timeIDX}]=DataConvert.cell2list(FRCell{condIDX,timeIDX});
        end
    end
end

%%
for timeIDX1=1:length(timeWindow)
    for timeIDX2=1:length(timeWindow)
        
        disp(timeIDX2)
        for format1=1:4
            for format2=1:4
                d1=FRTrue{format1,timeIDX1};
                
                d2=FRTrue{format2,timeIDX2};
                l1=LabelTrue{format1,timeIDX1};
                l2=LabelTrue{format2,timeIDX1};
                clear l2_pred l2_test
                for rep=1:size(d1,1)
                    trainIDX=setdiff(1:size(d1,1),rep);
                    
                    d1_test=d1(rep,:);
                    d1_train=d1(trainIDX,:);
                    
                    l1_test=l1(rep,:);
                    l1_train=l1(trainIDX,:);
                    if rep>size(d2,1)
                        break
                    end
                    d2_test=d2(rep,:);
                    l2_test(rep)=l2(rep,:);
                    
                    [FR,L]=DataConvert.list2cell(d1_train,l1_train);
                    
                    opts.dec.Train('TrainingData',{FR,[1:4],{'1','2','3','4'}});
                    l2_pred(rep)=opts.dec.PredictBatch(d2_test);
                end
                %                 CompleteStats=confusionmatStats(l2_test,l2_pred);
                cvAccuracy{format1,format2}(timeIDX1,timeIDX2)=(nnz(l2_test==l2_pred)/length(l2_pred))*100;
                %                 NShuffs=1000;
                %                   for shuffIDX=1:NShuffs
                %                 ChancePercent(shuffIDX)=(nnz(l2_test==Shuffle(l2_pred))/length(l2_pred))*100;
                %             end
            end
        end
        
        
    end
end

Results.cvAccuracy=cvAccuracy;

%%

% F=Analyze.returnFieldValues(PlotData,'DecodeTest_cvAccuracy');
FOut=cvAccuracy;
% for sessionIDX=1:length(F)
%     for format1=1:4
%         for format2=1:4
%             if sessionIDX==1
%                 FOut{format1,format2}=[];
%             end
%             FOut{format1,format2}=cat(3,FOut{format1,format2},F{sessionIDX}{format1,format2});
%         end
%     end
% end
%%
for format1=1:4
    for format2=1:4
        
        FOutMu{format1,format2}=mean(FOut{format1,format2},3);
    end
end


%% To plot as heatmap style

Conditions={'Fc','Fs','Oc','Os'};
plt.fig('units','inches','width',8,'height',6,'font','Arial','fontsize',12);
pnl = panel();  pnl.margin=15; pnl.pack(4,4);

allVals=cat(1,FOutMu{:}); mx=max(allVals(:)); mn=min(allVals(:));
pk95=prctile(allVals(:),[5 95]);



for format1=1:4
    for format2=1:4
        
        pnl(format1,format2).select();
        if 1==1
            imagesc(FOutMu{format1,format2},[20 pk95(2)]); colorbar;
        else
            allVals=FOutMu{format1,format2}; mx=max(allVals(:)); mn=min(allVals(:));
            pk95=prctile(allVals(:),[30 95]);
            imagesc(FOutMu{format1,format2},[20 pk95(2)]);
        end
        colormap(parula);
        %         colorbar
        axis tight
        axis image
        set(gca, 'XTick',1:5:length(opts.timeWindow));
        set(gca, 'YTick',1:5:length(opts.timeWindow));
        
        %         set(gca, 'XTickLabel',opts.timeWindow(1:5:end)-opts.winSize/2)
        %         set(gca, 'YTickLabel',opts.timeWindow(1:5:end)-opts.winSize/2)
        
        set(gca, 'XTickLabel',opts.timeWindow(1:5:end));
        xtickangle(45);
        set(gca, 'YTickLabel',opts.timeWindow(1:5:end));
        
        xlabel(sprintf('From: %s ',Conditions{format1}));
        ylabel(sprintf('To: %s ',Conditions{format2}));
        
        plt.vline(3,{'y','linewidth',1});
        plt.hline(3,{'y','linewidth',1});
    end
end

%% To plot ERA style

Conditions={'Fc','Fs','Oc','Os'};
plt.fig('units','inches','width',8,'height',6,'font','Arial','fontsize',12);
pnl = panel();  pnl.margin=15; pnl.pack(4,4);

allVals=cat(1,FOutMu{:}); mx=max(allVals(:)); mn=min(allVals(:));
pk95=prctile(allVals(:),[5 95]);

FOutMuRescaled= cellfun(@(x)rescale(x,20, pk95(2)),FOutMu,'UniformOutput',false);

clr=[];
clr([1,2],:) = repmat([237 0 38]./255,2,1);
clr([3,4],:) = repmat([31 120 180]./255,2,1);

for format1=1:4
    for format2=1:4
        
        pnl(format1,format2).select();
        
        cData = FOutMuRescaled(format1,format2);
        Analyze.plotEventRelatedAverage(cData,{'a'},'TimeVec',timeWindow,'useBootstrap',...
            'Colors',clr(format1,:));
        
        ax1 = gca;
        set(ax1,'box','off',StyleArgs{:});
        ax1.XRuler.Axle.LineStyle = 'none';
        set(gca,'TickDir','out');
        set(ax1,'TickLength',[0 0]);
        axis tight;
        xlabel('Time (s)',StyleArgs{:});
        ylabel('From Baseline',StyleArgs{:});
        title(sprintf('%s to %s',Conditions{format1},Conditions{format2}));
    end
end

count=1;
for format1 = 1:4
    for format2  =1:4
        
        pnl(format1,format2).select()
        ylimtmp1(count,:) = ylim();
        count=count+1;
    end
end

for format1=1:4
    for format2=1:4
        
        pnl(format1,format2).select()
        ylim([min([ylimtmp1(:,1)]) max(ylimtmp1(:,2))]);
        
        tEvents = [0];
        plt.vline(tEvents,'k');
        ylims = ylim();
        ytext = ylims(2) + .05*diff(ylims);
        text(tEvents,ytext,'Go',StyleArgs{:});
    end
end

%%
filename = sprintf('CrossDecodeDynamic-%s-vs-%s',DimNames{1},DimNames{2});
if PBPSpec
    var = 'PBPSpec';
    filename = [filename '-PBPSpec'];
    suptitle(sprintf('PBPSpec units only (%d)',size(ASF,1)));
elseif ActSpec
    var = 'ActSpec';
    filename = [filename '-ActSpec'];
    suptitle(sprintf('ActSpec units only (%d)',size(ASF,1)));
elseif sigUnits
    var = 'sigUnits';
    filename = [filename '-sigUnits'];
    suptitle(sprintf('All significant units (%d)',size(ASF,1)));
else
    var = 'allUnits';
    filename = [filename '-allUnits'];
    suptitle(sprintf('All units (%d)',size(ASF,1)));
    
end

if ~isempty(FileTag)
    filename = [filename '-' FileTag];
end

FigsDir = fullfile(env.get('results'),'FaceScratch3','SUAnal',[Tag '-Go'],'PopData','Mine');
filename = [filename];
plt.SaveFigure(saveFig,FigsDir,[filename],'PNG','SVGI');
% plt.SaveFigure(saveFig,FigsDir,[filename],'SVGI');

end

