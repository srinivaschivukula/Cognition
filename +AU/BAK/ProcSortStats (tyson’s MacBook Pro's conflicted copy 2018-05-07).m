classdef ProcSortStats < handle
    
    % This class is designed to support analyses revolving around the sorting
    % statisitics. This includes:
    % Reurning unit statistics from a set of converted data files (across days)
    % and saving out the file using the subselect data.
    %
    % Ability to take saved data file and return relevant sort stats for any
    % arbitrary set of units.
    %
    % Plot sort statistics summary for any arbitrary set of units.
    
    
    % Useage example
    %     Opts.Experiment='OMASemantic'
    %             Opts.TaskName='ObsConcept';
    %             Opts.Subject='p1';
    %             Opts.TaskDates={'20170426','20170428','20170503','20170505','20170508','20170512','20170526','20170605','20170606','20170607','20170608','20170710', '20170711'};
    %     this=ProcSortStats(Opts);
    %%
    properties
        SortData
        TaskDates
        TaskName
        Subject
        Experiment
    end
    
    methods
        
        function this=ProcSortStats(Opts,varargin)
            % init
            this.SortData=[];
            this.Experiment=[];
            this.TaskName=[];
            this.TaskDates=[];
            this.Subject=[];
            
            fn = fieldnames(Opts);
            for i=1:length(fn)
                this.(fn{i})=Opts.(fn{i});
            end
            
            [varargin, this.SortData]   = Utilities.ProcVarargin(varargin,'SortData',this.TaskDates);
            [varargin, this.TaskDates]   = Utilities.ProcVarargin(varargin,'TaskDates',this.TaskDates);
            [varargin, this.TaskName]   = Utilities.ProcVarargin(varargin,'TaskName',this.TaskName);
            [varargin, this.Experiment]   = Utilities.ProcVarargin(varargin,'Experiment',this.Experiment);
            [varargin, this.Subject]   = Utilities.ProcVarargin(varargin,'TaskName',this.Subject);
            [varargin, Overwrite]   = Utilities.ProcVarargin(varargin,'Overwrite');
            
            
            fileloc=fullfile(env.get('result'),this.Experiment,'SortData');
            filename=sprintf('%s.%s.mat',this.TaskName,this.Subject);
            filepath=fullfile(fileloc,filename);
            
            if~exist(fileloc)
                mkdir(fileloc);
            end
            
            if Overwrite
                this.CreateDataStruct();
                SortData=this.SortData;
                save(filepath,'SortData')
            else
                % check if it exists
                if~exist(filepath)
                    this.CreateDataStruct();
                    SortData=this.SortData;
                    save(filepath,'SortData')
                else
                    tmp=load(filepath);
                    this.SortData=tmp.SortData;
                end
            end
            % check if it exists
        end
        
        function this=CreateDataStruct(this, varargin)
            
            [varargin,this.TaskDates]   = Utilities.ProcVarargin(varargin,'TaskDates',this.TaskDates);
            [varargin,this.TaskName]   = Utilities.ProcVarargin(varargin,'TaskName',this.TaskName);
            %%
            for dateIDX=1:length(this.TaskDates)
                TaskDate=this.TaskDates{dateIDX};
                AllTrialData=Analyze.LoadConvertedData(this.TaskName,TaskDate);
                
                ASF = Analyze.returnUniqueFieldValues(AllTrialData,'UnitIDsSorted');
                
                this.SortData{dateIDX}.ASF=ASF';
                this.SortData{dateIDX}.unit=ASF(:,1:3)';
                this.SortData{dateIDX}.DateID=ASF(:,4)';
                %                 this.SortData.DateStr=repmat(TaskDate,);
                this.SortData{dateIDX}.Quality = Analyze.getNeuralData(AllTrialData,ASF,'Quality');
                this.SortData{dateIDX}.PeakSNR = Analyze.getNeuralData(AllTrialData,ASF,'SNRPeak');
                this.SortData{dateIDX}.ISIViol = Analyze.getNeuralData(AllTrialData,ASF,'ISIViol');
                this.SortData{dateIDX}.CV2 = Analyze.getNeuralData(AllTrialData,ASF,'CV2');
                this.SortData{dateIDX}.IsolDist = Analyze.getNeuralData(AllTrialData,ASF,'IsolDist');
                this.SortData{dateIDX}.LRatio = Analyze.getNeuralData(AllTrialData,ASF,'LRatio');
                this.SortData{dateIDX}.NDimsSort = Analyze.getNeuralData(AllTrialData,ASF,'NDimsSort');
                this.SortData{dateIDX}.StdEstimate = Analyze.getNeuralData(AllTrialData,ASF,'StdEstimate');
                this.SortData{dateIDX}.MeanRate = nanmean(Analyze.getNeuralData(AllTrialData,ASF,'MeanRate')');
                
                clear WF WFm WFpr WFm
                tmp=Analyze.getNeuralData(AllTrialData,ASF,'Waveform');
                for j=1:size(tmp,1)
                    isValid=~cell2mat(cellfun(@(x)any(any(isnan(x))),tmp(j,:),'UniformOutput',false));
                    WF{j}=double(cat(2,tmp{j,isValid}));
                    
                    WFm(:,j)=mean(WF{j}');
                    %                    WFci{j}=bootci(10,@mean,WF{j}');
                    
                    WFpr{j}(1,:)=prctile(WF{j}',10);
                    WFpr{j}(2,:)=prctile(WF{j}',90);
                    
                end
                
                %                 this.SortData.Waveform = WF;
                this.SortData{dateIDX}.WaveformMean = WFm;
                %                 this.SortData{dateIDX}.WaveformCI = WFci;
                this.SortData{dateIDX}.WaveformPerc = WFpr;
                
                %                 SilVals = Analyze.getNeuralData(AllTrialData,ASF,'SilVal');
                %                 UnitType = Analyze.getNeuralData(AllTrialData,ASF,'UnitType');
            end
            
        end
        
        function returnSortStats(this,SortMeasure)
            
        end
        
        function plotSortStats(this, PlotType)
            %%
            Quality = Analyze.returnFieldValues(this.SortData,'Quality');
% SilVals =Analyze.returnFieldValues(this.SortData,'DateID');
PeakSNR = Analyze.returnFieldValues(this.SortData,'PeakSNR');
ISIViol = Analyze.returnFieldValues(this.SortData,'ISIViol');
CV2 =Analyze.returnFieldValues(this.SortData,'CV2');
IsolDist = Analyze.returnFieldValues(this.SortData,'IsolDist');
LRatio = Analyze.returnFieldValues(this.SortData,'LRatio');
NDimsSort = Analyze.returnFieldValues(this.SortData,'NDimsSort');
StdEstimate = Analyze.returnFieldValues(this.SortData,'StdEstimate');
            
            plt.fig('units','inches','width',10,'height',6,'font','Arial','fontsize',13);
pnl = panel(); pnl.margin = 15; pnl.pack(2,3);

% panelIdx = [1 2; 1 3; 2 1; 2 2; 2 3];
% vals = {ISIViolSUMU,PeakSNRSUMU,CV2SUMU,LRatioSUMU,IsolDistSUMU};
panelIdx = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3];
vals = {NDimsSort,ISIViol,PeakSNR,CV2,LRatio.*StdEstimate,log10(IsolDist)};
StatLabels = {'Sort Dimensions','%ISI < 3ms','Peak SNR','CV2','LRatio','Isolation Distance'};

bins={1:.1:5, 0:.1:6, 0:.25:15, 0:.05:2,0:1:40,.8:.04:3.5};
for i = 1:length(vals)
    pnl(panelIdx(i,1),panelIdx(i,2)).select();
%     histogram(vals{i}{1},bins{i}); hold on;
%     histogram(vals{i}{2},bins{i});
    histogram(vals{i},bins{i});
    xlabel(StatLabels{i});
    ylabel('# of units');
    xlim([bins{i}(1),bins{i}(end)]);
%     legend(Labels{:});
end
%%
        end
    end
end





% ConfigureSubject(opts.Subject);





%
% function PlotSortingStatsMUSU(varargin)
% %PLOTSORTINGSTATSMUSU Plot silhouette dists split by SU and MU
% %   SU = quality 1/2, MU = quality 3/4
% %   5/9/2017: Now also looks at more stats (using Ueli's code).
% TaskDates = {'20150904','20150909','20150914','20150923'};
% [varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');
%
%
% baseDir='C:\Users\cyzhang\Google Drive\CareyImagAttempt\';
%
% AllTrialData = Analyze.LoadConvertedData('MixedRep2',TaskDates);
% ASF = Analyze.returnUniqueFieldValues(AllTrialData,'UnitIDsSorted');
%
% Quality = Analyze.getNeuralData(AllTrialData,ASF,'Quality');
% SilVals = Analyze.getNeuralData(AllTrialData,ASF,'SilVal');
% PeakSNR = Analyze.getNeuralData(AllTrialData,ASF,'SNRPeakPooled');
% ISIViol = Analyze.getNeuralData(AllTrialData,ASF,'ISIViolPooled');
% CV2 = Analyze.getNeuralData(AllTrialData,ASF,'CV2Pooled');
% IsolDist = Analyze.getNeuralData(AllTrialData,ASF,'IsolDistPooled');
% LRatio = Analyze.getNeuralData(AllTrialData,ASF,'LRatioPooled');
% NDimsSort = Analyze.getNeuralData(AllTrialData,ASF,'NDimsSort');
% StdEstimate = Analyze.getNeuralData(AllTrialData,ASF,'StdEstimatePooled');
%
% QualClass = {[1 2],[3 4]};
% % QualClass = {[1],[3 4]};
% % QualClass = {[1:4],[1:4]};
% Labels = {'SU','MU'};
% for i = 1:length(Labels)
% %     cUnitIdx = ismember(Quality,QualClass{i});
%     if i == 1
%         cUnitIdx = IsolDist >= 10^1.6;
%     else
%         cUnitIdx = IsolDist < 10^1.6;
%     end
%     SilValsSUMU{i} = SilVals(cUnitIdx);
%     PeakSNRSUMU{i} = PeakSNR(cUnitIdx);
%     ISIViolSUMU{i} = ISIViol(cUnitIdx);
%     CV2SUMU{i} = CV2(cUnitIdx);
%     IsolDistSUMU{i} = log10(IsolDist(cUnitIdx));
%     LRatioSUMU{i} = LRatio(cUnitIdx);
% end
%
% %% Run some stats comparing the distributions
% [~,pKS] = kstest2(SilValsSUMU{1},SilValsSUMU{2});
% [~,pT,ciT,tStats] = ttest2(SilValsSUMU{1},SilValsSUMU{2});
%
% %% Plot
% plt.fig('units','inches','width',6,'height',5,'font','Arial','fontsize',13);
% pnl = panel(); pnl.margin = 15; pnl.pack(1,1);
% histogram(SilValsSUMU{1}); hold on;
% histogram(SilValsSUMU{2});
% title('Silhouette Value distributions');
% xlabel('Silhouette value');
% ylabel(sprintf('ttest2: p = %e',pT));
% legend(Labels{:});
%
% plt.SaveFigure(saveFig,fullfile(baseDir,'Figs2'),'SUMUSilValDist','PNG','SVGI');
%
% %% Plot other sorting stats
% plt.fig('units','inches','width',10,'height',6,'font','Arial','fontsize',13);
% pnl = panel(); pnl.margin = 15; pnl.pack(2,3);
%
% % panelIdx = [1 2; 1 3; 2 1; 2 2; 2 3];
% % vals = {ISIViolSUMU,PeakSNRSUMU,CV2SUMU,LRatioSUMU,IsolDistSUMU};
% panelIdx = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3];
% vals = {NDimsSort,ISIViol,PeakSNR,CV2,LRatio.*StdEstimate,log10(IsolDist)};
% StatLabels = {'Sort Dimensions','%ISI < 3ms','Peak SNR','CV2','LRatio','Isolation Distance'};
%
% bins={1:.1:5, 0:.1:6, 0:.25:15, 0:.05:2,0:1:40,.8:.04:3.5};
% for i = 1:length(vals)
%     pnl(panelIdx(i,1),panelIdx(i,2)).select();
% %     histogram(vals{i}{1},bins{i}); hold on;
% %     histogram(vals{i}{2},bins{i});
%     histogram(vals{i},bins{i});
%     xlabel(StatLabels{i});
%     ylabel('# of units');
%     xlim([bins{i}(1),bins{i}(end)]);
% %     legend(Labels{:});
% end
% plt.vline(1.6,'k');
%
% plt.SaveFigure(saveFig,fullfile(baseDir,'Figs2'),'SUMUStatsDist','PNG','SVGI');
%
% %%
% figure;
%
% X=[SilVals; ISIViol; PeakSNR ;CV2;LRatio;log10(IsolDist)]';
%
% plotmatrix(X)
%
% %
% % title('Silhouette Value distributions');
% % xlabel('Silhouette value');
% % ylabel(sprintf('ttest2: p = %e',pT));
% % legend(Labels{:});
%
% end
%
