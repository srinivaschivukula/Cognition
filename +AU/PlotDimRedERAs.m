function PlotPopulationERAs(config,varargin)
%%
if strcmp(mfilename,'') || nargin==0
    opts=Analyze.Libet.opts_BasicLinearModel;
else
    opts=feval(config,varargin{:});
end


if ~exist(opts.ERADir,'dir'); mkdir(opts.ERADir); end

% Get necessary fields from opts and stuff.
PanelCondNames = opts.ERAPanelCondNames; % Conditions to be separated into each panel.
PanelCondVals = opts.ERAPanelCondVals;
PlotConds = opts.ERAPlotConds; % Conditions to plot in each panel.
SmoothingKernel=Utilities.getInputField(opts,'ERASmoothingKernel', struct('mode','MinJerk','value',.35));

if nargin < 4 % By default, plot all panel conditions found in Config.
    CondIDX = ones(1,length(PanelCondNames));
end
%%
Phases=fieldnames(opts.ERAPlotBins);

for i=1:length(Phases)
    Dur(i)=opts.ERAPlotBins.(Phases{i})(end)-opts.ERAPlotBins.(Phases{i})(1);
end
PhaseSizes=Dur/sum(Dur);
for i=1:length(PhaseSizes); hPack{i}=PhaseSizes(i); end

for dateIDX=1:length(opts.Dates)
    TaskDate=opts.Dates{dateIDX};
    
    AllTrialData=Analyze.LoadConvertedData(opts.TaskFileName,TaskDate);
    
    % remove all trials that don't meet specified criteria
    if isfield(opts,'ERAGlobalTrialCriteria')
        AllTrialData=Analyze.SubSelectTrials(AllTrialData,opts.ERAGlobalTrialCriteria{:});
    end
    %%
    clear unitsAboveThresh notANan goodSNR
    if isnumeric(opts.unitType)
        ASF=opts.unitType;
    else
        % Get (A)ll (S)pike (F)eatures
        ASF=Analyze.returnUniqueFieldValues(AllTrialData,opts.unitType);
        
        if isfield(opts,'limitWFSNR') && opts.limitWFSNR
            clear WF
            WFtmp=Analyze.getNeuralData(AllTrialData,ASF,'Waveform');
            if size(WFtmp,2)>1
                for i=1:size(WFtmp,1)
                    a=WFtmp(i,:);
                    a(cellfun(@length,a)==1)=[];
                    WF{i}=cat(2,a{:});
                end
            else
                WF=WFtmp;
            end
            [snrval,FR,goodSNR]=computeSNR(WF','RateThresh',.05,'SNRThresh',.5);
            fprintf('%d/%d units meet waveform SNR criteria\n',nnz(goodSNR),length(goodSNR));
            goodSNR=goodSNR(:);
        else
            goodSNR=double(ones(size(ASF,1),1));
            goodSNR=goodSNR(:);
        end
        
        MR=mean(Analyze.getNeuralData(AllTrialData,ASF,'MeanRate'),2);
        unitsAboveThresh=MR>opts.FRthreshold;
        fprintf('%d/%d units meet firing rate criteria\n',nnz(unitsAboveThresh),length(unitsAboveThresh));
        
        notANan=~isnan(MR);
        fprintf('%d/%d are not NaNs\n',nnz(notANan),length(notANan));
        ASF=ASF(unitsAboveThresh & notANan & goodSNR,:);
    end
    %%
    for PanelCondIDX = 1:length(PanelCondNames)
        if ~CondIDX(PanelCondIDX)
            continue;
        end
        h1=feval(opts.ERAPlotFunction{PanelCondIDX});
        
        cPanelCondName = PanelCondNames{PanelCondIDX};
        cPanelCondVals = PanelCondVals{PanelCondIDX};
        cPlotConds = PlotConds{PanelCondIDX};
        
        PlotType = getPlotType(cPanelCondName,cPlotConds);
        if strcmp(cPanelCondName,'*')
            PlotType=PlotType(2:end);
        end
        warning off;
        mkdir(opts.ERADir,PlotType);
        warning on;
        
        
        % Set up figure
        figure(h1);
        clf; pnl = panel();
        if isfield(opts,'ERAMargins')
            pnl.margin=opts.ERAMargins;
        else
            pnl.margin=15;
        end
        pnl.pack(length(cPanelCondVals),hPack);
        
        % Plot
        plotERA(AllTrialData,ASF,cPanelCondName,...
            cPanelCondVals,cPlotConds,Phases,pnl,opts.ERAPlotBins,SmoothingKernel,opts);
        
        Analyze.AddFigLabel(AllTrialData,ASF(cUnitINDX,:),PlotType);
        if isfield(opts,'plotRaster') && opts.plotRaster
            saveName= sprintf('%s-%s-%s-RAST',opts.Subject,TaskDate,Feature);
        else
            saveName= sprintf('%s-%s-%s-ERA',opts.Subject,TaskDate,Feature);
        end
        
        if isfield(opts,'dejitter')
            saveName=[saveName '-DeJit'];
        end
        plt.SaveFigure(2,fullfile(opts.ERADir,PlotType),saveName,'PNG','SVGI')
        
    end
end
end


function plotERA(AllTrialData,cUnitIDX,PanelCondName,PanelCondVals,PlotConds,Phases,pnl,timeWindows,SmoothingKernel,opts)


if isfield(opts,'plotArgs')
    Args=opts.plotArgs;
else
    Args={'NoLine'};
end


TimeWindowFields = fields(timeWindows);

[nx,ny]=size(PlotConds);

if nx>1;
    assert(length(PanelCondVals)==nx)
end

% Get neural data for each of the traces.
for kk = 1:length(PanelCondVals)
    if  nx>1;
        PlotCondsRow=kk;
    else
        PlotCondsRow=1;
    end
    
    if iscell(PanelCondVals)
        cPanelCondVal = PanelCondVals{kk};
    else
        cPanelCondVal = PanelCondVals;
        
    end
    for ii = 1:length(Phases)
        cPhase = Phases{ii};
        twField = TimeWindowFields{ii};
        
        if SmoothingKernel.value==0
            SmoothData=0;
        else
            SmoothData=1;
        switch SmoothingKernel.mode
            case 'MinJerk'
                smoothArg=[SmoothingKernel.value mean(diff(timeWindows.(twField)))];
            case 'Box'
                smoothArg= SmoothingKernel.value;
            case 'BoxWin'
                smoothArg= round(SmoothingKernel.value/mean(diff(timeWindows.(twField))));
        end
        end
        
        for jj = 1:length(PlotConds)
            if strcmp(cPanelCondVal,'*')
                varargin = [{'Phase',cPhase} PlotConds{PlotCondsRow,jj}];
            else
                varargin = [{PanelCondName,cPanelCondVal,'Phase',cPhase} PlotConds{PlotCondsRow,jj}];
            end
            cTrialData = Analyze.SubSelectTrials(AllTrialData,varargin{:});
            if SmoothData
            [BinnedData{kk,ii,jj} SpikeData{kk,ii,jj}]= Analyze.getNeuralData(cTrialData,cUnitIDX,...
                timeWindows.(twField),'smooth',smoothArg);
            else
                [BinnedData{kk,ii,jj} SpikeData{kk,ii,jj}]= Analyze.getNeuralData(cTrialData,cUnitIDX,...
                timeWindows.(twField));
            end
        end
    end
end

% Plot
for kk = 1:length(PanelCondVals)
    if  nx>1;
        PlotCondsRow=kk;
    else
        PlotCondsRow=1;
    end
    for ii = 1:length(Phases)
        tW = timeWindows.(TimeWindowFields{ii});
        timeCenters = tW(1:end-1)+diff(tW);
        pnl(kk,ii).select();
        
            PlotData=squeeze(BinnedData(kk,ii,:));
            NTrials=cellfun(@(x)size(x,1),PlotData);
            
            
            % covnert data to gpfa format
            AllData=cat(1,PlotData{:});
            for i=1:size(AllData,1)
               AllDataStruct(i).trialId=i;
               AllDataStruct(i).spikes=squeeze(AllData(i,:,:))';
            end
            
            Labels = getLegendLabels(PlotConds(PlotCondsRow,:),PlotData);
            
            if isfield(opts,'dejitter')
                maxShift=round(opts.dejitter/mean(diff(timeWindows.(twField))));
                for i=1:length(PlotData)
                    PlotData{i}=align2Mean(PlotData{i},'maxShift',maxShift);
                    PlotData{i}=PlotData{i}(:,maxShift:end-maxShift);
                end
                timeCenters=timeCenters(maxShift:end-maxShift);
            end
            
            
           %%

            % Select method to extract neural trajectories:
            % 'gpfa' -- Gaussian-process factor analysis
            % 'fa'   -- Smooth and factor analysis
            % 'ppca' -- Smooth and probabilistic principal components analysis
            % 'pca'  -- Smooth and principal components analysis
            method = 'ppca';
            
            % Select number of latent dimensions
            xDim = 4;
            % NOTE: The optimal dimensionality should be found using
            %       cross-validation (Section 2) below.
            
            % If using a two-stage method ('fa', 'ppca', or 'pca'), select
            % standard deviation (in msec) of Gaussian smoothing kernel.
            kernSD = 75;
            % NOTE: The optimal kernel width should be found using
            %       cross-validation (Section 2) below.
             runIdx = 14;
            % Extract neural trajectories
            result = neuralTraj(runIdx, AllDataStruct, 'method', method, 'xDim', xDim,...
                'kernSDList', kernSD);
            
            [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);

            
            plotEachDimVsTime(seqTrain, 'xorth', result.binWidth);

             for i=1:size(AllData,1)
                 for j=1:result.xDim               
               AllDataRed(i,:,j)=seqTrain(i).xorth(j,:)';
               PlotData2{j}(i,:)=seqTrain(i).xorth(j,:)';
                 end
            end
            figure
            Analyze.plotEventRelatedAverage(PlotData2,{'1','2','3','4','5','6'},'Legend','NoLine' ,  'useBootStrap');
            
                    figure
            Analyze.plotEventRelatedAverage(PlotData2,{'1','2','3','4','5','6'},'Legend','NoBound');
%                             Analyze.plotEventRelatedAverage(PlotData2,{'1','2','3','4'},'Legend',Args{:});
%%
            
            warning off;
            if ii == 1
                Analyze.plotEventRelatedAverage(PlotData,Labels,'TimeVec',timeCenters,'Legend',Args{:});
            else
                Analyze.plotEventRelatedAverage(PlotData,Labels,'TimeVec',timeCenters,Args{:});
            end
            warning on;
        
        if ii~=1
            set(gca,'yTick',[])
        end
        axis tight
        title(Phases{ii});
    end
end

% Update axis to be same scale across all plots
if ~PlotRasters;
    idx=1;
    for i=1:length(PanelCondVals)
        for j=1:length(Phases)
            Limits(idx,:)=axis(pnl(i,j).axis);
            idx=idx+1;
        end
    end
    for i=1:length(PanelCondVals)
        for j=1:length(Phases)
            %         ylim(pnl(i,j).axis,[0 max(Limits(:,4))]);
            ylim(pnl(i,j).axis,[min(Limits(:,3)) max(Limits(:,4))]);
        end
    end
end
% Labels for each row.
for i = 1:length(PanelCondVals)
    if  ny>1;
        PlotCondsRow=kk;
    else
        PlotCondsRow=1;
    end
    
    pnl(i,1).select();
    if length(PanelCondVals)>1
        if strcmp(PanelCondVals{i},'*')
            ylabel(PlotConds{PlotCondsRow,jj}{1})
        else
            ylabel(PanelCondVals{i});
        end
    end
end

end

function PlotType = getPlotType(cPanelCondName,cPlotConds )
PanelName = extractNameFromPlotCond(cPanelCondName);
PlotType = sprintf('%sBy%s',PanelName,extractNameFromPlotCond(cPlotConds{1}{1}));
for i = 3:2:length(cPlotConds{1})
    PlotType = [PlotType 'And' extractNameFromPlotCond(cPlotConds{1}{i})];
end
end

function Labels = getLegendLabels(PlotConds,PlotData)
Labels = {};
for jj = 1:length(PlotConds)
    tempCond = cell(1,length(PlotConds{jj})/2);
    for c = 1:length(PlotConds{jj})/2
        [~,op] = extractNameFromPlotCond(PlotConds{jj}{2*c-1});
        cond = PlotConds{jj}{2*c};
        if isnumeric(cond)
            cond = num2str(cond);
        end
        tempCond{c} = sprintf('%s%s',op,cond);
    end
    cLabel = sprintf('%s: %d',strcat(tempCond{:}),size(PlotData{jj},1));
    Labels = [Labels cLabel];
end
end

function [Name,op] = extractNameFromPlotCond(plotCond)


if length(plotCond)>1 && strcmp('>=',plotCond(1:2))
    Name = plotCond(3:end);
    op = '>=';
elseif length(plotCond)>1 && strcmp('<=',plotCond(1:2))
    Name = plotCond(3:end);
    op = '<=';
elseif length(plotCond)>1 && strcmp('>',plotCond(1))
    Name = plotCond(2:end);
    op = '>';
elseif length(plotCond)>1 && strcmp('<',plotCond(1))
    Name = plotCond(2:end);
    op = '<';
else
    Name = plotCond;
    op = '';
end

end
