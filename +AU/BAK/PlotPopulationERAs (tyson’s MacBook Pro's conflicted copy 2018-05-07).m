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
opts.plotWaveform=Utilities.getInputField(opts,'plotWaveform', 1);

% if nargin < 4 % By default, plot all panel conditions found in Config.
CondIDX = ones(1,length(PanelCondNames));
% end
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
    
    if isfield(opts,'selectNSP')
        nspUnq = unique(ASF(:,1));
        if length(nspUnq) > 1
            ASF1 = ASF(ASF(:,1)==nspUnq(1,1),:);
            ASF2 = ASF(ASF(:,1)==nspUnq(2,1),:);
            fprintf('\n\n two nsps found \n\n')
            ASF = ASF1;
            if opts.selectNSP == 1
                ASF = ASF1;
            elseif opts.selectNSP == 2
                ASF = ASF2;
            end
        end
    end
    
    if isfield(opts,'limitWFSNR') && opts.limitWFSNR
        clear WF
        WFtmp=Analyze.getNeuralData(AllTrialData,ASF,'Waveform');
        if size(WFtmp,2)>1
            for i=1:size(WFtmp,1)
                a=WFtmp(i,:);
                a(cellfun(@length,a)==1)=[];
                WF{i}=cat(2,a{:});
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
        N=length(opts.ERAPlotFunction); if PanelCondIDX>N; ERAPlotidx=N; else, ERAPlotidx=PanelCondIDX; end
        h1=feval(opts.ERAPlotFunction{ERAPlotidx});
        
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
        for cUnitINDX=1:size(ASF,1);
            unit=ASF(cUnitINDX,:);
            opts.Quality=mean(Analyze.getNeuralData(AllTrialData,unit,'Quality'));
            Feature=sprintf('%d-%d-%d-%d',unit(1),unit(2),unit(3),opts.Quality);
            disp(Feature)
            
                if opts.plotWaveform
                  WFtmp=Analyze.getNeuralData(AllTrialData,unit,'Waveform');
                  opts.Waveform=cat(2,WFtmp{:});                  
                end
    
            TrialResults=[];
            
            
            % Set up figure
            figure(h1);
            clf; pnl = panel();
            if isfield(opts,'ERAMargins')
                pnl.margin=opts.ERAMargins;
            else
                pnl.margin=10;
            end
            if isfield(opts,'ERASideMargin')
                pnl.marginright=opts.ERASideMargin;
                pnl.marginleft=opts.ERASideMargin;
            end
            
            
            
            pnl.pack(length(cPanelCondVals),hPack);
            
            % Plot
            plotERA(AllTrialData,unit,cPanelCondName,...
                cPanelCondVals,cPlotConds,Phases,pnl,opts.ERAPlotBins,SmoothingKernel,opts);
            
            
            
            Analyze.AddFigLabel(AllTrialData,ASF(cUnitINDX,:),PlotType);
            if isfield(opts,'plotRaster') && opts.plotRaster
                saveName= sprintf('%s-%s-%s-RAST',opts.Subject,TaskDate,Feature);
            else
                SK=sprintf('%s-%0.2f',SmoothingKernel.mode,SmoothingKernel.value);
                saveName= sprintf('%s-%s-%s-ERA-%s',opts.Subject,TaskDate,Feature,SK);
            end
            
            if isfield(opts,'dejitter')
                saveName=[saveName '-DeJit'];
            end
            
            if isfield(opts, 'ERABaselinePhase')
                
                saveName=[saveName '-No' opts.ERABaselinePhase];
            end
            
            plt.SaveFigure(2,fullfile(opts.ERADir,PlotType),saveName,'PNG','SVG')
            
        end
    end
end
end

function plotERA(AllTrialData,cUnitIDX,PanelCondName,PanelCondVals,PlotConds,Phases,pnl,timeWindows,SmoothingKernel,opts)
%%
if isfield(opts,'plotRaster') && opts.plotRaster
    PlotRasters=1;
else
    PlotRasters=0;
end

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
        
        switch SmoothingKernel.mode
            case 'MinJerk'
                smoothArg=[SmoothingKernel.value mean(diff(timeWindows.(twField)))];
            case 'Box'
                smoothArg= SmoothingKernel.value;
            case 'BoxWin'
                smoothArg= round(SmoothingKernel.value/mean(diff(timeWindows.(twField))));
        end
        tmp=~cellfun(@isempty,PlotConds(PlotCondsRow,:));
        for jj = 1:length(PlotConds(PlotCondsRow,tmp))
            if strcmp(cPanelCondVal,'*')
                varargin = [{'Phase',cPhase} PlotConds{PlotCondsRow,jj}];
            else
                varargin = [{PanelCondName,cPanelCondVal,'Phase',cPhase} PlotConds{PlotCondsRow,jj}];
            end
            cTrialData = Analyze.SubSelectTrials(AllTrialData,varargin{:});
            [BinnedData{kk,ii,jj} SpikeData{kk,ii,jj}]= Analyze.getNeuralData(cTrialData,cUnitIDX,...
                timeWindows.(twField),'smooth',smoothArg);
            
            if isfield(opts, 'ERABaselinePhase') 
               A=varargin; A{find(strcmp(A,'Phase'))+1}=opts.ERABaselinePhase;
            cTrialData = Analyze.SubSelectTrials(AllTrialData,A{:});
            [BaseLine]= Analyze.getNeuralData(cTrialData,cUnitIDX,...
                opts.ERABaselineWindow);
            if isfield(opts,'ERABaselineZScore')
                if std(BaseLine) == 0
                    BinnedData{kk,ii,jj} = nan(size(BinnedData{kk,ii,jj}));
                    warning('Standard Deviation was found to be equal to 0...skipping and setting to NaN')
                    continue
                else
                BinnedData{kk,ii,jj}=bsxfun(@minus,BinnedData{kk,ii,jj},mean(BaseLine));
                BinnedData{kk,ii,jj} = BinnedData{kk,ii,jj}/std(BaseLine);
                end
            else
                BinnedData{kk,ii,jj}=bsxfun(@minus,BinnedData{kk,ii,jj},BaseLine);
            end
            end
            
        end
    end
end
%%
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
        
        if ~PlotRasters;
            PlotData=squeeze(BinnedData(kk,ii,:));
            tmp=~cellfun(@isempty,PlotConds(PlotCondsRow,:));
            Labels = getLegendLabels(PlotConds(PlotCondsRow,tmp),PlotData);
            
            if isfield(opts,'dejitter')
                maxShift=round(opts.dejitter/mean(diff(timeWindows.(twField))));
                for i=1:length(PlotData)
                    PlotData{i}=align2Mean(PlotData{i},'maxShift',maxShift);
                    PlotData{i}=PlotData{i}(:,maxShift:end-maxShift);
                end
                timeCenters=timeCenters(maxShift:end-maxShift);
            end
            
            
            
            warning off;
            if ii == 1
                Analyze.plotEventRelatedAverage(PlotData,Labels,'TimeVec',timeCenters,'Legend',Args{:});
            else
                Analyze.plotEventRelatedAverage(PlotData,Labels,'TimeVec',timeCenters,Args{:});
            end
            warning on;
        else
            PlotData=squeeze(SpikeData(kk,ii,:));
            Labels = getLegendLabels(PlotConds(PlotCondsRow,:),PlotData);
            
            if isfield(opts,'dejitter')
                BinData=squeeze(BinnedData(kk,ii,:));
                deltaT=mean(diff(timeWindows.(twField)));
                maxShift=round(opts.dejitter/deltaT);
                for i=1:length(PlotData)
                    [BinData{i},ShiftVals]=align2Mean(BinData{i},'maxShift',maxShift);
                    for trialIDX=1:length(PlotData{i})
                        PlotData{i}{trialIDX}=PlotData{i}{trialIDX}-ShiftVals(trialIDX)*deltaT;
                    end
                end
            end
            
            warning off;
            if ii == 1
                Analyze.plotEventRelatedRasters(PlotData,Labels,'Legend',Args{:});
            else
                Analyze.plotEventRelatedRasters(PlotData,Labels,Args{:});
            end
            warning on;
            
        end
        if ii~=1
            set(gca,'yTick',[])
        end       
        
        axis tight
        title(Phases{ii});
    end
end
%%
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
%             ylim(pnl(i,j).axis,[1 16]);
        end
    end
end

if isfield(opts,'vline');
    idx=1;
    for i=1:length(PanelCondVals)
        for j=1:length(Phases)
            %         ylim(pnl(i,j).axis,[0 max(Limits(:,4))]);
%             plt.vline(opts.vline{j},{'k'},'axesHandle',pnl(i,j).axis)
            plt.vline(opts.vline{j},{'k'},[],[],[],pnl(i,j).axis);
%             ylim(,[min(Limits(:,3)) max(Limits(:,4))]);
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

        if opts.plotWaveform
           axes('Position',[.85 .9 .15 .1]-0.01);
           Analyze.plotEventRelatedAverage({opts.Waveform'/4},1,'usePercentile','Colors',[1 0 0]); %'usePercentile' useBootStrap
           axis tight
           handle=title(['Q=' num2str(opts.Quality)]);
%            v=axis; 
%            set(handle,'Position',[10 v(4)*.6 0]);
            set(handle,'Position',[10 15 0]);
           set(gca,'XTickLabel',[]);
           %%
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
end