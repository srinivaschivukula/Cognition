% Address the basic question, when x is the salient category, what
% percentage of units are tuned to x.


function PlotPopulationAnalThroughTime(config,varargin)

if strcmp(mfilename,'') || nargin==0
    opts=Analyze.Libet.Config_BasicLinearModel;
else
    opts=feval(config,varargin{:});
end
%%

BaseDir=opts.ResultsDir;
TaskType=opts.TaskName;
TaskDate=opts.PlotDates;
Phases=opts.Phases;


if opts.PlotShuffle
    Prefix=sprintf('Shuffle-*-%s-', opts.PlotSubject);
else
    % Filename will include name passed to analysis function and subject ID
    Prefix=sprintf('%s-%s-',opts.AnalysisFCN{2}, opts.PlotSubject);
end

if ~isfield(opts,'PlotSaveType')
    PlotSaveType={'PNG','PDF'};
else
    PlotSaveType=opts.PlotSaveType;
end

%%
if strcmp(TaskDate,'*')
    fprintf('Loading from: %s\n',fullfile(BaseDir,[Prefix '*.mat']))
    FitList=dir(fullfile(BaseDir,['*' Prefix '*.mat']));
    Dates='AllDays';
else
    if iscell(TaskDate)
        FitList=[];
        for i=1:length(TaskDate)
            FitList=[FitList; dir(fullfile(BaseDir,[Prefix '' TaskDate{i} '*.mat']))];
            TaskDateSave{i}=[TaskDate{i}(3:end) '-'];
        end
        Dates=cat(2,TaskDateSave{:});
    else
        
        fprintf('Loading from: %s\n',fullfile(BaseDir,[Prefix '' TaskDate ,'*.mat']))
        
        FitList=dir(fullfile(BaseDir,[Prefix '' TaskDate ,'*.mat']));
        Dates=TaskDate;
    end
end
PopulationData=fullfile(BaseDir,'PopData'); if ~exist(PopulationData), mkdir(PopulationData); end

if strfind(opts.PlotSubject,'*')
    PrefixName=strrep(Prefix,'*','All');
else
    PrefixName=Prefix;
end

if isfield(opts,'Array') && opts.Array
    array=num2str(opts.Array);
    AllFiles={FitList.name}
    validFile=cellfun(@(x)~isempty(x),strfind(AllFiles,['Array' array]));
    FitList=FitList(validFile);
    ResultsFile=[PrefixName Dates '-Array' array];
    
else
    ResultsFile=[PrefixName Dates '-AllArrays'];
    
    end
    
    
    %% Gather data
    if exist(fullfile(PopulationData,[ResultsFile(3:end) '.mat']),'file') && ~opts.PopulationOverwrite
        fprintf('Population Data Already Generated: loading...\n');
        load(fullfile(PopulationData,ResultsFile(3:end)))
    else
        % Note Out is nPhases x nTime bins
        idx=1;
        fprintf('Loading %d files similar to %s \n', numel(FitList), (FitList(1).name))
        for fileIDX=1:length(FitList)
            clear TrialResults_NoTime TrialResults
            load(fullfile(BaseDir,FitList(fileIDX).name));
            
            for phaseIDX=1:length(Phases)
                cTrialResults=Analyze.SubSelectTrials(TrialResults,'Phase',Phases(phaseIDX));
                
                % get output field names
                allFields=fieldnames(TrialResults{1});
                tmp=strfind(allFields,[TrialResults{1}.Info.AnalysisFCN,'_']);
                opts.PlotFields=allFields(cell2mat(cellfun(@(x)isequal(x,1),tmp,'UniformOutput',false)));
                if isempty(opts.PlotFields)
                    tmp=allFields; tmp(ismember(tmp,{'Phase'    'Date'    'Unit'    'PhaseStart'    'PhaseEnd'    'Info'}))=[];
                    opts.PlotFields=tmp;
                end
                for dataIDX=1:length(opts.PlotFields)
                    PlotData{phaseIDX}.(opts.PlotFields{dataIDX})(:,fileIDX)=Analyze.returnFieldValues(cTrialResults,opts.PlotFields{dataIDX});
                end
                
                if exist('TrialResults_NoTime')
                    PlotFields_NoTime=fieldnames(TrialResults_NoTime{1});
                    for dataIDX=1:length(PlotFields_NoTime)
                        PlotData{phaseIDX}.(PlotFields_NoTime{dataIDX})(:,fileIDX)=Analyze.returnFieldValues(TrialResults_NoTime(1),PlotFields_NoTime{dataIDX});
                    end
                end
                
                if isfield(cTrialResults{phaseIDX}.Info,'UnitQuality')
                    PlotData{phaseIDX}.UnitQuality{fileIDX}=cTrialResults{phaseIDX}.Info.UnitQuality;
                end
                PlotData{phaseIDX}.UnitWaveform{fileIDX}=cTrialResults{phaseIDX}.Info.UnitWaveform;
                PlotData{phaseIDX}.TaskDate{fileIDX}=cTrialResults{phaseIDX}.Info.TaskDate;
                PlotData{phaseIDX}.UnitmeanRate{fileIDX}=cTrialResults{phaseIDX}.Info.UnitmeanRate;
                PlotData{phaseIDX}.Unit{fileIDX}=cTrialResults{phaseIDX}.Unit{1};
                
                PlotData{phaseIDX}.Phase{1,fileIDX}=Phases{phaseIDX};
                %%
                %     for i=1:length(FitResults)
                %         Population(fileIDX).CoefEst(i,:)=FitResults(i).Coefficients.Estimate;
                %         Population(fileIDX).CoefSE(i,:)=FitResults(i).Coefficients.SE;
                %         Population(fileIDX).CoeftStat(i,:)=FitResults(i).Coefficients.tStat;
                %         Population(fileIDX).CoefpValue(i,:)=FitResults(i).Coefficients.pValue;
                %     end
                Time{phaseIDX}.WindowStarts=Analyze.returnFieldValues(cTrialResults,'PhaseStart');
                Time{phaseIDX}.WindowStops=Analyze.returnFieldValues(cTrialResults,'PhaseEnd');
                Time{phaseIDX}.WindowMid=[Time{phaseIDX}.WindowStarts+Time{phaseIDX}.WindowStops]/2;
            end
            
        end
        save(fullfile(PopulationData,ResultsFile(3:end)),'PlotData','Time')
    end
    
    
    %%
    for i=1:length(Phases)
        Dur(i)=Time{i}.WindowStops(end)-Time{i}.WindowStarts(1);
    end
    PhaseSizes=Dur/sum(Dur);
    for i=1:length(PhaseSizes); hPack{i}=PhaseSizes(i); end
    
    %%
    for plotFCNidx=1:size(opts.PlotFCN,1)
        try
            
            pltFCNArgs=opts.PlotFCN{plotFCNidx,3};
            try
                % if a functon, this should run
                FigHandels = feval(opts.PlotFCN{plotFCNidx,1}, PlotData, pltFCNArgs{:});
            catch
                % but if it is a script, it will fail out and come here
                feval(opts.PlotFCN{plotFCNidx,1})
            end
            
            if opts.PlotResults && ~isempty(opts.PlotFCN{plotFCNidx,2})
                
                if ~exist('FigHandels') || length(FigHandels)==1 || isempty(FigHandels)
                    if ~isfield(opt,'Array')
                        array=num2str(opts.Array);
                    else
                        array = {'ArraysAll'};
                    end
                    plt.SaveFigure(1,PopulationData,[opts.PlotFCN{plotFCNidx,2} '-' PrefixName Dates array], PlotSaveType{:})
                else
%                     array=num2str(opts.Array);
                    for figIDX=1:length(FigHandels)
                        figure(FigHandels(figIDX))
                        % Filename will include name passed to plot function
                        plt.SaveFigure(1,PopulationData,[opts.PlotFCN{plotFCNidx,2}{:}, '-', PrefixName, Dates, '-', int2str(figIDX)], PlotSaveType{:})
                    end
                end
            end
        catch err
            if isfield(opts, 'FailLoudly') && opts.FailLoudly
                rethrow(err)
            else
                array = 'ArraysAll';
                
                plt.SaveFigure(1,PopulationData,[opts.PlotFCN{plotFCNidx,2} '-' PrefixName Dates array],PlotSaveType{:})
                warning(['Trouble with ' [opts.PlotFCN{plotFCNidx,2}, '-', PrefixName, Dates]])
                disp(getReport(err))
            end
        end
    end

    
    %
    %%
    % plt.fig('units','inches','width',10,'height',6,'font','Arial','fontsize',9);
    %
    % p = panel(); p.pack('h', hPack);
    %
    % for phaseIDX=1:length(Phases)
    %     p(phaseIDX).select(); hold on
    %     cPlotData=Analyze.SubSelectTrials(PlotData,'Phase',Phases{phaseIDX});
    % %     plotData=Analyze.returnFieldValues(cPlotData,'ChoiceP');
    % %     plotData1=Analyze.returnFieldValues(cPlotData,'InterP');
    % % PThresh=double(plotData<0.05 | plotData1<0.05)';
    %     plotData=Analyze.returnFieldValues(cPlotData,'ActionComp');
    % PThresh=double(plotData<0.05)';
    % PPerc=sum(PThresh,2)/size(PThresh,2);
    %     plot(Time{phaseIDX}.WindowMid,PPerc)
    %     ylim([0 .35])
    % end
    % disp('done!!')
    % %% Percent of population tuned
    % windowStarts=Population(1).windowStarts;
    % windowMid=windowStarts+Population(1).Duration/2;
    %
    % AllP=[Population.Pval];
    % PThresh=double(AllP<0.05);
    % PPerc=sum(PThresh,2)/size(PThresh,2);
    %
    % AllActionComp=[Population.ActionComp];
    % PThreshAllActionComp=double(AllActionComp<0.05 & logical(PThresh));
    % PPercAllActionComp=sum(PThreshAllActionComp,2)/size(PThreshAllActionComp,2);
    %
    % dur=windowStarts(end)-windowStarts(1);
    % h2=plt.fig('units','centimeters','width',dur,'height',8,'font','Arial','fontsize',9);
    % hold on;
    % area(windowMid,smooth(PPerc*100,3));
    % area(windowMid,smooth(PPercAllActionComp*100,3),'FaceColor',[0.5 0.9 0.6])
    % % area(windowMid,(PPerc*100));
    % % area(windowMid,(PPercAllActionComp*100),'FaceColor',[0.5 0.9 0.6])
    % axis tight
    % xlabel('Time (secs)')
    % ylabel('Population (%)')
    % title(sprintf('Population tuning'))
    % legend('Effector General','Effector Specific')
    % ylim([0 45])
    % FigName=sprintf('%sPopPerc-%s-%s',Prefix,Phase,Dates);
    % plt.SaveFigure(1,PopulationData,FigName, 'PNG','FIG','SVGI')
    % %%
    % % CoeftStat=[Population.CoeftStat];
    % % figure; hold on
    % % %
    % % A=max(CoeftStat(:,2:3:end)');
    % % B=max(CoeftStat(:,3:3:end)');
    % % plot(windowStarts,max(CoeftStat(:,2:3:end)'))
    % % plot(windowStarts,max(CoeftStat(:,3:3:end)'))
    % %
    % % plot(windowStarts,smooth(mean([A;B])),'linewidth',5)
    % % %%
    % % CoeftStat=[Population.CoeftStat];
    % % figure; hold on
    % %
    % % A=mean(abs(CoeftStat(:,2:3:end))');
    % % B=mean(abs(CoeftStat(:,3:3:end)'));
    % % plot(windowStarts,A')
    % % plot(windowStarts,B')
    % %
    % plot(windowStarts,smooth(mean([A;B])),'linewidth',5)
    % %
    %%
    % timeIDX=2
    % figure; hold on
    %
    % plt.ecdfPlot(FitP(:,timeIDX),'r')
    %
    % x = 0:.01:1;
    % f = unifcdf(x,0,1);
    % plot(x,f,'m')
    % axis tight; axis equal
    % grid minor
    % %%
    % timeIDX=2
    % figure; hold on
    %
    % plt.ecdfPlot(GH.pR(:,timeIDX),'r')
    % plt.ecdfPlot(GT.pR(:,timeIDX),'g')
    % plt.ecdfPlot(HT.pR(:,timeIDX),'b')
    %
    % x = 0:.01:1;
    % f = unifcdf(x,0,1);
    % plot(x,f,'m')
    % axis tight; axis equal
    % grid minor
    %
    %
    % %%
    % figure;
    % hold on
    % x = 0:.01:1;
    % wins=opts.windowStarts
    % for windowIDX=5:13
    %     subplot(3,3,windowIDX-4); hold on
    %     % cdfplot(Pval(:,windowIDX))
    %     % cdfplot(HLp(:,windowIDX))
    %     % cdfplot(ELp(:,windowIDX))
    %     % cdfplot(INTp(:,windowIDX))
    %
    %     [f,x,flo,fup]=ecdf(Pval(:,windowIDX),'bounds','on');
    %     fup(isnan(fup))=f(isnan(fup)); flo(isnan(flo))=f(isnan(flo)); ci=[abs(f-flo), abs(f-fup)];
    %     [hl,cl]=Utilities.boundedline(x,f,ci,'transparency',.5,'alpha','cmap',[1 0 0]);
    %     delete(hl)
    %
    %     x = 0:.01:1;
    %     f = unifcdf(x,0,1);
    %     plot(x,f,'m')
    %     title(sprintf('T:%0.2f(%0.2f)',wins(windowIDX),opts.Duration))
    % end
    % legend({'Pval','Null'},'Location','NW')
    %
    % plt.SaveFigure(2,PopulationData,'FullModel', 'PNG','SVG')
    %
    % %%
    % figure;
    % hold on
    % x = 0:.01:1;
    % wins=opts.windowStarts
    % for windowIDX=5:13
    %     subplot(3,3,windowIDX-4); hold on
    %
    %     [f,x,flo,fup]=ecdf(HLp(:,windowIDX),'bounds','on');
    %     fup(isnan(fup))=f(isnan(fup)); flo(isnan(flo))=f(isnan(flo)); ci=[abs(f-flo), abs(f-fup)];
    %     [hl,cl]=Utilities.boundedline(x,f,ci,'transparency',.5,'alpha','cmap',[0 1 0]);
    %     delete(hl)
    %
    %
    %     x = 0:.01:1;
    %     f = unifcdf(x,0,1);
    %     plot(x,f,'m')
    %     title(sprintf('T:%0.2f(%0.2f)',wins(windowIDX),opts.Duration))
    % end
    % legend({'HLp','Null'},'Location','NW')
    % plt.SaveFigure(2,PopulationData,'Hand', 'PNG','SVG')
    % %%
    % %%
    % figure;
    % hold on
    % x = 0:.01:1;
    % wins=opts.windowStarts
    % for windowIDX=5:13
    %     subplot(3,3,windowIDX-4); hold on
    %
    %
    %     [f,x,flo,fup]=ecdf(ELp(:,windowIDX),'bounds','on');
    %     fup(isnan(fup))=f(isnan(fup)); flo(isnan(flo))=f(isnan(flo)); ci=[abs(f-flo), abs(f-fup)];
    %     [hl,cl]=Utilities.boundedline(x,f,ci,'transparency',.5,'alpha','cmap',[0 0 1]);
    %     delete(hl)
    %
    %
    %     x = 0:.01:1;
    %     f = unifcdf(x,0,1);
    %     plot(x,f,'m')
    %     title(sprintf('T:%0.2f(%0.2f)',wins(windowIDX),opts.Duration))
    % end
    % legend({'ELp','Null'},'Location','NW')
    %
    % plt.SaveFigure(2,PopulationData,'Eye', 'PNG','SVG')
    % %%
    % figure;
    % hold on
    % x = 0:.01:1;
    % wins=opts.windowStarts
    % for windowIDX=5:13
    %     subplot(3,3,windowIDX-4); hold on
    %
    %
    %     [f,x,flo,fup]=ecdf(INTp(:,windowIDX),'bounds','on');
    %     fup(isnan(fup))=f(isnan(fup)); flo(isnan(flo))=f(isnan(flo)); ci=[abs(f-flo), abs(f-fup)];
    %     [hl,cl]=Utilities.boundedline(x,f,ci,'transparency',.5,'alpha','cmap',[.5 .5 .5]);
    %     delete(hl)
    %
    %     x = 0:.01:1;
    %     f = unifcdf(x,0,1);
    %     plot(x,f,'m')
    %     title(sprintf('T:%0.2f(%0.2f)',wins(windowIDX),opts.Duration))
    % end
    % legend({'INTp','Null'},'Location','NW')
    %
    % plt.SaveFigure(2,PopulationData,'Interaction', 'PNG','SVG')
    % % %%
    % % figure;
    % % hold on
    % % x = 0:.01:1;
    % %
    % % for windowIDX=1:length(opts.Phases2Process)
    % %
    % % % [f_,x_] = ecdf(Pval(:,windowIDX));
    % % % f(:,windowIDX)=interp1(x_(2:end),f_(2:end),x);
    % % cdfplot(INTp(:,windowIDX))
    % % end
    % %
    % % %%
    % % x = 0:.01:1;
    % % f = unifcdf(x,0,1);
    % plot(x,f,'m')
    % legend({opts.Phases2Process{:},'Null'},'Location','NW')
    %
