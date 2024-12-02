% Address the basic question, when x is the salient category, what
% percentage of units are tuned to x.


function SlidingTimeDPrimeMotorVersionStatic(config,varargin)

[varargin,DateName]   = Utilities.ProcVarargin(varargin,'DateName',[]);
[varargin,Pmethod]   = Utilities.ProcVarargin(varargin,'Pmethod','fdr');
[varargin,saveFig]=Utilities.ProcVarargin(varargin,'saveFig');
[varargin,DatesAll]   = Utilities.ProcVarargin(varargin,'Dates',{});
[varargin,trialdata]   = Utilities.ProcVarargin(varargin,'AllTrialData', {});
[varargin,ImagTrialData]   = Utilities.ProcVarargin(varargin,'ImagTrialData', {});
[varargin,FigSize] = Utilities.ProcVarargin(varargin,'FigSize',[]);
[varargin,SmoothingKernel]=Utilities.ProcVarargin(varargin,'SmoothingKernel', struct('mode','MinJerk','value',.5));
[varargin,TimeStep]   = Utilities.ProcVarargin(varargin,'TimeStep',0.1);
[varargin,StepSize]   = Utilities.ProcVarargin(varargin,'StepSize',0.3);
[varargin,SubtractBaseline]   = Utilities.ProcVarargin(varargin,'SubtractBaseline',0);
[varargin,NumSamples]   = Utilities.ProcVarargin(varargin,'NumSamples',15);



if isa(config,'function_handle') || isstr(config)
    opts=feval(config,varargin{:});
else
    opts=config;
end

if ~isfield(opts,'PlotSaveType')
    PlotSaveType={'PDF','PNG','PDF'};
else
    PlotSaveType=opts.PlotSaveType;
end
%%
TaskLabel = opts.GlobalTrialCriteria{2};
BaseDir=opts.ResultsDir;
TaskType=opts.TaskName;
TaskDate=opts.PlotDates;
Phases=opts.Phases;

Prefix=sprintf('%s-', opts.PlotSubject);

array='';

%%
% fprintf('Searching %s \n',BaseDir)
if strcmp(TaskDate,'*')
    FitList=dir(fullfile(BaseDir,[Prefix '*' array '-*.mat']));
    Dates='AllDays';
else
    if iscell(TaskDate)
        FitList=[];
        for i=1:length(TaskDate)
            FitList=[FitList; dir(fullfile(BaseDir,[Prefix '' TaskDate{i} array '-*.mat']))];
            TaskDateSave{i}=[TaskDate{i}(3:end) '-'];
        end
        Dates=cat(2,TaskDateSave{:});
    else
        FitList=dir(fullfile(BaseDir,[Prefix '' TaskDate , array,'*.mat']));
        Dates=TaskDate;
    end
end
PopulationData=fullfile(BaseDir,'PopData'); if ~exist(PopulationData), mkdir(PopulationData); end

if strfind(opts.PlotSubject,'*')
    PrefixName=strrep(Prefix,'*','All');
else
    PrefixName=Prefix;
end

if isempty(DateName)
    ResultsFile=[PrefixName Dates array];
else
    ResultsFile=[PrefixName DateName  array];
    
end

%% Gather data
load(fullfile(PopulationData,ResultsFile))


%%
try
    for i=1:length(Phases)
        Dur(i)=Time{i}.WindowStops(end)-Time{i}.WindowStarts(1);
    end
    PhaseSizes=Dur/sum(Dur);
    for i=1:length(PhaseSizes); hPack{i}=PhaseSizes(i); end
    
catch
    warning('WindowStops and WindowStarts not defined. do your own damn sizing.')
end

%% If option enabled, only look at well tuned units
if isfield(opts,'UseSparse') && opts.UseSparse
    PlotData = Analyze.SubSelectTrials(PlotData,'>=testdPBW',2);
end

%%
if isfield(opts,'unitWFDPThresh') && ~isempty(opts.unitWFDPThresh)
    if isnan(opts.unitWFDPThresh)
        unitWFDP=Analyze.returnFieldValues(PlotData(1),'unitWFDP');
        figure; histogram(unitWFDP,round(length(unitWFDP)/10))
        opts.unitWFDPThresh=prctile(unitWFDP,50);
        title('unitWFDP')
        plt.vline(opts.unitWFDPThresh);
    end
    PlotData=Analyze.SubSelectTrials(PlotData,'>=unitWFDP',opts.unitWFDPThresh);
    
    PrefixName=[PrefixName sprintf('WFDP%0.1f', opts.unitWFDPThresh) ,'-'];
end

plotFCNidx=1;

pltFCNArgs=opts.PlotFCN{plotFCNidx,3};

[~,Labels]   = Utilities.ProcVarargin(pltFCNArgs,'Labels','');
[~,PThresh]   = Utilities.ProcVarargin(pltFCNArgs,'PThresh',0.05);
SavePre=[opts.PlotFCN{plotFCNidx,2} '-' PrefixName Dates array];

if isfield(opts,'UseSparse') && opts.UseSparse
    SavePre=[opts.PlotFCN{plotFCNidx,2} '-' PrefixName 'sparse-' Dates];
else
    SavePre=[opts.PlotFCN{plotFCNidx,2} '-' PrefixName Dates];
end
SavePre;
[~,Phase]   = Utilities.ProcVarargin(pltFCNArgs,'Phase','Go');
PlotData=Analyze.SubSelectTrials(PlotData,'Phase',Phase);

ActionPvals=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'Pvals'));
ActionCoef=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'Coef'));
NC=size(ActionPvals,2);

pAnova=Analyze.returnFieldValues(PlotData,'pAnova');
[H,effAlphaVal,adj_p]=Utilities.MultipleComparisonsCorrection(pAnova,'method','fdr');
[Hhb,effAlphaValHB]=Utilities.MultipleComparisonsCorrection(pAnova,'method','holm-bonferroni');
R2=Analyze.returnFieldValues(PlotData,'R2adj');

dPrime=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'ComdPrime'));
AUC=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'ComAUC'));
AUCvAll=DataConvert.embCell2mat(Analyze.returnFieldValues(PlotData,'AUCvAll'));


ModelP=Analyze.returnFieldValues(PlotData,'pFull');
PopTuned=ModelP<PThresh;

%% Plotting stuff

pVals=cat(2,PlotData{1}.Pvals{:})';
if strcmp(Pmethod, 'holm-bonferroni')
    [Hhb,effAlphaValHB]=Utilities.MultipleComparisonsCorrection(pVals,'method',Pmethod);
    PThresh = effAlphaValHB;
elseif strcmp(Pmethod, 'fdr')
    [H,effAlphaVal,adj_p]=Utilities.MultipleComparisonsCorrection(pVals,'method',Pmethod);
    PThresh = effAlphaVal;
else
    [H,effAlphaVal,adj_p]=Utilities.MultipleComparisonsCorrection(pVals,'method',Pmethod);
    PThresh = min([0.001 effAlphaVal]);
end

ActionTuned=ActionPvals<PThresh;

%% Calculating the time-resolved distance measures for 100 randomly sampled units, 100 times

%% this section is for att vs imagined data
for sampIdx=1:NumSamples
    disp(sampIdx)
    %% first do all the att/imagined stuff
    
TimeWindows = [0.5 3];
%     TimeWindows = [-5:TimeStep:5];
    %     TimeWindows = [-4.5:TimeStep:5.5];
    %     TimeWindows = [-1.5:TimeStep:7.5];
    %     TimeCenters = TimeWindows+TimeStep;
%     TimeCenters=[];
%     for timeIdx = 1:length(TimeWindows)
%         cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
%         TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
%     end
    
    %% this portion is for att vs imagined data
    AllUnits = vertcat(PlotData{1}.Unit{:});
    randSamp = randsample(length(AllUnits),100);
    cUnits = AllUnits(randSamp,:);
    
    %@@@@%@@@@%
    %     dateCode = Blackrock.Helper.date2unitId(DatesAll{sampIdx});
    %     cUnits = AllUnits(AllUnits(:,4)==dateCode,:);
    %@@@@%@@@@%
    
    
    condNames = {'AttemptHandRight','AttemptShoulderRight','ImagineHandRight','ImagineShoulderRight'};
    
    Attempt_hand = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'AttemptHandRight'});
    Attempt_shoulder = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'AttemptShoulderRight'});
    
    Imag_hand = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'ImagineHandRight'});
    Imag_shoulder = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'ImagineShoulderRight'});
    
    BaselineData_attHand = Analyze.SubSelectTrials(trialdata, 'Phase', 'ITI','Condition',{'AttemptHandRight'});
    BaselineData_attShoulder = Analyze.SubSelectTrials(trialdata, 'Phase', 'ITI','Condition',{'AttemptShoulderRight'});
    BaselineData_imagHand = Analyze.SubSelectTrials(trialdata, 'Phase', 'ITI','Condition',{'ImagineHandRight'});
    BaselineData_imagShoulder = Analyze.SubSelectTrials(trialdata, 'Phase', 'ITI','Condition',{'ImagineShoulderRight'});
    %     BaselineDataFI_null = Analyze.SubSelectTrials(trialdata, 'Phase', 'ITI','Condition',{'XX'});
    
%         condNames = {'AttemptHandLeft','AttemptShoulderLeft','ImagineHandLeft','ImagineShoulderRight'};
%     
%     Attempt_hand = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'AttemptHandLeft'});
%     Attempt_shoulder = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'AttemptShoulderLeft'});
%     
%     Imag_hand = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'ImagineHandLeft'});
%     Imag_shoulder = Analyze.SubSelectTrials(trialdata, 'Phase', 'Go', 'Condition', {'ImagineShoulderRight'});
%     
%     BaselineData_attHand = Analyze.SubSelectTrials(trialdata, 'Phase', 'ITI','Condition',{'AttemptHandLeft'});
%     BaselineData_attShoulder = Analyze.SubSelectTrials(trialdata, 'Phase', 'ITI','Condition',{'AttemptShoulderLeft'});
%     BaselineData_imagHand = Analyze.SubSelectTrials(trialdata, 'Phase', 'ITI','Condition',{'ImagineHandLeft'});
%     BaselineData_imagShoulder = Analyze.SubSelectTrials(trialdata, 'Phase', 'ITI','Condition',{'ImagineShoulderRight'});
%     %     BaselineDataFI_null = Analyze.SubSelectTrials(trialdata, 'Phase', 'ITI','Condition',{'XX'});

    
    %% do the sliding time calculation of distances/dprimes
    
    for timeIdx = 1:length(TimeWindows)
        
cTW = [0 3];        
% cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        
        %% get the att/imag data
        % add prevTrialLabel
%         for unitIDX=1:size(AllUnits,1)
%                      cTrialData=Analyze.SubSelectTrials(trialdata,'Phase',Phase);
% 
%                     FR = [squeeze(Analyze.getNeuralData(cTrialData, AllUnits(unitIDX,:), cTW))];
%                     L=Analyze.returnFieldValues(cTrialData,'Condition');
%                     
%                     validIDX=~isnan(FR);
%                     FR=FR(validIDX);
%                     L=L(validIDX);
%                     
%                    tmp=table(L,FR);    
%                    tmp.L=nominal(tmp.L);
%                    mdl2{unitIDX}=fitlm(tmp);
%                   p(unitIDX)= mdl2{unitIDX}.coefTest;
%         end
%              
%     end
%                   Res= mdl.Residuals.Raw;
%                   
%                   tmp.Res=Res;
%                   
%                   FR_attHand(unitIDX,:)=Res(tmp.L==AttemptHandRight)
% 
%         end
        
        
        FR_attHand = [squeeze(Analyze.getNeuralData(Attempt_hand, cUnits, cTW))];
        FR_attShoulder = [squeeze(Analyze.getNeuralData(Attempt_shoulder, cUnits, cTW))];
        
        FR_imagHand = [squeeze(Analyze.getNeuralData(Imag_hand, cUnits, cTW))];
        FR_imagShoulder = [squeeze(Analyze.getNeuralData(Imag_shoulder, cUnits, cTW))];
        
%         if cTW( 1)>0
%             baseFR = [squeeze(Analyze.getNeuralData(BaselineData_attHand, cUnits, [0 3]));squeeze(Analyze.getNeuralData(BaselineData_attShoulder, cUnits, [0 3]));...
%                 squeeze(Analyze.getNeuralData(BaselineData_imagHand, cUnits, [0 3]));squeeze(Analyze.getNeuralData(BaselineData_imagShoulder, cUnits, [0 3]))];
%         else
            baseFR = [squeeze(Analyze.getNeuralData(BaselineData_attHand, cUnits, [2 3]));squeeze(Analyze.getNeuralData(BaselineData_attShoulder, cUnits, [2 3]));...
                squeeze(Analyze.getNeuralData(BaselineData_imagHand, cUnits, [2 3]));squeeze(Analyze.getNeuralData(BaselineData_imagShoulder, cUnits, [2 3]))];
            %             baseFR = [squeeze(Analyze.getNeuralData(BaselineData_attHand, cUnits, cTW));squeeze(Analyze.getNeuralData(BaselineData_attShoulder, cUnits, cTW));...
            %                 squeeze(Analyze.getNeuralData(BaselineData_imagHand, cUnits, cTW));squeeze(Analyze.getNeuralData(BaselineData_imagShoulder, cUnits, cTW))];
%         end
        
        newFR_attHand =[]; newFR_attShoulder =[]; newFR_imagHand = []; newFR_imagShoulder = []; newbaseFR=[];
        for i=1:size(FR_attHand,2)
            tmp = FR_attHand(:,i);
            abc = tmp(~isnan(tmp));
            newFR_attHand(:,i) = [abc; repmat(mean(abc),12-length(abc),1)];
            
            tmp = FR_attShoulder(:,i);
            abc = tmp(~isnan(tmp));
            newFR_attShoulder(:,i) = [abc; repmat(mean(abc),12-length(abc),1)];
            
            tmp = FR_imagHand(:,i);
            abc = tmp(~isnan(tmp));
            newFR_imagHand(:,i) = [abc; repmat(mean(abc),12-length(abc),1)];
            
            tmp = FR_imagShoulder(:,i);
            abc = tmp(~isnan(tmp));
            newFR_imagShoulder(:,i) = [abc; repmat(mean(abc),12-length(abc),1)];
            
            tmp = baseFR(:,i);
            abc = tmp(~isnan(tmp));
            newbaseFR(:,i) = [abc; repmat(mean(abc),48-length(abc),1)];
        end
        
        FR_attHand = newFR_attHand; FR_attShoulder = newFR_attShoulder; FR_imagHand=newFR_imagHand; FR_imagShoulder=newFR_imagShoulder; baseFR=newbaseFR;
        
        %% calculate all the distances/dprimes
        
        ResultsAtt_hand = Analyze.FitPLSLM([FR_attHand;baseFR],[ones(size(FR_attHand,1),1);zeros(size(baseFR,1),1)],'CrossValidate',1);
        ResultsAtt_shoulder = Analyze.FitPLSLM([FR_attShoulder;baseFR],[ones(size(FR_attShoulder,1),1);zeros(size(baseFR,1),1)],'CrossValidate',1);
        ResultsImag_hand = Analyze.FitPLSLM([FR_imagHand;baseFR],[ones(size(FR_imagHand,1),1);zeros(size(baseFR,1),1)],'CrossValidate',1);
        ResultsImag_shoulder = Analyze.FitPLSLM([FR_imagShoulder;baseFR],[ones(size(FR_imagShoulder,1),1);zeros(size(baseFR,1),1)],'CrossValidate',1);
        
        DataAtt_hand{sampIdx}{timeIdx} = ResultsAtt_hand.Dist;
        DataAtt_shoulder{sampIdx}{timeIdx} = ResultsAtt_shoulder.Dist;
        DataImag_hand{sampIdx}{timeIdx} = ResultsImag_hand.Dist;
        DataImag_shoulder{sampIdx}{timeIdx} = ResultsImag_shoulder.Dist;
        
    end
    
end

%% post processing and plotting

%% plot felt vs imagined

Att_hand=[]; Att_shoulder=[]; Imag_hand=[]; Imag_shoulder=[];
for i=1:length(DataAtt_hand)
    
    tmp = [DataAtt_hand{i}{:}];
    Att_hand = cat(1,Att_hand,tmp);
    
    tmp = [DataAtt_shoulder{i}{:}];
    Att_shoulder = cat(1,Att_shoulder,tmp);
    
    tmp = [DataImag_hand{i}{:}];
    Imag_hand = cat(1,Imag_hand,tmp);
    
    tmp = [DataImag_shoulder{i}{:}];
    Imag_shoulder = cat(1,Imag_shoulder,tmp);
end
PlotData = {{Att_hand},{Imag_hand},{Att_shoulder},{Imag_shoulder}};

StyleArgs = {'FontName','Arial','FontSize',9,'FontWeight','bold'};

figH = plt.fig('units','inches','width',14,'height',8,'font','Arial','fontsize',9);
figure(figH);
clf; pnl = panel(); pnl.margin=15; pnl.pack(3,4);

condNames = {'Att Hand','Imag Hand','Att Shoulder','Imag Shoulder'};

for condIdx=1:4
    
    cData = PlotData{condIdx};
    
    clr=[];
    clr(1,:) = repmat([237 0 38]./255,1,1);
    clr(2,:) = repmat([31 120 180]./255,1,1);
    clr(3,:) = repmat([237 0 38]./255,1,1);
    clr(4,:) = repmat([1 108 89]./255,1,1);
    
    cLabel = condNames{condIdx};
    
    TimeWindows = [-5:TimeStep:5];
    TimeCenters=[];
    for timeIdx = 1:length(TimeWindows)
        cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
        TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
    end
    
    % plot the discriminability of touch from baseline
    pnl(1,condIdx).select()
    
    Analyze.plotEventRelatedAverage(cData,cLabel,'TimeVec',TimeCenters,'usePercentile',...
        'Colors',clr(condIdx,:));
    title(cLabel);
    ax1 = gca;
    set(ax1,'box','off',StyleArgs{:});
    ax1.XRuler.Axle.LineStyle = 'none';
    set(gca,'TickDir','out');
    set(ax1,'TickLength',[0 0]);
    axis tight;
    xlabel('Time (s)',StyleArgs{:});
    ylabel('From Baseline',StyleArgs{:});
    
    xl = xlim;
    xlim([-4.5 4])
    
end

for plotIdx = 1:4
    
    pnl(1,plotIdx).select()
    ylimtmp1(plotIdx,:) = ylim();
    
end

for plotIdx = 1:4
    
    pnl(1,plotIdx).select()
    ylim([min([ylimtmp1(:,1);0]) max(ylimtmp1(:,2))]);
    
    tEvents = [0];
    plt.vline(tEvents,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEvents,ytext-2,'Go',StyleArgs{:});
    
    tEvents = [-4];
    plt.vline(tEvents,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEvents,ytext-2,'Cue',StyleArgs{:});
    
    tEventsCue = [-1.5];
    plt.vline(tEventsCue,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEventsCue,ytext-2,'Delay',StyleArgs{:});
    
end

%% plot all felt, observed, imagined, averaged across body parts

Attempt=[]; Imagine=[];
Attempt = [Att_hand+Att_shoulder]./2; Imagine = [Imag_hand+Imag_shoulder]./2;
PlotData = {{Attempt},{Imagine}};

StyleArgs = {'FontName','Arial','FontSize',9,'FontWeight','bold'};

figH = plt.fig('units','inches','width',14,'height',8,'font','Arial','fontsize',9);
figure(figH);
clf; pnl = panel(); pnl.margin=15; pnl.pack(3,4);

condNames = {'Attempted','Imagined'};

for condIdx=1:length(condNames)
    
    cData = PlotData{condIdx};
    
    clr=[];
    clr(1,:) = repmat([237 0 38]./255,1,1);
    clr(2,:) = repmat([31 120 180]./255,1,1);
    clr(3,:) = repmat([237 0 38]./255,1,1);
    clr(4,:) = repmat([1 108 89]./255,1,1);
    
    cLabel = condNames{condIdx};
    
    if condIdx==1 || condIdx==2
        
        TimeWindows = [-5:TimeStep:5];
        TimeCenters=[];
        for timeIdx = 1:length(TimeWindows)
            cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
            TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
        end
        
    elseif condIdx==3 || condIdx==4
        TimeWindows = [-5:TimeStep:5];
        TimeCenters=[];
        for timeIdx = 1:length(TimeWindows)
            cTW = [TimeWindows(timeIdx) TimeWindows(timeIdx)+StepSize];
            TimeCenters(timeIdx) = cTW(1)+(cTW(2)-cTW(1))/2;
        end
    end
    % plot the discriminability of touch from baseline
    pnl(1,condIdx).select()
    
    Analyze.plotEventRelatedAverage(cData,cLabel,'TimeVec',TimeCenters,'usePercentile',...
        'Colors',clr(condIdx,:));
    
    title(cLabel)
    ax1 = gca;
    set(ax1,'box','off',StyleArgs{:});
    ax1.XRuler.Axle.LineStyle = 'none';
    set(gca,'TickDir','out');
    set(ax1,'TickLength',[0 0]);
    axis tight;
    xlabel('Time (s)',StyleArgs{:});
    ylabel('From Baseline',StyleArgs{:});
    
    if condIdx==1 || condIdx==2
        xl = xlim;
        xlim([xl(1) 6])
    elseif condIdx==3 || condIdx==4
        xl = xlim;
        xlim([-2.8 5])
    end
end

for plotIdx = 1:length(condNames)
    
    pnl(1,plotIdx).select()
    ylimtmp1(plotIdx,:) = ylim();
    
end

for plotIdx = 1:length(condNames)
    
    pnl(1,plotIdx).select()
    ylim([min([ylimtmp1(:,1);0]) max(ylimtmp1(:,2))]);
    
    tEvents = [0];
    plt.vline(tEvents,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEvents,ytext-2,'Go',StyleArgs{:});
    
    tEvents = [-4];
    plt.vline(tEvents,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEvents,ytext-2,'Cue',StyleArgs{:});
    
    tEventsCue = [-1.5];
    plt.vline(tEventsCue,'k');
    ylims = ylim();
    ytext = ylims(2) + .05*diff(ylims);
    text(tEventsCue,ytext-2,'Delay',StyleArgs{:});
    
end

FigsDir=fullfile(BaseDir, 'PopData', 'Mine');
plt.SaveFigure(saveFig,FigsDir,['AttemptedImagined'],'SVG','PNG');
% plt.SaveFigure(saveFig,FigsDir,['AttemptedImagined'],'SVG');

end
