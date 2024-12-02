function UnitIDs = FindSignificantUnits(TaskLabel,varargin)
%FindSignificantUnits Returns list of significant units, by pval

% if only mirror units desired, include the following in the varargin in
% the function call

% ('Mirror','AllTrialData',AllTrialData);



[varargin,PThresh]=Utilities.ProcVarargin(varargin,'PThresh','fdr');
[varargin,Labels]=Utilities.ProcVarargin(varargin,'Labels',{});
[varargin,Labels2Incl]=Utilities.ProcVarargin(varargin,'Labels2Incl',{});
[varargin,BPSpecOnly]=Utilities.ProcVarargin(varargin,'BPSpecOnly');
[varargin,PSpecOnly]=Utilities.ProcVarargin(varargin,'PSpecOnly');
[varargin,Mirror]=Utilities.ProcVarargin(varargin,'Mirror');
[varargin,Dates]=Utilities.ProcVarargin(varargin,'Dates',[]);
[varargin,All]=Utilities.ProcVarargin(varargin,'All');
[varargin,trialdata]=Utilities.ProcVarargin(varargin,'AllTrialData',{});

datafile = fullfile(env.get('results'),'FaceScratch3','SUAnal',[TaskLabel '-Go'],'PopData','p1-AllDays.mat');
load(datafile); % Loads PlotData
ASF = cell2mat(vertcat(Analyze.returnFieldValues(PlotData,'Unit')));
if ~BPSpecOnly & ~PSpecOnly & ~Mirror & ~All
    PVals = horzcat(PlotData{1}.Pvals{:})';
    if strcmpi(PThresh,'fdr') || strcmpi(PThresh,'holm-bonferroni')
        [IsSig,alpha] = Utilities.MultipleComparisonsCorrection(PVals,'method',PThresh);
    else
        IsSig = PVals < 0.001;
    end
    %     IsSig = PVals<1;
    %
    %     MdlP = [];
    %     for i = 1:length(PlotData{1}.Model_p)
    %         MdlP = [MdlP PlotData{1}.Model_p{i}{1}];
    %     end
    %     MdlP = MdlP';
    %     [IsSig,alpha2] = Utilities.MultipleComparisonsCorrection(MdlP,'method','fdr');
    
    
    if isempty(Labels) || isempty(Labels2Incl)
        idx2incl = 1:length(Labels);
    else
        idx2Incl = find(ismember(Labels,Labels2Incl));
    end
    GoodIdx = any(IsSig(:,idx2Incl),2);
    %     GoodIdx = any(IsSig,2);
    UnitIDs = ASF(GoodIdx,:);
elseif All
    UnitIDs = vertcat(PlotData{1}.Unit{:});
    
elseif Mirror
    BIC = [];
    BIC = vertcat(PlotData{1}.ModelCompare_BIC{:});
    
    MdlP = [];
    for i = 1:length(PlotData{1}.Model_p)
        % First model should be 2x2 anova figure out body spec
        MdlP = [MdlP PlotData{1}.Model_p{i}{1}];
    end
    MdlP = MdlP';
    [~,effAlpha] = Utilities.MultipleComparisonsCorrection(MdlP,'method',PThresh);
    isSig = MdlP < effAlpha;
    PureA = (isSig(:,1)==1 & isSig(:,2)==0 & isSig(:,3)==0); % Pure Action
    PureB=(isSig(:,1)==0 & isSig(:,2)==1 & isSig(:,3)==0); % Pure Action
    AandB=(isSig(:,1)==1 & isSig(:,2)==1 & isSig(:,3)==0); % Pure Action
    %         GoodIdx = PureB | AandB;
    GoodIdx = PureB;
    
    models2use = logical([0 0 1 1 1 0 0 0]);
    
    BIC = BIC(:,models2use);
    [~,BestModelIdx] = min(BIC,[],2);
    
    NUnits = length(BestModelIdx(GoodIdx==3));
    cUnits = ASF(GoodIdx(BestModelIdx==3),:);
    TimeWindow = [1 4];
    
    testPvals = [];
    testPvals = vertcat(PlotData{1}.CheekShoulderPval{:});
    testPvals = vertcat(testPvals{:});
    
    order =[];
    order = vertcat(PlotData{1}.CheekShoulderOrder{:});
    order = vertcat(order{:});
    
    unitIDX = (GoodIdx & BestModelIdx==3);
    testPvals = testPvals(unitIDX);
    order = order(unitIDX);
    
    Selective = [~(testPvals<0.05);zeros(length(ASF)-length(~(testPvals<0.05)),1)];
    
    mirrorUnits = logical(zeros(length(ASF),1));
    mirrorUnits(ismember(ASF,cUnits,'rows')) = ~(testPvals<0.05);
    
    UnitIDs = ASF(mirrorUnits,:);
else
    MdlP = [];
    for i = 1:length(PlotData{1}.Model_p)
        MdlP = [MdlP PlotData{1}.Model_p{i}{1}];
    end
    MdlP = MdlP';
    
    if strcmpi(PThresh,'fdr') || strcmpi(PThresh,'holm-bonferroni')
        isSig = Utilities.MultipleComparisonsCorrection(MdlP,'method',PThresh);
    else
        isSig = MdlP < 0.001;
    end
    %     isSig = MdlP < PThresh;
    
    if size(MdlP,2)==3
        % 3 terms means terms are: [person, bodypart, person:bodypart]
        PureA=(isSig(:,1)==1 & isSig(:,2)==0 & isSig(:,3)==0); % Pure P
        PureB=(isSig(:,1)==0 & isSig(:,2)==1 & isSig(:,3)==0); % Pure BP
        AandB=(isSig(:,1)==1 & isSig(:,2)==1 & isSig(:,3)==0); % Pure BP + Person, no itneraction
        if BPSpecOnly
            InclIdx = PureB | AandB;
        elseif PSpecOnly
            InclIdx = PureA | AandB;
        end
    else
        error('Wrong number of terms?');
        %         if size(MdlP,2)==5
        %         % means terms are: [person, bodypart1, bodypart2, person:bodypart1, person:bodypart2]
        %         isSig(:,[2 3]) = MdlP(:,[2 3]) < PThresh/2;
        %         PureB=(isSig(:,1)==0 & isSig(:,2)==1 & isSig(:,3)==0); % Pure BP
        %         AandB=(isSig(:,1)==1 & isSig(:,2)==1 & isSig(:,3)==0); % Pure BP + Person, no itneraction
        %         InclIdx = PureB | AandB;
    end
    
    UnitIDs = ASF(InclIdx,:);
end


if ~isempty(Dates)
    unitDates = Blackrock.Helper.unitId2date(UnitIDs);
    goodIdx = ismember(unitDates,Dates);
    UnitIDs = UnitIDs(goodIdx,:);
end

end