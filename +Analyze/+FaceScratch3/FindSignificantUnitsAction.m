function UnitIDs = FindSignificantUnitsAction(TaskLabel,varargin)
%FindSignificantUnits Returns list of significant units, by pval

[varargin,PThresh]=Utilities.ProcVarargin(varargin,'PThresh',.001);
[varargin,Labels]=Utilities.ProcVarargin(varargin,'Labels',{});
[varargin,Labels2Incl]=Utilities.ProcVarargin(varargin,'Labels2Incl',{});
[varargin,UnitGroup]=Utilities.ProcVarargin(varargin,'UnitGroup','All');
[varargin,Dates]=Utilities.ProcVarargin(varargin,'Dates',[]);

datafile = fullfile(env.get('results'),'FaceScratch3','SUAnal',[TaskLabel '-Go'],'PopData','p1-AllDays.mat');
load(datafile); % Loads PlotData
ASF = cell2mat(vertcat(Analyze.returnFieldValues(PlotData,'Unit')));
if strcmp(UnitGroup,'All')
    PVals = horzcat(PlotData{1}.Pvals{:})';
%     IsSig = PVals < PThresh;
    IsSig = Utilities.MultipleComparisonsCorrection(PVals,'method','fdr');

    IsSig = PVals<1;
    if isempty(Labels) | isempty(Labels2Incl)
        idx2incl = 1:length(Labels);
    else
        idx2Incl = find(ismember(Labels,Labels2Incl));
    end
    GoodIdx = any(IsSig(:,idx2Incl),2);
    UnitIDs = ASF(GoodIdx,:);
    InclIdx = GoodIdx;
else
    MdlP = [];
    for i = 1:length(PlotData{1}.Model_p)
        MdlP = [MdlP PlotData{1}.Model_p{i}{1}];
    end
    MdlP = MdlP';
%     isSig = MdlP < PThresh;
    isSig = Utilities.MultipleComparisonsCorrection(MdlP,'method','fdr');
    
    % Terms are: Person, BP, Act, P*BP, P*Act, BP*Act
    if strcmp(UnitGroup,'PBPSpec')
        PBP = isSig(:,1) | isSig(:,2);
        Intxn = any(isSig(:,4:6),2);
        InclIdx = PBP & ~Intxn;
    elseif strcmp(UnitGroup,'ActSpec')
        Act = isSig(:,3);
        Intxn = any(isSig(:,4:6),2);
        InclIdx = Act & ~Intxn;
    elseif strcmp(UnitGroup,'PBPSpecOnly')
        PBP = (isSig(:,1) | isSig(:,2)) & ~isSig(:,3);
        Intxn = any(isSig(:,4:6),2);
        InclIdx = PBP & ~Intxn;
    elseif strcmp(UnitGroup,'ActSpecOnly')
        Act = ~isSig(:,1) & ~isSig(:,2) & isSig(:,3);
        Intxn = any(isSig(:,4:6),2);
        InclIdx = Act & ~Intxn;
    elseif strcmp(UnitGroup,'NotActSpec')
        Act = isSig(:,3);
        InclIdx = ~Act;
    elseif strcmp(UnitGroup,'NotActObs')
        PVals = horzcat(PlotData{1}.Pvals{:})';
%         isSig = PVals < PThresh;
        isSig = Utilities.MultipleComparisonsCorrection(PVals,'method','fdr');
        ActObs = ~any(isSig(:,1:8),2) & any(isSig(:,9:16),2);
        InclIdx = ~ActObs;
    elseif strcmp(UnitGroup,'sigUnits')
        InclIdx = any(isSig,2);
    end
    
%     if size(MdlP,2)==3
%         % 3 terms means terms are: [person, bodypart, person:bodypart]
%         PureB=(isSig(:,1)==0 & isSig(:,2)==1 & isSig(:,3)==0); % Pure BP
%         AandB=(isSig(:,1)==1 & isSig(:,2)==1 & isSig(:,3)==0); % Pure BP + Person, no itneraction
%         InclIdx = PureB | AandB;
%     else
%         error('Wrong number of terms?');
% %         if size(MdlP,2)==5
% %         % means terms are: [person, bodypart1, bodypart2, person:bodypart1, person:bodypart2]
% %         isSig(:,[2 3]) = MdlP(:,[2 3]) < PThresh/2;
% %         PureB=(isSig(:,1)==0 & isSig(:,2)==1 & isSig(:,3)==0); % Pure BP
% %         AandB=(isSig(:,1)==1 & isSig(:,2)==1 & isSig(:,3)==0); % Pure BP + Person, no itneraction
% %         InclIdx = PureB | AandB;
%     end
    
    UnitIDs = ASF(InclIdx,:);
end


if ~isempty(Dates)
    unitDates = Blackrock.Helper.unitId2date(UnitIDs);
    goodIdx = ismember(unitDates,Dates);
    UnitIDs = UnitIDs(goodIdx,:);
end

fprintf('Including %d/%d units\n',nnz(InclIdx),size(ASF,1));

end