function [UnitIDs,UnitLabels] = FindInterestingUnitsAction(TaskLabel,varargin)
%ISINTERESTINGUNIT Returns list of "interesting" units.

[varargin,PThresh]=Utilities.ProcVarargin(varargin,'PThresh',.001);

datafile = fullfile(env.get('results'),'FaceScratch','SUAnal',[TaskLabel '-Go'],'PopData','p1-AllDays.mat');
load(datafile); % Loads PlotData
Coefs = horzcat(PlotData{1}.Coef{:})';
PVals = horzcat(PlotData{1}.Pvals{:})';
IsSig = PVals < PThresh;
[IsSig,effa] = Utilities.MultipleComparisonsCorrection(PVals,'method','fdr');
effa
ASF = cell2mat(vertcat(Analyze.returnFieldValues(PlotData,'Unit')));

for m = 1:length(PlotData{1}.Model_p{1})
    MdlP{m} = [];
    for i = 1:length(PlotData{1}.Model_p)
        MdlP{m} = [MdlP{m} PlotData{1}.Model_p{i}{m}];
    end
    MdlP{m} = MdlP{m}';
    MdlSig{m} = MdlP{m} < PThresh;
end

% Terms are: Person, BP, Act, P*BP, P*Act, BP*Act

BPOnly = ismember(MdlSig{1},[0 1 0 0 0 0],'rows');
POnly = ismember(MdlSig{1},[1 0 0 0 0 0],'rows');
ActOnly = ismember(MdlSig{1},[0 0 1 0 0 0],'rows');

code = [1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];
PiC = ismember(IsSig,code,'rows');

code = [0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
PrC = ismember(IsSig,code,'rows');

code = [0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
RbC = ismember(IsSig,code,'rows');

code = [0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0];
TpC = ismember(IsSig,code,'rows');

code = [0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0];
PiS = ismember(IsSig,code,'rows');

code = [0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0];
PrS = ismember(IsSig,code,'rows');

code = [0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0];
RbS = ismember(IsSig,code,'rows');

code = [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0];
TpS = ismember(IsSig,code,'rows');

good = BPOnly | POnly | ActOnly;
UnitIDs = ASF(good,:);
UnitLabels = cell(size(ASF,1),1);
UnitLabels(BPOnly) = {'BodyPart'};
UnitLabels(POnly) = {'Person'};
UnitLabels(ActOnly) = {'Action'};
UnitLabels(PiC) = {'PiC'};
UnitLabels(PrC) = {'PrC'};
UnitLabels(RbC) = {'RbC'};
UnitLabels(TpC) = {'TpC'};
UnitLabels(PiS) = {'PiS'};
UnitLabels(PrS) = {'PrS'};
UnitLabels(RbS) = {'RbS'};
UnitLabels(TpS) = {'TpS'};
UnitLabels = UnitLabels(good);

%     switch TaskLabel
%         case {'XBodySpec','XBodySpecC'}
%             % Tuned, BPSpec,PersonSpec,Both
%             NC = 1; NS = 3; TC = 4; TS = 6;
%             idx = [NC NS TC TS];
%             
%             %         Tuned = any(IsSig(:,idx),2);
%             BOnly = ~MdlSig{1}(:,1) & MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
%             POnly = MdlSig{1}(:,1) & ~MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
%             BandP = MdlSig{1}(:,1) & MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
%             %         Intxn = MdlSig{1}(:,3);
%             
%             good = BOnly | POnly | BandP;
%             
%             UnitIDs = ASF(good,:);
%             UnitLabels = cell(size(ASF,1),1); % Empty by default
%             for i = 1:length(UnitLabels)
%                 if BOnly(i)
%                     UnitLabels{i} = 'BodyPart';
%                 elseif POnly(i)
%                     UnitLabels{i} = 'Person';
%                 elseif BandP(i)
%                     UnitLabels{i} = 'Both';
%                 end
%                 % empty if none of these cases.
%             end
%             UnitLabels = UnitLabels(good);


end

