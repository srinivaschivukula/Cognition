function [UnitIDs,UnitLabels] = FindInterestingUnits(TaskLabel,varargin)
%ISINTERESTINGUNIT Returns list of "interesting" units.

[varargin,PThresh]=Utilities.ProcVarargin(varargin,'PThresh',.001);
[varargin,mode]=Utilities.ProcVarargin(varargin,'mode','Anova');

datafile = fullfile(env.get('results'),'FaceScratch3','SUAnal',[TaskLabel '-Go'],'PopData','p1-AllDays.mat');
load(datafile); % Loads PlotData
Coefs = horzcat(PlotData{1}.Coef{:})';
PVals = horzcat(PlotData{1}.Pvals{:})';
IsSig = PVals < PThresh;
ASF = cell2mat(vertcat(Analyze.returnFieldValues(PlotData,'Unit')));

if strcmp(mode,'Anova')
    for m = 1:length(PlotData{1}.Model_p{1});
        MdlP{m} = [];
        for i = 1:length(PlotData{1}.Model_p)
            MdlP{m} = [MdlP{m} PlotData{1}.Model_p{i}{m}];
        end
        MdlP{m} = MdlP{m}';
        MdlSig{m} = MdlP{m} < PThresh;
    end
    
    switch TaskLabel
        case {'XBodySpec','XBodySpecC'}
            % Tuned, BPSpec,PersonSpec,Both
            NC = 1; NS = 3; TC = 4; TS = 6;
            idx = [NC NS TC TS];
            
            %         Tuned = any(IsSig(:,idx),2);
            BOnly = ~MdlSig{1}(:,1) & MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
            POnly = MdlSig{1}(:,1) & ~MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
            BandP = MdlSig{1}(:,1) & MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
            Intxn = MdlSig{1}(:,3);
            
            good = BOnly | POnly | BandP | Intxn;
            
            UnitIDs = ASF(good,:);
            UnitLabels = cell(size(ASF,1),1); % Empty by default
            for i = 1:length(UnitLabels)
                if BOnly(i)
                    UnitLabels{i} = 'BodyPart';
                elseif POnly(i)
                    UnitLabels{i} = 'Person';
                elseif BandP(i)
                    UnitLabels{i} = 'Both';
                else
                    UnitLabels{i} = 'Intxn';
                end
                % empty if none of these cases.
            end
            UnitLabels = UnitLabels(good);
            
            %         ChSpec = IsSig(:,NC) & IsSig(:,TC) & ~IsSig(:,NS) & ~IsSig(:,TS);
            %         ShSpec = ~IsSig(:,NC) & ~IsSig(:,TC) & IsSig(:,NS) & IsSig(:,TS);
            
            %         BPSpec = ChSpec | ShSpec;
            %         MdlBPSpec = ~MdlSig{1}(:,1) & MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
            
            %         BPSpecIdx = find(BPSpec);
            %         MdlBPSpecIdx = find(MdlBPSpec);
            
            %         size(intersect(BPSpecIdx,MdlBPSpecIdx))
            %         size(setdiff(BPSpecIdx,MdlBPSpecIdx))
            %         size(setdiff(MdlBPSpecIdx,BPSpecIdx))
            
            %         LabelTypes = {'Lin1','Lin2','Both'};
            %         good = BPSpec | MdlBPSpec;
            %         UnitLabels = LabelTypes(BPSpec(good) + MdlBPSpec(good)*2);
            %         UnitIDs = ASF(good,:);
            
            %         UnitIDs = ASF(MdlBPSpec,:);
            %         UnitLabels = {};
            
        case 'XFixation'
            % BP Spec units that show mirroring only with free gaze
            BPSpec = MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
            
            FixSig = MdlSig{2}(:,1);
            
            
            LabelTypes = {'FixInv','FixTuned'};
            good = BPSpec;
            UnitLabels = LabelTypes(FixSig(good)+1);
            UnitIDs = ASF(good,:);
        case 'XFSTL'
            % BPSpec and modality??
        case 'XActSense'
            % BP Spec and invariant to self vs other doing the action
            BPSpec = MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
            
            ActorSig = MdlSig{2}(:,1);
            
            LabelTypes = {'ActorInv','ActorTuned'};
            good = BPSpec;
            UnitLabels = LabelTypes(ActorSig(good)+1);
            UnitIDs = ASF(good,:);
        case 'XVideoLiveXBodySpecC'
            % BP Spec for live vs video
            MdlBPSpecLive = MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
            MdlBPSpecVideo = MdlSig{3}(:,2) & ~MdlSig{3}(:,3);
            
            LabelTypes = {'Live','Video','Both'};
            good = MdlBPSpecLive | MdlBPSpecVideo;
            UnitLabels = LabelTypes(MdlBPSpecLive(good) + MdlBPSpecVideo(good)*2);
            UnitIDs = ASF(good,:);
        case 'XActSenseLive'
            BPSpec = MdlSig{1}(:,2) & ~MdlSig{1}(:,3);
            IsSig = any(PVals < PThresh,2);
            
            LabelTypes = {'Sig','BPSpecOnly','BPSpec'};
            good = BPSpec | IsSig;
            UnitLabels = LabelTypes(IsSig(good) + BPSpec(good)*2);
            UnitIDs = ASF(IsSig | BPSpec,:);
    end
elseif strcmp(mode,'ModelCompare')
    BIC = [];
    MdlP = [];
    for i = 1:length(PlotData{1}.ModelCompare_BIC)
        BIC = [BIC; PlotData{1}.ModelCompare_BIC{i}{1}];
        MdlP = [MdlP; PlotData{1}.ModelCompare_PVal{i}{1}];
    end
    [~,BestModelIdx] = min(BIC,[],2);
    BPSpec = (BestModelIdx==3) & any(MdlP<PThresh,2);
    
    if strcmp(TaskLabel,'XBodySpec')
        UnitIDs = ASF(BPSpec,:);
        UnitLabels = [];
        return;
    end
    
    BIC = [];
    for i = 1:length(PlotData{1}.ModelCompare_BIC)
        BIC = [BIC; PlotData{1}.ModelCompare_BIC{i}{2}];
    end
    if strcmp(TaskLabel,'XFixation')
        models2use = [5 7];
        LabelTypes = {'FixInv','Other'};
    elseif strcmp(TaskLabel,'XVideoLiveXBodySpecC')
        models2use = [3 5];
        LabelTypes = {'FormatInv','Other'};
    elseif strcmp(TaskLabel,'XActSense')
        models2use = [1 2];
        LabelTypes = {'ActorInv','Other'};
    else
        models2use = 1:size(BIC,2);
    end
    BIC = BIC(BPSpec,models2use);
    [~,BestModelIdx] = min(BIC,[],2);
%     invIdx = BestModelIdx==1;
    
    
    UnitLabels = LabelTypes(BestModelIdx);
    UnitIDs = ASF(BPSpec,:);
end


end

