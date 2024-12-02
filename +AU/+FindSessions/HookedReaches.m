function [filteredTrialData]=HookedReaches(varargin)
%% HookedReaches: Find Center-out data with the largest hook at the
% beginning.
% Focus on first 500ms, before subject has time to compensate.
% Note: throughout each session's day, subject tends to compensate over
% time.

% Goal is to find some of these reaches to display for Grant 2020
% application

%% Grab task files
subject = env.get('subject');
% DateVec = {{'20190601', '20200101'}};
% [finalList,dates,RunData]=getTaskFiles('ParamFile', 'Task.UnityCursorControl.Parameters_CenterStart',...
% [finalList,dates,RunData]=AU.getTaskFiles(...
%     'ParamFile', 'Task.UnityCursorControl.Parameters_CenterOut',...
%     'Subject', subject,...
%     'TaskWild', '*UnityCursorControl*.mat',...
%     'DateVec', DateVec,...
%     varargin{:});

%% Plotting parameters
targetDist = 0.4;


%%
% TaskName = 'CenterOut'; Dates = {'20190718', '20190122'};
% TaskName = 'CenterOut'; Dates = {'20190910', '20190820', '20190718', '20190122', '20190507'};
TaskName = 'CenterOut'; Dates = {'20190122', '20190125', '20190215', '20190314', '20190402', '20190419', '20190813', '20190827', '20191008'};
% TaskName = 'NeuralDynamics'; Dates = {'20190412', '20190517', '20190528'};
% TaskName = 'SpeedTracking'; Dates = {'20190719', '20190726', '20190806', '20190813'};
AllTrialData = Analyze.LoadConvertedData(TaskName, Dates);

%% Try this since server is down
% AllTrialData = Analyze.ProcessRawData(TaskName,...
%     'TaskFile', '20191029-133514-142031-UnityCursorControl',...
%     'TaskConfig', 'Analyze.CenterOut.Config',...
%     'Run', 1,...
%     'TaskLabel', 'test label',...
%     'NoNeural',...
%     'overwrite', true);
AllTrialData=Analyze.CenterOut.CenterOut_ComputeTargetAngles(AllTrialData);
AllTrialData=Analyze.CenterOut.CenterOut_AppendContinuousStats(AllTrialData);

% Double check that this is 
TaskName
% AllTrialData{1}.Info.Task.parameterName
AllTrialData{1}.Info.Options.taskConfig

% Or load data, then do Task.params.parameterFcn
% Data.endComment

%% Filter for outward reaches
filteredTrialData = Analyze.SubSelectTrials(AllTrialData, 'ReturnHome', false, 'Phase', 'Go');

%% Look for curved reaches
windowOffset = [0.1, -0.1];
windowType = 'TrialDuration';
if true
for sessionIdx = 1:numel(filteredTrialData)
    Analyze.CenterOut.PlotBehavior(filteredTrialData(sessionIdx), 'windowOffset', windowOffset, 'WindowType', windowType);
    try
        suptitle(filteredTrialData{sessionIdx}.Info.TaskID);
    catch
        warning('skipping title');
    end
end
return
end

%% Find reaches with the highest angular error in the reaction/reach-start-up time
if false
tWindow = [0.1, 0.5];

for sessionIdx = 1:numel(filteredTrialData)
    cTrialData = filteredTrialData{sessionIdx};
    if ~isempty(cTrialData)
        angularError = Analyze.getContinuousData2(cTrialData, tWindow, 'AngularError', 'takeMean');
        state = Analyze.getContinuousData2(cTrialData, tWindow(end), 'state');
        stateMatrix = cell2mat(state.');
        x = stateMatrix(:,1);
        dx = stateMatrix(:,2);
        y = stateMatrix(:,3);
        dy = stateMatrix(:,4);
        
        figure(sessionIdx);
        quiver(x, y, dx, dy);
        xlim([-targetDist, targetDist]);
        ylim([-targetDist, targetDist]);
        axis('equal')
        
        title({cTrialData.TaskLabel{1},...
               cTrialData.Info.Task.userEndComment,...
               sprintf('Average angular error: %0.1f', nanmean(cell2mat(angularError)))});
    end
end
end

%%


%% Process each session separately?
% Other option: measure the orthogonal distance at time
% Time at which we are measuring the orthogonal distance from trajectory
t = 0.5;

% Or, filter for highest angular error in first targetDist after seeing
% target

for sessionIdx = 1:numel(filteredTrialData)
    cTrialData = filteredTrialData{sessionIdx};
    if ~isempty(cTrialData)
        orthogDist = Analyze.getContinuousData2(cTrialData, t, 'OrthogDistTemp');
        state = Analyze.getContinuousData2(cTrialData, t, 'state');
        stateMatrix = cell2mat(state.');
        x = stateMatrix(:,1);
        y = stateMatrix(:,3);
        
        figure(sessionIdx);
        scatter(x, y);
        xlim([-targetDist, targetDist]);
        ylim([-targetDist, targetDist]);
        axis('equal')
        
        title({cTrialData.TaskLabel{1},...
            cTrialData.Info.Task.userEndComment,...
            sprintf('orthogDist: %0.3f', nanmean(cell2mat(orthogDist)))});
    end
end