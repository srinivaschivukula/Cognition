function Options=Config_text(varargin)


[varargin,IncludeTouchSensor] = Utilities.ProcVarargin(varargin,'IncludeTouchSensor');
[varargin,TaskLabel]   = Utilities.ProcVarargin(varargin,'TaskLabel',{});

% Config file that determines how to process the Symbolic Task
Options.TaskType='experiment2';
Options.TaskVersion='DateTimeTime';

% Info

if strcmpi(TaskLabel,'XMotor')
    Options.phases2Process={'ITI','CueTarget','Delay','Go'};  % Phases to Analyze
else
    Options.phases2Process={'CueTarget','Delay','Go'};  % Phases to Analyze
end

Options.useSorted=1; % use crossings or single units

% Options.NSP = 1:2;  % NSPs to analyze
% Options.activeChannels = {[1:96]',[1:96]'};  % channels to analyze
Options.NSP = {'NSP1'};  % NSPs to analyze\
if IncludeTouchSensor
    fprintf('\n\nIncluding Touch Sensor!');
    Options.loadRawAInp.loadData=1;
    Options.loadRawAInp.nsp='NSP1';
    Options.loadRawAInp.channelLabel={'ainp2'};
    Options.loadRawAInp.Fs = 2000; % Configured in central
end
Options.activeChannels = {[1:96]'};  % channels to analyze

Options.saveDataCBmexUnits='seconds';

if strcmpi(TaskLabel,'XMotor')
    Options.TrialDataConstructor=@Analyze.FaceScratch3.TrialDataConstructorMotor;
else
    Options.TrialDataConstructor=@Analyze.FaceScratch3.TrialDataConstructorText;
end
Options.ResultsDir=fullfile(env.get('result'), 'FaceScratch3')


