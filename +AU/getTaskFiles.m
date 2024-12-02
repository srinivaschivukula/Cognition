function [finalList,dates,RunData]=getTaskFiles(varargin)


ParamFile='';
TaskWild='*.mat';
Subject='p1';
BaseDirs=env.get('data');
DateVec={{'20170303','20180314'}};
% ParamFil 
[varargin,ParamFile]=Utilities.ProcVarargin(varargin,'ParamFile',ParamFile);
[varargin,TaskWild]=Utilities.ProcVarargin(varargin,'TaskWild',TaskWild);
[varargin,Subject]=Utilities.ProcVarargin(varargin,'Subject',Subject);
[varargin,BaseDirs]=Utilities.ProcVarargin(varargin,'BaseDirs',BaseDirs);
[varargin,DateVec]=Utilities.ProcVarargin(varargin,'DateVec',DateVec);
 
% dates = hst.getSessionDates(Subject,BaseDirs,DateVec{1},DateVec{2});
% dates = Analyze.getSessionDates(Subject,'Y:\',DateVec{1},DateVec{2});
dates = Analyze.getSessionDates(Subject,'Y:\data',DateVec);


for i=1:length(dates)
    fileList{i} = Analyze.getSessionFiles(dates{i},TaskWild,Subject,BaseDirs);
end

dates(cellfun(@isempty,fileList))=[];
fileList(cellfun(@isempty,fileList))=[];

for i=1:length(fileList)
    isValid=~cellfun(@isempty,strfind(fileList{i},'Task'));
    fileList{i}= fileList{i}(isValid);
end
%%

for i=1:length(fileList)
    tmpFiles=fileList{i};
    clear isMatch EndCommentTMP EndTimeTMP NTrialsTmp
    if ~isempty(tmpFiles)
        for j=1:length(tmpFiles)
            Block=load(tmpFiles{j});
            if ~isfield(Block,'saveData')
                if ~isempty(ParamFile)
                if iscell(Block.Options.taskConfig)
                    isMatch(j)=~isempty(strfind(Block.Options.taskConfig{1},ParamFile));
                else
                    isMatch(j)=~isempty(strfind(Block.Options.taskConfig,ParamFile));
                end
                else
                    isMatch(j)=1;
                end
                try
                EndCommentTMP{j}=Block.Data.comments.Comment{end};
                catch
                     try
                EndCommentTMP{j}=Block.saveData.endComment;
                catch
                    EndCommentTMP{j}='NoEndComment';
                end
                end
                EndTimeTMP(j)=Block.Options.timerPeriod*Block.Task.frameId;
                NTrialsTmp(j)=length(Block.Task.TrialData);
            else
                isMatch(j)=logical(1);
                try
                EndCommentTMP{j}=Block.saveData.endComment;
                catch
                  try
                EndCommentTMP{j}=Block.Data.comments{end}{end};
                catch
                    EndCommentTMP{j}='NoEndComment';
                end
                end
                
                
                EndTimeTMP(j)=Block.saveData.RunEnd_secs;
                NTrialsTmp(j)=length(Block.saveData.Trials);
            end
            
        end
    finalList{i}=fileList{i}(isMatch);
    EndComment{i}=EndCommentTMP(isMatch);
    EndTime{i}=EndTimeTMP(isMatch);
    NTrials{i}=NTrialsTmp(isMatch);
    end
end


noFiles=cellfun(@isempty,finalList);


EndComment(noFiles)=[];
EndTime(noFiles)=[];
NTrials(noFiles)=[];
dates(noFiles)=[];
finalList(noFiles)=[];

RunData.EndComment=EndComment;
RunData.EndTime=EndTime;
RunData.NTrials=NTrials;

outfile = ['Task_log.txt'];
fileID = fopen(outfile,'w');

out={};
idx = 1;
out{idx,1} = sprintf('~~~%s~~~~~~~~~~~~~~~',Subject);idx = idx+1;

for d = 1:length(dates)

out{idx,1} = sprintf('~~~%s~~~~~~~~~~~~~~~',dates{d});idx = idx+1;
 out{idx,1} = ''; idx = idx+1;

    for f = 1:length(finalList{d})
        if isempty(finalList{d}{f})
            continue;
        end
        out{idx,1} = sprintf('  %s',finalList{d}{f}); idx = idx+1;
        out{idx,1} = sprintf('      NTrials= %d; Time = %0.2f',RunData.NTrials{d}(f),RunData.EndTime{d}(f));idx = idx+1;
        out{idx,1} = sprintf('      %s',RunData.EndComment{d}{f});idx = idx+1;
        out{idx,1} = ''; idx = idx+1;
    end
    out{idx,1} = ''; idx = idx+1;
end

Utilities.writecell(outfile,out);

% fclose(fileID);
winopen(outfile);
%%