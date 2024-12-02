% Based on Tyson/MLfw/+tests/FindCenterOutFiles.m
N1=hst.Subject('n1')

today = datestr(now, 'yyyymmdd');
o=N1.getSessions('RANGE',{'20190401',today});
AllOpenLoop=table()
for sesID=1:length(o)
    cSes=N1.getSessionObject(o{sesID});
    
    cSes.TaskInfo=cSes.TaskInfo(cSes.TaskInfo.NumTrials>=48 & strcmp(cSes.TaskInfo.TaskConfig,'Parameters_CenterOut'),:);
    
    for taskID=1:size(cSes.TaskInfo,1)
%         try
        task=hst.Task(cSes.TaskInfo.TaskFile{taskID});
        
        t=cSes.TaskInfo(taskID,[2 4 6 9 10]);
        t.Date=task.idString(1:8);
        tmp=strsplit(task.taskString,'-');
        t.Time= tmp{3};
        
        if isfield(task.predictor,'hDecoderNDF')
            t.OpenLoop=task.predictor.hDecoderNDF.Params.assistLevel;
        elseif isfield(task.predictor,'DecoderIDX')
                    t.OpenLoop=task.predictor.hDecoder{task.predictor.DecoderIDX}.Params.assistLevel;
        elseif isfield(task.predictor,'hDecoder')
            t.OpenLoop=task.predictor.hDecoder.runtimeParams.assistLevel;
        end
        
        if t.OpenLoop~=1
            continue
        end
       
        AllOpenLoop=[AllOpenLoop;t];
%         catch
%             lasterr
%         end
    end
end
    % TaskFiles=cSes.TaskInfo.TaskFile 
