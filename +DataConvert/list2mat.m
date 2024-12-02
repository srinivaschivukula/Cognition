function [cData,cLabel]=list2mat(Data,Label)
% Turn a cell array into a list with labels
% Note : list2cell does the reverse.


[Labels,ia,ib]=unique(Label);
for i=1:length(Labels)
    
    cData(:,i)=Data(ib==i);
    if iscell(Labels)
        cLabel{i}=Labels{i};
    else
        cLabel{i}=num2str(Labels(i));
    end
end
