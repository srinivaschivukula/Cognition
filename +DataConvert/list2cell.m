function [cData,cLabel,LabelList]=list2cell(Data,Label)
% Turn a list with labels into a cell array, grouped by label
% Note : cell2list does the reverse.

% Handle case of numerical list
% if ~iscell(Data)
if ndims(Data)==2
    [Labels,ia,ib]=unique(Label);
    for i=1:length(Labels)
        
        cData{i}=Data(ib==i,:);
         LabelList{i}=Label(ib==i,:);
        if iscell(Labels)
            cLabel{i}=Labels{i};
        else
            cLabel{i}=Labels(i);
        end
    end
else
    
    [Labels,ia,ib]=unique(Label);
    for i=1:length(Labels)
        
        cData{i}=Data(ib==i,:,:);
          LabelList{i}=Label(ib==i,:);
        if iscell(Labels)
            cLabel{i}=Labels{i};
        else
            cLabel{i}=Labels(i);
        end
    end
    
end

% else
%     for i=1:length(Label)
%         INDS=strcmp(Data,Label{i});
%         INDSStore(:,i)=INDS;
%     end

% end
