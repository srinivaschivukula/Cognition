function [LabelIDX,uniqY]=labels2IDX(Labels,labelOrder)

if nargin<2
 uniqY=unique(Labels);
else
   uniqY=labelOrder; 
end

LabelIDX=nan(size(Labels));

       for i=1:length(uniqY); 
       tmp=strcmp(Labels,uniqY{i});
       LabelIDX(tmp)=i;
    end