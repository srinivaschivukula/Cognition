function [uniqY,LabelName,N]=getLabelAssociations(Y,Labels)

 uniqY=unique(Y);
    
    for i=1:length(uniqY); 
       tmp=find(Y==uniqY(i));
       N(i)=length(find(Y==uniqY(i)));
       tmp=tmp(1);
       LabelName{i}=Labels{tmp};
       
    end