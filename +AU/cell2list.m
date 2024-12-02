function [Data,Inds]=cell2list(CellData,dim)
% Turn a cell array into a list with labels
% Note : list2cell does the reverse. 

if nargin<2
    dim=1;
end
%%
Data=cat(dim,CellData{:});

N=cellfun(@(x)size(x,dim),CellData);

Inds=[];
for i=1:length(CellData)
    Inds=[Inds;ones(N(i),1)*i];
end


end
