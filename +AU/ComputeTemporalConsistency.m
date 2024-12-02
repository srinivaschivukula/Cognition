function [CC]=ComputeTemporalConsistency(Data,Labels,NumTests)

% Test to see if one condition (column) is reliably the larger then all
% other columns across random 50/5-0 NumTests in the data

% Data : Number of samples X number of conditions. 
% labels Optional condition names
% NumTests  Number of NumTests
%
% Reliability : fraction of time splits result in same optimal condition.
% PrefCond : The condition that is most relaible
% Counts : Number of times each condition came out as preferred
%%
if nargin<2 || isempty(Labels); Labels=1:size(Data,2); end
if nargin<3; NumTests=500; end

N=size(Data,1);

% If the amount of data is reasonable, explicitly generate all possible
% combinations, and draw from the list (otherwise just sample randomly)
if N<=20
    useExplicit=1;
    C = nchoosek(1:N,ceil(N/2));
    C=C(randperm(size(C,1)),:);
    % can't do more NumTests then unique combinations.
    if size(C,1)<NumTests 
        NumTests=size(C,1);
    end
else
    useExplicit=0;
end

%%

for i=1:NumTests
    if useExplicit
        TrainInds=C(i,:); % grab 50% of the data
    else
        TrainInds=randperm(N,N/2); % grab 50% of the data
    end
    TestInds=setdiff(1:N,TrainInds);  % grab remaining data data
    
    train=mean(Data(TrainInds,:),1);
    test=mean(Data(TestInds,:),1);
    
    CC(i)=corr(train',test');
    
%     RAW(i).train=train;
%     RAW(i).test=test;
%     RAW(i).trainValsRAW=trainValsRAW(i,:);
%     RAW(i).testValsRAW=testValsRAW(i,:);
end
