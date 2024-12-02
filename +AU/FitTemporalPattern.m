function [FitCC,CC,coef]=FitTemporalPattern(Data,Labels,NumTests)

% How well is the neural response to an "action" that combines two
% movements explained by the neural response to those movements alone?
%%
X1=Data{1};
X2=Data{2};
Y=Data{3};

if nargin<2 || isempty(Labels); Labels=1:size(X1,2); end
if nargin<3; NumTests=500; end

N=size(X1,1);

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
    
    X1train=mean(X1(TrainInds,:),1);
    X1test=mean(X1(TestInds,:),1);
    
    X2train=mean(X2(TrainInds,:),1);
    X2test=mean(X2(TestInds,:),1);
    
    Ytrain=mean(Y(TrainInds,:),1);
    Ytest=mean(Y(TestInds,:),1);
    
    coef(i,:)=Ytrain'\[X1train' X2train'];
        FitCC(i)=corr(Ytest',(coef(i,:)*[X1test' X2test']')');
    
    CC(i)=corr(Ytrain',Ytest');
    
%     RAW(i).train=train;
%     RAW(i).test=test;
%     RAW(i).trainValsRAW=trainValsRAW(i,:);
%     RAW(i).testValsRAW=testValsRAW(i,:);
end
