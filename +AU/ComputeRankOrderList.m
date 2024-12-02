function [testVals,trainVals,order,RAW,AUC]=ComputeRankOrderList(Data,Condition,NumTests,varargin)

% Test to see if one condition (column) is reliably the larger then all
% other columns across random 50/5-0 NumTests in the data

% Data : Number of samples X number of conditions.
% labels Optional condition names
% NumTests  Number of NumTests
[varargin,ComputeAUC]   = Utilities.ProcVarargin(varargin,'ComputeAUC');
%%
DefaultNumTests=100;
if nargin<3; NumTests=[]; end
%%
CondIDs=unique(Condition);
N=(histc(Condition,CondIDs)');
NMax=max(N);
NConds=length(CondIDs);
% If the amount of data is reasonable, explicitly generate all possible
% combinations, and draw from the list (otherwise just sample randomly)
for i=1:NConds
    if N(i)<=12
        
        useExplicit(i)=1;
        C{i} = nchoosek(1:N(i),ceil(N(i)/2));
        C{i}=C{i}(randperm(size(C{i},1)),:);
        
        if ~isempty(NumTests) && size(C{i},1)<NumTests
            C{i}=repmat(C{i},ceil(NumTests/size(C{i},1)),1);
        end
        % can't do more NumTests then unique combinations.
    else
        useExplicit(i)=0;
    end
    
end
%%
if isempty(NumTests)
    tmp=cell2mat(cellfun(@(x)size(x,1),C,'UniformOutput',false))
    tmp(tmp==0)=[];
    if isempty(tmp)
        NumTests=DefaultNumTests;
    else
        NumTests=min(tmp);
    end
end
%%
AUC=nan(NConds,NConds,NumTests);
for i=1:NumTests
    for j=1:NConds
        if useExplicit(j)
            TrainInds=C{j}(i,:); % grab 50% of the data
        else
            TrainInds=randperm(N(j),ceil(N(j)/2)); % grab 50% of the data
        end
        TestInds=setdiff(1:N(j),TrainInds);  % grab remaining data data
        
        condIDX=find(Condition==CondIDs(j));
        
        testRawData{j}=Data(condIDX(TestInds),:);
        train(j)=mean(Data(condIDX(TrainInds),:),1);
        test(j)=mean(Data(condIDX(TestInds),:),1);
    end
    
    [trainValsRAW(i,:),order(i,:)]=sort(train,'descend');
    testValsRAW(i,:)=test(order(i,:));
    
   if ComputeAUC
       
    testRawDataOrdered=testRawData(order(i,:));
    for ii=1:(NConds-1)
        for jj=ii+1:NConds
            Labels=[testRawDataOrdered{ii}*0+1;testRawDataOrdered{jj}*0];
            [~,~,~,AUC(ii,jj,i)] = perfcurve(Labels,[testRawDataOrdered{ii};testRawDataOrdered{jj}],1);
%             [~,~,~,AUCshuf(ii,jj,i)] = perfcurve(Shuffle(Labels),[testRawDataOrdered{ii};testRawDataOrdered{jj}],1);
        end
    end
    
   end
    
    
    
    RAW(i).train=train;
    RAW(i).test=test;
    RAW(i).trainValsRAW=trainValsRAW(i,:);
    RAW(i).testValsRAW=testValsRAW(i,:);
end
%%
testVals=mean(testValsRAW,1);
trainVals=mean(trainValsRAW,1);

% N=nan(NConds,NConds);
% for ii=1:(NConds-1)
%         for jj=ii+1:NConds
%             mean(squeeze(AUC(ii,jj,:)))>prctile(AUCshuf(ii,jj,:),[95]);
%             N(ii,jj)=nnz(mean(squeeze(AUC(ii,jj,:)))>squeeze(AUCshuf(ii,jj,:)))/NumTests;
% %             [N(ii,jj)]=signrank(squeeze(AUC(ii,jj,:)),.5)
%             
%         end
%     end