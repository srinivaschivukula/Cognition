function [FRout,Lout,NTrials,NTimeBins,NUnits]=CollapseNeuralData(FR,L)

% Collapses neural data output by Analyze.getNeuralData from NTrials x
% NTimebins x Nunits -> NTimebins*NTrials x Nunits
% useful for dimensionality reduction etx


if~iscell(FR)
    FR={FR};
    if nargin>1
        L={L};
    end
    cb=1;
else
    cb=0;
end

for i=1:length(FR)
    
    NTrials(i)=size(FR{i},1);
    NTimeBins(i)=size(FR{i},2);
    NUnits(i)=size(FR{i},3);
    FRout{i}=reshape(permute(FR{i},[2 1 3]),NTrials(i)*NTimeBins(i),NUnits(i));
    if nargin>1
        tmpL=repmat(L{i},1,NTimeBins(i));
        Lout{i}=reshape(permute(tmpL,[2 1]),NTrials(i)*NTimeBins(i),1);
    end
    
end


if cb
    FRout=FRout{1};
    if nargin>1
        Lout=Lout{1};
    end
end
if nargin<2
    Lout=[];
end
    
    
%
% % & back
% FR3=permute(reshape(FR2,NTimeBins,NTrials,NUnits),[2 1 3]);
%
%
% [FR{i}(1,1:16,1),FR{i}(2,1:16,1)]
% FR2(1:32,1)'