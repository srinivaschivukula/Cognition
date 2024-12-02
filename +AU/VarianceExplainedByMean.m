function [R2,SNR]=VarianceExplainedByMean(X)

mu=mean(X,2);



y=repmat(mu,size(X,2),1);


R2=AU.ComputeR2(X(:),y);

v=X(:)-y;
dom=max(mu)-min(mu);

SNR=dom/std(v);