function [R2,p,R2shuff,ErrRedCoef,pErrRed,MSE]=ComputeR2(X,Y,nShuffles)
[R2,SSresid]=DoIt(X,Y);

n=length(X);
MSE=sum((X-Y).^2)./n;

if nargin>2
    
    for i=1:nShuffles
        [R2shuff(i),SSresidShuff(i)]=DoIt(X,Shuffle(Y));
    end
   p=1-(nnz(R2> R2shuff)/nShuffles);
   
   ErrRedCoef=1-SSresid/mean(SSresidShuff);
   pErrRed=1-(nnz(SSresid< SSresidShuff)/nShuffles);
end


function [R2,SSresid]=DoIt(X,Y)
    SStot=sum( (X-mean(X)).^2 );
%     SSreg=sum( (Y-mean(X)).^2 );
    SSresid=sum((X-Y).^2);
    
    
    R2=1-SSresid./SStot;