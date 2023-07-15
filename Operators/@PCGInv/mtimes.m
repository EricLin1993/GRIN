function res = mtimes(A,X)

    if A.adjoint == 0 % M*X
         OPs = A.OPs;  
         mu = A.mu;
         x0 = A.x0;
         pcgtol = 0.001;
         pcgmaxi = 5000;
         temp = 2*ones(size(x0))+1./(mu*x0.^2);
         p1 = 1./(temp);
     [res,pflg,prelres,pitr,presvec] = ...
            pcg(@AXfun,X,pcgtol,pcgmaxi,@Mfun,...
                [],x0,OPs,mu,p1);%?
    else % Mh*X
        
    end
    
return    
%------------------------------------------------------------
%       COMPUTE AX (PCG)
%------------------------------------------------------------
function [y] = AXfun(x,OPs,mu,p1)
    y = 0;
    for it = 1:length(OPs)
       Atemp = OPs{it};
       y = y+mu*(Atemp'*(Atemp*x));
    end
   


%------------------------------------------------------------
%       COMPUTE P^{-1}X (PCG)
%------------------------------------------------------------
function [y] = Mfun(x,OPs,mu,p1)
    
    y = p1.*x;

