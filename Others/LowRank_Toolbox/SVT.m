function [ Y ] = SVT( X,tau )
% tau 

   [U,S,V] = svd(X,'econ');
   sv = diag(S);
   thr = tau  ;   
   Y = U*diag(max(sv-thr,0))*V';

end

