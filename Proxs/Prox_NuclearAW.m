function [ Y ] = Prox_Nuclear( beta,X )
%
   [u,s,v] = svd(X,'econ');
   sv = abs(diag(s));
   sv = max(sv-beta*(1./(sv+eps)),0);
%    n = sum(sv>0);
%    Y = u(:,1:n)*diag(sv(1:n))*(v(:,1:n))';
   Y = u*diag(sv)*v';
end

