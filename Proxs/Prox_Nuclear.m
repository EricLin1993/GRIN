function [ Y ] = Prox_Nuclear( beta,X )
%
   [u,s,v] = svd(X,'econ');
   sv = diag(s);
   sv = max(sv-beta,0);
   n = sum(sv>0);
   Y = u(:,1:n)*diag(sv(1:n))*(v(:,1:n))';
end

