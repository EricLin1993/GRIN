function [ Nu ] = NuclearNorm( X )
%  return Nu ,the nuclear norm of X.

   [~,sv,~] = svd(X,'econ');
    Nu = sum(diag(sv));

end

