function [ Y ] = Prox_LogBarrier( beta,X )
%
    Y = X/2.*sqrt(X.^2+4/beta);
end

