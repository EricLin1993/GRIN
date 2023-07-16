function [ Y ] = Prox_L1( beta,X )
%
    Y = sign(X).*max(abs(X)-beta,0);
%     Y(Y<0)=0;

end

