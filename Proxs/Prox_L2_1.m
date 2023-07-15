function [ Y ] = Prox_L2_1( beta,X )

%     Y = sign(X).*max(abs(X)-beta,0);
%     Y(Y<0)=0;
%   20220505 Enping Lin
    Y = X;
    for it =1:size(X,2)
        ft = norm(X(:,it));
        Y(:,it) = X(:,it)/ft*max(ft-beta,0);
    end    
end

