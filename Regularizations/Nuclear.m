function [ y ] = Nuclear( x )

    [~,s,~]=svd(x,'econ');
    y = sum(diag(s));

end

