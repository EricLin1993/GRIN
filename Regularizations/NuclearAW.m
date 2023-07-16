function [ y ] = NuclearAW( x )

    [~,s,~]=svd(x,'econ');
    y = sum(diag(s));

end

