function [ res ] = S2DMat2Vec_Inv( xc,c )


    uflag=triu(ones(c),0);
    X_Symetric = zeros(c);
    X_Symetric(logical(uflag)) = xc;
    Temp = (X_Symetric+X_Symetric.')/4;
    res = Temp+diag(diag(Temp));




end

