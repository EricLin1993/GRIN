function [ X_Symetric ] = Vec2Sym2DPer( xc,c )

    uflag=triu(ones(c),0);
    X_Symetric = zeros(c);
    X_Symetric(logical(uflag)) = xc;
    Temp = (X_Symetric+X_Symetric.');
    X_Symetric = Temp-0.5*diag(diag(Temp));

end

