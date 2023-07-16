function [ res ] = Vec2S2DMat_Inv( XC )

    uflag=triu(ones(size(XC)),0);
    X = (XC+XC.')./2;
    res = X(logical(uflag));
end

