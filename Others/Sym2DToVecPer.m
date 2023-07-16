function [ res ] = Sym2DToVecPer( xc )

    uflag=triu(ones(size(xc)),0);
    flag=uflag-diag(0.5*ones(size(xc,1),1));
    X = (xc+xc.').*flag;
    res = X(logical(uflag));

end

