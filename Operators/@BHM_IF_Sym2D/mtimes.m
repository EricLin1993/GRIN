function res = mtimes(H,x)

if H.adjoint == 0 % H*x
   
    S = Vec2Sym2DPer( x,H.n );
    [r,c]=size(S);
    X = ifft2(S).*sqrt(r*c);   
%     X = S;

    r = round(H.n/2);
    c = r;
    HBX =  Matrix2BHankel(X,r,c);
    res = HBX;
   
else % H'*x
    r = round(H.n/2);
    c = H.n+1-r;
    [ X ] = BHankel2Matrix( x,r,c);
    [r,c]=size(X);
    S = fft2(X)./sqrt(r*c);
%     S = X;
    x = Sym2DToVecPer( S );
    res = x;
    
    
end
