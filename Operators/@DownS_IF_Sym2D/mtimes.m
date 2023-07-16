function res = mtimes(A,x)

if isa(A,'DownS_IF_Sym2D') == 0
    error('In  A.*x only A can be DownS_IF_Sym2D');
end

if A.adjoint
    XT=zeros(size(A.mask));
    XT(logical(A.mask))=x;
    [r,c]=size(XT);
    XS = fft2(XT)./sqrt(r*c);
%     XS =XT;
    
    res = Sym2DToVecPer( XS );
    
else
    XS = Vec2Sym2DPer( x,A.column );
    [r,c]=size(XS);
    XT=ifft2(XS).*sqrt(r*c);
%     XT = XS;
    
    res=XT(logical(A.mask));
end