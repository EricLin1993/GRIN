function res = mtimes(a,b)

if isa(a,'Fu_downsample') == 0
    error('In  A.*B only A can be Fu_downsample');
end

if a.adjoint
    kspace=zeros(size(a.mask));
    kspace(logical(a.mask))=b;
    res=ifft2_norm(kspace);
    
else
    b=reshape(b,a.row,a.column);
    kspace=fft2_norm(b);
    res=kspace(logical(a.mask));
end