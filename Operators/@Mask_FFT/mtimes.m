function res = mtimes(a,b)

if isa(a,'Mask_FFT') == 0
    error('In  A*x only A can be Mask_FFT');
end

if a.adjoint
    sz = size(a.mask);
    kspace=zeros(sz);
    kspace(logical(a.mask))=b;
    res=fftn(kspace)./sqrt(prod(sz));
    res = res(:);
else
    sz = size(a.mask);
    b=reshape(b,sz);
    kspace=ifftn(b)*sqrt(prod(sz));
    res=kspace(logical(a.mask));
end