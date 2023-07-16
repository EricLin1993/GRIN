function res = mtimes(a,b)

if isa(a,'Rice_Wavelet') == 0
    error('In  A.*B only A can be Rice_Wavelet');
end

if a.adjoint
    b = reshape(b,a.row,a.col*(3*a.wavescale+1));
    col=size(b,2)./(3*a.wavescale+1);
    yl=b(:,1:col);
    yh=b(:,col+1:end);
    if a.complex_yes==1
        res=mirdwt(real(yl),real(yh),a.h,a.wavescale)...
        + 1j * mirdwt(imag(yl),imag(yh),a.h,a.wavescale);
    else
        res=mirdwt(yl,yh,a.h,a.wavescale);
    end
else
    b = reshape(b,a.row,a.col);
    if a.complex_yes
        [yl_real,yh_real]=mrdwt(real(b),a.h,a.wavescale);
        [yl_imag,yh_imag]=mrdwt(imag(b),a.h,a.wavescale);
        res=[yl_real+1j*yl_imag yh_real+1j*yh_imag];
    else
        [yl,yh]=mrdwt(b,a.h,a.wavescale);
        res=[yl yh];
    end
end