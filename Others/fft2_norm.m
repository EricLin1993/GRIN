function [ y ] = fft2_norm( x )
%
     [r,c]=size(x);
     y = fftshift(fft2(x)/sqrt(r*c));

end

