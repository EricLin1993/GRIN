function [ y ] = ifft2_norm( x )
%
     [r,c]=size(x);
     y = ifft2(ifftshift(x))*sqrt(r*c);

end

