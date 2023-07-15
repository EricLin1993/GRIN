function [y ] = Minus_Entropy( x )
%  -entropy
     x =  x(:);
     x(x<=0) =eps;
     logx = log(x);
     y = sum( x.*logx );
end

