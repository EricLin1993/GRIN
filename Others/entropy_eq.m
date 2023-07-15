function [y] = entropy_eq(x,beta,b)
%  this function is  the derivative of fx = 0.5*(x-b).^2+lambda*x*log(x)
%  which is the individual proximal function of -entropy(x) at the point b
    y = x-b+beta*(log(abs(x))+1);

end