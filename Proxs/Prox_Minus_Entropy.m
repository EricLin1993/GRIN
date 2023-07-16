function [ Y ] = Prox_Minus_Entropy( beta,X )
%
   fun = @entropy_eq;
   tol = 1e-5;
   nmax = 1000;
   Y = zeros(size(X));
   for it = 1:length(X(:))
        % the zero point lies between exp(-1) and X(it) and satisfies >0;
        t(1) = exp(-1);
        t(2) = max(X(it),eps);
        if t(1)<= t(2)
           a = t(1);
           b = t(2);
        else
           a = t(2);
           b = t(1);
        end    
        [x0,res,niter] = bisection(fun,a,b,tol,nmax,beta,X(it));
        Y(it) = x0;
   end
end

