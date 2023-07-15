function [ res ] = ImInv( mask,mu )

   res.adjoint = 0;
   res.mask = mask;
   res.mu = mu;
   
   % Register this variable as a ImInv class
   res = class(res,'ImInv');

end

