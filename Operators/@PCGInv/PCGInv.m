function [ res ] = PCGInv( OPs,x0,mu )

   res.adjoint = 0;
   res.OPs = OPs;
   res.mu = mu;
   res.x0 = x0; 
   % Register this variable as a PCGInv class
   res = class(res,'PCGInv');

end

