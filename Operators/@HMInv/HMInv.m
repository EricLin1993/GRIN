function [ res ] = HMInv( mask,mu,Hr,Hc )

   res.adjoint = 0;
   res.mask = mask;
   res.mu = mu;
   res.Hr = Hr;
   res.Hc = Hc;
   % Register this variable as a HMInv class
   res = class(res,'HMInv');

end

