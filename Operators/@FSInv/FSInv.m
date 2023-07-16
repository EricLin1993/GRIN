function [ res ] = FSInv( mask,mu )
% FIInv serves as the inversion of  (F'M'MF+mu)
   res.adjoint = 0;
   res.mask = mask;
   res.mu = mu;
   % Register this variable as a FSInv class
   res = class(res,'FSInv');

end

