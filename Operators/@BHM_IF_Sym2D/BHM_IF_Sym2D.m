function [ res ] = BHM_IF_Sym2D( n )

   res.adjoint = 0;
   res.n = n;
   % Register this variable as a HM_IF_Sym2D class
   res = class(res,'BHM_IF_Sym2D');

end

