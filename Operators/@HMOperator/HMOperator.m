function [ res ] = HMOperator( r,c )

   res.adjoint = 0;
   res.r = r;
   res.c = c;
   % Register this variable as a HankelMatrix class
   res = class(res,'HMOperator');

end

