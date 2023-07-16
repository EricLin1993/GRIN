function [ res ] = UFS_Inv( W )

   res.adjoint = 0;
   res.W = W;
   % Register this variable as a UFS_Inv class
   res = class(res,'UFS_Inv');

end

