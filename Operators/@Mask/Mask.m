function [ res ] = Mask( M )

   res.adjoint = 0;
   res.size = size(M);
   res.indice = find(M~=0);
   res.M = M;
   % Register this variable as a Mask class
   res = class(res,'Mask');

end

