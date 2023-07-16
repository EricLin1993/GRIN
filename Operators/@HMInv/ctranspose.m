function  res = ctranspose(mask)

mask.adjoint = xor(mask.adjoint,1);
res = mask;