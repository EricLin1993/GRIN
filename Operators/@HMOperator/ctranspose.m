function  res = ctranspose(H)

 H.adjoint = xor(H.adjoint,1);
 res = H;