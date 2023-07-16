function res = ctranspose(a)
a.adjoint = xor(a.adjoint,1);% if a.adjoint==1, xor(a.adjoint,1)=0; if a.adjoint==0, xor(a.adjoint,1)=1;
res = a;

