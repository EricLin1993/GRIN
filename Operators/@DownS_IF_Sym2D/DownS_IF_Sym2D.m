function res =DownS_IF_Sym2D(mask)

res.adjoint = 0;
res.mask=mask;
res.row=size(mask,1);
res.column=size(mask,2);
res = class(res,'DownS_IF_Sym2D');
