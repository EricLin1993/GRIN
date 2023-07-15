function res =Fu_downsample(mask,row,column)

res.adjoint = 0;
res.mask=mask;
res.row=row;
res.column=column;
res = class(res,'Fu_downsample');
