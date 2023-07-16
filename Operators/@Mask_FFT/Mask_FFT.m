function res =Mask_FFT(mask)

res.adjoint = 0;
res.mask=mask;
% res.row=size(mask,1);
% res.column=size(mask,2);
res = class(res,'Mask_FFT');
