function res =Rice_Wavelet(row,col,wavescale,filter_type,filterSize,complex_yes)

res.adjoint = 0;
res.row = row;
res.col = col;
res.wavescale = wavescale;
res.complex_yes = complex_yes;
res.h = MakeONFilter(filter_type,filterSize);
res = class(res,'Rice_Wavelet');
