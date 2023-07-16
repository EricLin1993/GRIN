function res = mtimes(mask,X)

if mask.adjoint == 0 % M.*X
    res = X;
    res = res(mask.indice);
else % Mh.*X
    res = zeros(mask.size);
    res(mask.indice) = X(:);
end
