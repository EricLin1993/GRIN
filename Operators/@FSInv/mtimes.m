function res = mtimes(A,X)

    if A.adjoint == 0 % M*X   
        sz = size(A.mask);
        X = reshape(X,sz);
        res = ifftn(X)./sqrt(prod(sz));
        temp = 1./(A.mask+A.mu);
        res = fftn(temp.*res)*sqrt(prod(sz));
        res = res(:);
    else % Mh*X
        
    end   
end
