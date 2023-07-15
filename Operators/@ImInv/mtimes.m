function res = mtimes(A,X)

    if A.adjoint == 0 % M*X
        
        res = fft2_norm(X);
        temp = 1./(A.mask+A.mu);
        res = ifft2_norm(temp.*res);
    else % Mh*X
        
    end   
end
