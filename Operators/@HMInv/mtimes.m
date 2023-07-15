function res = mtimes(A,X)

    if A.adjoint == 0 % M*X
        
      
        h = Hankel2vec( ones(A.Hr,A.Hc) );
        h = h(:); 
        mask = A.mask;
        mask = mask(:);
        temp = 1./(A.mu*mask+h);
        res = temp.*X;
    else % Mh*X
        
    end   
end
