function res = mtimes(A,x)

    if A.adjoint == 0 % M*X
        c = size(A.W,2);
        Temp = S2DMat2Vec_Inv( x,c );    %Sh inv
        [r,c] = size(Temp);
        Temp = ifft2(Temp).*sqrt(r*c);
        Temp = Temp./(A.W);
        Temp = fft2(Temp)./sqrt(r*c);
        res = Vec2S2DMat_Inv( Temp);
      
     
    else % Mh*X
        
    end   
end
