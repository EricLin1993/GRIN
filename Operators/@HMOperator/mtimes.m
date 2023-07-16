function res = mtimes(H,x)

if H.adjoint == 0 % H*x
   
   if H.r+H.c-1 ~= length(x)
        error(['The size of Hankel Matrix Operator and the manuplated vector do not agree.']);
   end
   n = H.r;   
   res = hankel(x(1:n),x(n:end));
   
else % H'*x
   
   [r,c] = size(x);
   if r ~= H.r || c ~= H.c
      error(['The size of the adjoint Hankel Matrix Operator and the manuplated matrix do not agree.']); 
   end    
   a = min(r,c);
   b = max(r,c);
   s = 0;

   for it1 = 1:a 
       s = 0;
       for it2 =1:it1
          s =s + x(it1+1-it2,it2);
       end
       vec(it1) = s; 
%              vec(it1) = s/it1;
   end  
 
  for it1 = a+1:b 
       s = 0;
       if r <= c
           for it2 =it1-r+1:it1
              s = s + x(it1+1-it2,it2);
           end
       else
           for it2 =1:c
              s =s + x(it1+1-it2,it2);
           end
       end
       vec(it1) = s; 
%        vec(it1) = s/a;
  end   
  for it1 = b+1:c+r-1 
       s = 0;
       for it2 =it1-r+1:c
          s = s + x(it1+1-it2,it2);
       end
       vec(it1) = s;
%        vec(it1) = s/(c+r-it1);
   end    
   res = vec(:);
    
end
