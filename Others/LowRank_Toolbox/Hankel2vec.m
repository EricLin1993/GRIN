function [ vec ] = Hankel2vec( H )
% Transform Hankel matrix H into a vector vec  
% Authored by Enping Lin
% 2020.6.29
%%--------------------------------------------------------------
 
   [r,c] = size(H);
   a = min(r,c);
   b = max(r,c);
   s = 0;

   
   if c > 1     %% testing Hankel formation of H
       Ht = hankel(H(:,1),H(r,:));  
       if norm(H-Ht)/(r*c)> 1e-4
%            warning('Input Matirx might not be strict Hankel Matrix !');
       end    
   end
   
   for it1 = 1:a 
       s = 0;
       for it2 =1:it1
          s =s + H(it1+1-it2,it2);
       end
       vec(it1) = s; 
%              vec(it1) = s/it1;
   end  
 
  for it1 = a+1:b 
       s = 0;
       if r <= c
           for it2 =it1-r+1:it1
              s = s + H(it1+1-it2,it2);
           end
       else
           for it2 =1:c
              s =s + H(it1+1-it2,it2);
           end
       end
       vec(it1) = s; 
%        vec(it1) = s/a;
  end   
  for it1 = b+1:c+r-1 
       s = 0;
       for it2 =it1-r+1:c
          s = s + H(it1+1-it2,it2);
       end
       vec(it1) = s;
%        vec(it1) = s/(c+r-it1);
   end    
   vec = vec(:);
   
end

