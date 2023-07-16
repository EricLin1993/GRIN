function [ BH ] = Matrix2BHankel(X,n1,n2)
% Transform Matrix X into a Block Hankel matrix BH  
% Input :  
%         X : the matirx to be permuted in Block Hankel matrix
%         n1 : the row number of the sub Hankel matrix
%         n2 : the number of the sub Hankel matrix in the row of the Block Hankel matrix
% Authored by Enping Lin
% 2021.3.29
%%--------------------------------------------------------------
   [r,c] = size(X);
   cl = 1:c;
%    if exist('n','var') == 0
%        n1 = ceil(r/2);
%        n2 = ceil(c/2);
%    else
%        n1 = n(1);n2 = n(2);
%    end
   hcl = r-n1+1; shcl = c-n2+1;
   BH = zeros(n1*n2,hcl*shcl);
   for it = 1:c
       Hc(:,:,it) = hankel(X(1:n1,it),X(n1:end,it));
   end    
   Hcl = hankel(cl(1:n2),cl(n2:end));
   for it1 = 1:size(Hcl,1)
      for it2 = 1:size(Hcl,2) 
         BH( (it1-1)*n1+1:it1*n1,(it2-1)*hcl+1:it2*hcl ) = squeeze(Hc(:,:,Hcl(it1,it2)));
      end
   end  
end

