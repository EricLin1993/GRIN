function [ X ] = BHankel2Matrix( BH,r,c)
% Transform the Block Hankel matrix into a  matrix    
%  Input: BH : the Block Hankel matrix
%         r : the number of the sub Hankel matrix in the row of the Block Hankel matrix
%         c : the number of the sub Hankel matrix in the column of the Block Hankel matrix
%  Output: X : the constructed matrix

% Authored by Enping Lin
% 2021.3.29
%%--------------------------------------------------------------

   a = min(r,c);   
   b = max(r,c);
   [shs] = size(BH)./[r,c];
   sr = shs(1);sc = shs(2);
   
%    if c > 1     %% testing Hankel formation of H
%        Ht = hankel(H(:,1),H(r,:));  
%        if norm(H-Ht)/(r*c)> 1e-4
% %            warning('Input Matirx might not be strict Hankel Matrix !');
%        end    
%    end
   
   for it1 = 1:a 
       s = zeros(sr,sc);
       for it2 =1:it1
          br = it1+1-it2; 
          bc = it2;
          s =s + BH((br-1)*sr+1:br*sr,(bc-1)*sc+1:bc*sc);
       end
         Hc(:,:,it1) = s;
%          Hc(:,:,it1) = s/it1;
   end  
 
  for it1 = a+1:b 
       s = zeros(sr,sc);
       if r <= c
           for it2 =it1-r+1:it1
              br =  it1+1-it2;
              bc = it2;
              s = s + BH((br-1)*sr+1:br*sr,(bc-1)*sc+1:bc*sc);
           end
       else
           for it2 =1:c
              br =  it1+1-it2;
              bc = it2;
              s = s + BH((br-1)*sr+1:br*sr,(bc-1)*sc+1:bc*sc);
           end
       end
       Hc(:,:,it1) = s;
%        Hc(:,:,it1) = s/a;
      
  end   
  for it1 = b+1:c+r-1 
       s = zeros(sr,sc);
       for it2 =it1-r+1:c
          br =  it1+1-it2;
          bc = it2;
          s = s + BH((br-1)*sr+1:br*sr,(bc-1)*sc+1:bc*sc);
       end
       Hc(:,:,it1) = s;
%        Hc(:,:,it1) = s/(c+r-it1);

  end    
% ------------ Sub-Hankel to the collumn of Matirx --------------
   for it = 1:size(Hc,3)
       X(:,it) = Hankel2vec( squeeze(Hc(:,:,it)) );
   end

end

