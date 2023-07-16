function [y] = L2_1(x)
% L2_1 norm 
%   
   for it = 1:size(x,2)
         y = y+norm(x(:,it));
   end    
end

