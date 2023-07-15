function [ H ] = vec2Hankel( x,n )
% Transform vector x into Hankel matrix H with n rows  
% Authored by Enping Lin
% 2020.6.29
%%--------------------------------------------------------------
   if exist('n','var') == 0
        n = round(length(x)/2);
   end
   n = round(n);
   if (n < 1) or (n > length(x))
        error('n must be no less than 1 and no greater than the length of x');
   end    
   H = hankel(x(1:n),x(n:end));
   
end

