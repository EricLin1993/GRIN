function [y] = LogBarrier(x)
% logrithm barrier 
%  Enping Lin 20220505
     y = sum(-log(x),'all');
end

