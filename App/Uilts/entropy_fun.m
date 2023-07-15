function f = entropy_fun(x,s)
%function f = entropy_fun(x,s)
%written by Chen Li 
%input:   x - PHC0 and PHC1
%         s - NMR Spectral data
%output   f - entropy value (Using the first derivative)

%initial parameters
stepsize=1; 
func_type=1;

%dephase
[N,L]=size(s);
phc0=x(1);phc1=x(2);
% Rc=i*(pi/180);
% phcorr=(exp(Rc*(phc0+phc1*[1:length(s)]/length(s))))';
% s0=s.*phcorr;
s0=dephase_fun(s,phc0,phc1);
s=real(s0);

% Calculation of first derivatives 
if (func_type == 1)
    ds1 = abs((s(3:L)-s(1:L-2))/(stepsize*2));
else
    ds1 = ((s(3:L)-s(1:L-2))/(stepsize*2)).^2;
end  
p1 = ds1./sum(ds1);

%Calculation of Entropy
[M,K]=size(p1);
for j0=1:M
    for j=1:K
        if (p1(j0,j)==0)%in case of ln(0)
           p1(j0,j)=1; 
        end
    end
end
h1  = -p1.*log(p1);
H1  = sum(h1);
%Calculation of penalty
Pfun	= 0.0;
as      = s - abs(s);
sumas   = sum(sum(as));
if (sumas < 0)
   Pfun = Pfun + sum(sum((as./2).^2));
end
P       = 1000*Pfun;

% The value of objective function

f = H1+P;

