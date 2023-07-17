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


% s0=dephase_fun(s,phc0,phc1);
phc0=phc0*pi/180;              %convert degree to radian
phc1=phc1*pi/180;              %convert degree to radian
[m,n]=size(s);
a_num=[1:m]./m;
a=phc0.*ones(1,m)+a_num.*phc1; % calculate a(i,j)
re=real(s);                % get the real part of complex numbers
im=imag(s);                % get the imaginary part of complex numbers
re_new=re.*cos(a)'-im.*sin(a)';  % get the new real part by phase correction
im_new=re.*sin(a)'+im.*cos(a)';  % get the new imaginary part by phase correction

s0=re_new+im_new*1i;         % dephased nmr spectral data


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

