function [x0,res,niter]=bisection(fun,a,b,tol,nmax,varargin)

    x=[a,(a+b)*0.5,b];
    fx(1)=feval(fun,x(1),varargin{:});
    fx(2)=feval(fun,x(2),varargin{:});
    fx(3)=feval(fun,x(3),varargin{:});
    if fx(1)*fx(3)>0
%         disp('the signs of the at the both ends of the interval are the same and the smaller one is returned');
        if fx(1) <= fx(3)
           x0 = x(1);
        else
           x0 = x(3);
        end 
        x0=a;res=0;niter=0;
        return;
    elseif fx(1)==0
        x0=a;res=0;niter=0;
        return
    elseif fx(3)==0
        x0=b;res=0;niter=0;
        return
    end
    niter=0;
    I=(b-a)*0.5;
    while I>=tol && niter<=nmax
        niter=niter+1;
        fx(1)=feval(fun,x(1),varargin{:});
        fx(2)=feval(fun,x(2),varargin{:});
        fx(3)=feval(fun,x(3),varargin{:});
        if sign(fx(1))*sign(fx(2))<0
            x(3)=x(2);x(2)= 0.5*(x(3)+x(1));
             fx(2)=feval(fun,x(2),varargin{:});
             fx(3)=feval(fun,x(3),varargin{:});
             I = I*0.5;
        elseif sign(fx(2))*sign(fx(3))<0
            x(1)=x(2);x(2)=0.5*(x(3)+x(1));
            fx(1)=feval(fun,x(1),varargin{:});
            fx(2)=feval(fun,x(2),varargin{:}); 
            I = I*0.5;
        else
            x0=x(2);
            I=0;
        end
    end
    x0=x(2);
    res = I;
    if niter>nmax
        fprintf('bisection stopped without converging to the desired tolerance','because the maxinum number of iterations was reached\n');
    end
end 