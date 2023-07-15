function [ Mask ] = Poisson( RandomSeed,Num_Indirect_Full,AccFactor,ExpFactor )
%
%  
rng('default');
rng(RandomSeed);
p=round(Num_Indirect_Full/AccFactor);% the length of sampled datas
v=zeros(Num_Indirect_Full,1);
ld=1.*Num_Indirect_Full/p;
w=2.0;
trial=0;
  while true
    k=0;n=0;trial=trial+1;
      while n<Num_Indirect_Full
            v(k+1)=n;
            k=k+1;
            n=n+1;
             %t=sin(n/(Num_Indirect_Full+1)*pi/2.0);
            t=exp(ExpFactor*n/(Num_Indirect_Full+1));
            t=(ld-1.)*w*t;        
            n=n+psn( t );
      end
    if (k==p) break;end
    if (k>p) w=w*1.02;else w=w/1.02;end
  end       
  v=v(1:p);
  v=v+1;
  Mask=zeros(Num_Indirect_Full,1);
  Mask(v)=1;
  

end

  function[k] = psn(lambda)
    u=rand();
    L=exp(-lambda);
    P=0;k=0;
    while P<u
       P=P+L*lambda^k/factorial(k);
       k=k+1;
    end
    k=k-1;
  end
  
  
  
  
