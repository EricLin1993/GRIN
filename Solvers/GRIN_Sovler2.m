function [ GR_Out ] = GRIN_Sovler2( GR_InPara)
%--------------------------------------------------------------------------
%  This program aims at  solving generally regularized inverse problem
%  formaluted as:
%      1/2*||K*x-y||+ lambda_i*Ri(Ai*x) for i=1,2,3,...,n 

%  Enping Lin  
%  enpinglin2022@126.com
%  last revision :2022.5.16
%
%--------------------------------------------------------------------------
 y = GR_InPara.y;
 K = GR_InPara.K;
 Kh = K'; 
 RegTerm = GR_InPara.RegTerm;
 RegTermNum = length(RegTerm);
 Khy = Kh*y;
 lambda(1:RegTermNum) = 1;
 if ~isfield(GR_InPara,'lambda')
    
 else
    lambda(1:length(GR_InPara.lambda)) = GR_InPara.lambda;
end
if ~isfield(GR_InPara,'mu')
    mu = 1;
 else
    mu = GR_InPara.mu;
end
% A setting
A = cell(RegTermNum,1);  
 
for it =1:RegTermNum  
       A{it} = 1; 
end 
if ~isfield(GR_InPara,'A')
 
else  
    for it =1:length(GR_InPara.A)  
       A{it} = GR_InPara.A{it};  
    end 
end 
for it = 1:RegTermNum   
      eval(['A',num2str(it),'=A{it};']);
      eval(['Ah',num2str(it),'=A',num2str(it),''';']);       
end 


% setting of inverse operator for x 

    if ~isfield(GR_InPara,'Lxp')
         PartLxp = Kh*K; 
         for it = 1:RegTermNum    
              eval(['PartLxp = PartLxp +', 'mu*Ah',num2str(it),'*A',num2str(it),'*eye(size(PartLxp));']);  % revised 20220314
         end 

         Lxp = inv(PartLxp);                         
    else
         Lxp = GR_InPara.Lxp;
    end

% the condition to end computation
if ~isfield(GR_InPara,'maxit')
    maxit = 1000;
 else
    maxit = GR_InPara.maxit;
end
if ~isfield(GR_InPara,'tol')  
    tol  = 1e-5;
 else
    tol  = GR_InPara.tol ;
end
if ~isfield(GR_InPara,'nnt')  
    nnt = 10;
 else
    nnt = GR_InPara.nnt ;
end
nn = 0; % counter of consecutive times where objval < tol 

%% -------------Initialization -----------------------  
%----------------------------------------------------------

if ~isfield(GR_InPara,'x0')
     x = ones(size(Kh*y));
     %  x = Kh*y; 
 else
     x = GR_InPara.x0;
end
  Z0 = K*x; 
  D0 = zeros(size(Z0)); 
  for it = 1:RegTermNum  
      eval(['Z',num2str(it),'= A',num2str(it),'*x;']);
      eval(['D',num2str(it),'= zeros(size(Z',num2str(it),'));']);             
  end     

  % calculating the objective value of initial x
    temp = 0.5*(norm(K*x-y)).^2;
    for it = 1:RegTermNum  
     eval(['temp=temp+lambda(',num2str(it),')*',RegTerm{it},'(A',num2str(it),'*x);']); 
    end   
    objval = zeros(maxit+1,1);
    objval(1) = temp; 
%======================Iteration=============================
%============================================================
GR_Out.restolYes = false;   % revised 20220516
tic  % start timer
    f = waitbar(0, 'Sovling:  0%', 'Name', 'Please Wait','CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); % revised 20220328 by chen
    setappdata(f,'canceling',0); % revised 20220405 by chen
    for itloop = 1:maxit
%         itloop    
        if getappdata(f,'canceling')
            break
        end
        str=['Sovling:  ',num2str(round(itloop/maxit*100)),'%']; % revised 20220404 by chen
        waitbar(itloop / maxit, f, str);  % revised 20220328 by chen
        
        if mod(itloop, 10) == 1
            Z0L = Z0;
            for it = 1:RegTermNum 
               eval(['Z',num2str(it),'L=Z',num2str(it),';']); 
            end    
        end
       % min x      
         Rxp =  Kh*(Z0-D0); %*mu % 
            for it = 1:RegTermNum  
                 eval(['Rxp = Rxp+(','Ah',num2str(it),'*(Z',num2str(it), '-D',num2str(it),'));']); % mu*  % 
            end

         xlast = x;

         x = Lxp*Rxp;

        % min Z
        Z0 = 1/(1+mu)*(y+mu*(K*x+D0));   % 
        for it = 1:RegTermNum     %  
            eval(['PX=','A',num2str(it),'*x+D',num2str(it),';']); %/mu  
            eval(['tau=lambda(it)/mu;']);  % 
            eval(['Z',num2str(it),'=Prox_',RegTerm{it},'(tau,PX);']);    % 
        end
        

        % updata D
        D0 = D0 + (K*x-Z0);  %  mu*  % 
        for it = 1:RegTermNum   % 
            eval(['D',num2str(it),'=D',num2str(it),'+(A',num2str(it),'*x-Z',num2str(it),');']); %mu*   %  
        end

%%     adjust the D and mu for a better convergence            
     if mod(itloop, 10) == 1
%         primal residual
       temp = K*x-Z0;  
       res_p2 = (norm(temp(:)))^2;
        for it = 1:RegTermNum 
            eval(['temp = Z',num2str(it),'-A',num2str(it),'*x;']); % 
            res_p2 = res_p2 + (norm(temp(:)))^2;
        end 
        res_p = sqrt(res_p2);

%         dual residual
       temp = Z0-Z0L; 
       res_d2 = (norm(temp(:)))^2;
        for it = 1:RegTermNum  
            eval(['temp = Z',num2str(it),'-Z',num2str(it),'L;']); 
            res_d2 = res_d2 + (norm(temp(:)))^2; 
        end 
        res_d = sqrt(res_d2);
        
        if res_p > 10*res_d %%
            mu = mu*2;
            D0 = D0/2; 
            for it = 1:RegTermNum 
              eval(['D',num2str(it),'=D',num2str(it),'/2;']);  
            end 

        elseif res_d > 10*res_p
            mu = mu/2;
            D0 = D0*2; 
            for it = 1:RegTermNum 
              eval(['D',num2str(it),'=D',num2str(it),'*2;']);  
            end  

       end
     end   
     % calculating the relative difference of two neighboring x
     xdif(itloop) = norm(x(:)-xlast(:))/norm(x(:));
     %  objective value curve
     temp = 0.5*(norm(K*x-y)).^2;
     for it = 1:RegTermNum  
         eval(['temp=temp+lambda(',num2str(it),')*',RegTerm{it},'(A',num2str(it),'*x);']); %
     end
     objval(itloop+1) = temp;
     rel_objval(itloop) = abs( (objval(itloop+1)-objval(itloop))/objval(itloop) );
     if rel_objval(itloop)<tol   
          nn = nn+1;
          if nn>=nnt 
             GR_Out.restolYes = true;
             break;
          end    
     else
           nn = 0;
     end    
%        if max(res_p,res_d)<res_tol   
%            GR_Out.restolYes = true;
%            break;
%        end    
%          
    end
    delete(f);  % revised 20220328 by chen
    GR_Out.ProcessingTime = toc/60;% end timer
    GR_Out.ProcessingTimeUnit = 'minute';
    if GR_Out.restolYes
      fprintf('The relative difference of objective value is less than %e for consecutive %d times and the processing is over.\n',tol,nnt)
      fprintf('The processing time is %f mins.\n',GR_Out.ProcessingTime )
    else
      fprintf('The maximum iteration number, %d, is reached and the processing is over.\n',maxit)
      fprintf('The processing time is %f mins.\n',GR_Out.ProcessingTime )
    end  
%     GR_Out.maxres = max(res_p,res_d); 
    GR_Out.x = x;
    GR_Out.xdif = xdif;
    GR_Out.objval = objval(1:itloop+1); % revised 20220516
    GR_Out.rel_objval = rel_objval(1:itloop); % revised 20220516
%     figure,plot(rel_objval),title('rel objval')
end

