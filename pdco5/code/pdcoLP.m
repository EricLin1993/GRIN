function [x,y,z,obj] = pdcoLP(A,b,bl,bu,c,scale,wait)

% [xpdco,ypdco,zpdco] = pdcoLP( A,b,bl,bu,c,scale,wait );
% Solves an LP problem using pdco.m.
% The LP problem is of the form
%   min c'x  st  Ax = b, bl < x < bu.
%
%    scale = 0  suppresses scaling (OK if A,b,c are well scaled)
%    scale = 1  requests scaling (default)
%    wait  = 0  prevents pdco from waiting (default)
%    wait  = 1  asks pdco to wait to allow parameters to be reset.
    
%-----------------------------------------------------------------------
% 03 Jul 2006: LPnetlib.m derived from LPnetlib3.m, which was
%              developed by Michael Saunders and Leo Tenenblat
%              (SOL, Stanford University)
%              for experimenting with the Zoom strategy.
%              LPnetlib.m solves a single LP, with scaling as an option.
%
% 04 Jul 2006: Solution unscaled at the end.
% 13 Nov 2013: (San Kim) The Problem structure now matches the
%              revised format used by Tim Davis in the LPnetlib group.
% 27 Apr 2018: pdcoLP.m derived from LPnetlib.m.
%-----------------------------------------------------------------------

  if nargin < 6
     scale = 1;   % Default is to request scaling
  end
  if nargin < 7
     wait  = 0;   % Default is not to wait in pdco
  end

  [m,n]   = size(A);

  if scale
    iprint  = 1;
    scltol  = 0.9;
    [cscale,rscale] = gmscale(A,iprint,scltol);

    C = spdiags(cscale,0,n,n);   Cinv = spdiags(1./cscale,0,n,n);
    R = spdiags(rscale,0,m,m);   Rinv = spdiags(1./rscale,0,m,m);
  
    A  = Rinv*A*Cinv;
    b  = b ./rscale;
    c  = c ./cscale;
    bl = bl.*cscale;
    bu = bu.*cscale;
  end

  fixed   = find(bl==bu);
  blpos   = find(bl< bu & bl>0);
  buneg   = find(bl< bu & bu<0);
  rhs     = b - A(:,fixed)*bl(fixed) ...
              - A(:,blpos)*bl(blpos) ...
              - A(:,buneg)*bu(buneg);

  bscale  = norm(rhs,inf);   bscale  = max(bscale,1);
  oscale  = norm(c,inf);     oscale  = max(oscale,1);

  if scale
    b       = b /bscale;
    bl      = bl/bscale;
    bu      = bu/bscale;
    c       = c /oscale;
    fprintf('\n\n  Final b and c scales:  %11.4e     %11.4e', bscale, oscale)
  end
  
  
  c10     = zeros(n,1);
  c20     = zeros(n,1);
  d0      = zeros(m,1);
  gamma   = 1e-3;           % Primal regularization
  delta   = 1e-3;           % 1e-3 or 1e-4 for LP;  1 for Least squares.
  d1      = gamma;          % Can be scalar if D1 = d1*I.
  d2      = delta*ones(m,1);

%-----------------------------------------------------------------------
% Solve LP accurately for comparison.
%-----------------------------------------------------------------------
  disp(' ');  disp(' ');  disp(' ');  disp(' ');  disp(' ');  disp(' ');
  disp('==============================================================')
  x0      = zeros(n,1);      % Initial x
  x0      = max(x0,bl);
  x0      = min(x0,bu);
  y0      = zeros(m,1);      % Initial y
  z0      = zeros(n,1);      % Initial z

  if scale
    xsize   = 1;               % Estimate of norm(x,inf) at solution
    zsize   = 1;               % Estimate of norm(z,inf) at solution
  else
    xsize   = bscale;
    zsize   = oscale;
  end

  options           = pdcoSet;
  options.mu0       = 1e-0;  % An absolute value
  options.Method    = 1;     % 1=chol  2=QR  3=LSQR
  options.LSQRatol1 = 1e-10; % For LPs, LSQR must solve quite accurately
  options.LSQRMaxIter = 100; % Default = 10 (*m)
  options.OptTol    = 1e-6;
  options.FeaTol    = 1e-6;
  options.wait      = wait;  % 1 allows options to be reviewed before solve
  options.MaxIter   = 100;

  [x,y,z,inform,PDitns,CGitns,time] = ...
    pdco( c,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize );

  if scale          % unscale x,y,z
    b  = b *bscale;
    bl = bl*bscale;
    bu = bu*bscale;
    c  = c *oscale;
    x  = x *bscale;
    y  = y *oscale;
    z  = z *oscale;

    A  = R*A*C;
    b  = b .*rscale;
    c  = c .*cscale;
    bl = bl./cscale;
    bu = bu./cscale;
    x  = x ./cscale;
    y  = y ./rscale;
    z  = z ./cscale;
  end

  obj = c'*x;
  fprintf('\n Unscaled linear objective = %12.7e\n', obj)

  if wait
    disp(' ');   disp('pdco solution:')

    fprintf('\n    j       bl           x            bu           z\n')
    for j=1:n
        fprintf('%5i %12.4e %12.4e %12.4e %12.4e\n', j, bl(j), x(j), bu(j), z(j))
    end

    fprintf('\n    abs(y) > 1e-4\n')
    Y = find(abs(y) > 1e-4);
    for k=1:length(Y)
        j = Y(k);
        fprintf('%5i %12.4e\n', j, y(j))
    end

    fprintf('\n Linear objective =%12.4e\n', c'*x)

    disp('Waiting in pdcoLP in case you want to look further at the solution')
    keyboard
end
