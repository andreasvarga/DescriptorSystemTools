function [sysx,mindist] = glinfldp(sys1,sys2,options)
%GLINFLDP  Solution of the least distance problem min||G1-X G2||_inf.
%
%   [SYSX,MINDIST] = GLINFLDP(SYS1,SYS2) returns the LTI descriptor 
%   realization SYSX of a stable solution X(lambda) of the  
%   2-block L-infinity least distance problem (LDP)
%
%      min ||G1(lambda)-X(lambda)  G2(lambda) ||_inf        (*)
%
%   given the generalized LTI realizations SYS1 and SYS2 of the 
%   transfer-function matrices G1(lambda) and G2(lambda), respectively. 
%   MINDIST is the achieved minimum distance corresponding to the optimal 
%   solution. If SYS2 = [], an 1-block LDP is solved. 
% 
%   [SYSX,MINDIST] = GLINFLDP(SYS,M2) specifies [G1(lambda) G2(lambda)] 
%   as the compound realization SYS = [SYS1 SYS2], where G2(lambda) has M2 
%   columns. If M2 = 0 (or M2 = []), an 1-block LDP is solved. 
%
%   [SYSX,MINDIST] = GLINFLDP(SYS1,SYS2,OPTIONS) uses OPTIONS structure 
%   to specify various user options, as follows:
%   OPTIONS.tol    - specifies the tolerance for rank determinations 
%                    (Default: [] - internally computed)
%   OPTIONS.gamma  - specifies gamma > norm(G2,inf), the desired 
%                    sub-optimality degree to solve the gamma-suboptimal 
%                    LDP 
%
%                    ||G1(lambda)-X(lambda)  G2(lambda) ||_inf <= gamma
%
%                    MINDIST is the achieved suboptimal distance. 
%                    (Default: [] - the optimal problem (*) is solved)
%   OPTIONS.reltol - specifies the relative tolerance RELTOL  
%                    for the desired accuracy of gamma-iteration. 
%                    The iterations are performed until the current 
%                    estimations of maximum gu and minimum gl of  
%                    the optimal distance go, gl =< go =< gu, satisfies
%                    gu-gl < RELTOL*(norm([G1 G2],inf)-norm(G2,inf)).
%                    (Default: RELTOL = 1.e-4)    
%

%  Author:      A. Varga, 15.01.2016.
%  Revision(s): A. Varga, 14.08.2016, 09.06.2017, 27.11.2017, 27.07.2018.
%
%  Method: The approach of [1] is used for the solution of the 2-block 
%          least distance problem.
%
%  References:
%  [1] C.-C. Chu, J. C. Doyle, and E. B. Lee
%      "The general distance problem in H-infinity optimal control theory",
%      Int. J. Control, vol 44, pp. 565-596, 1986.

narginchk(1,3)
nargoutchk(0,2)

if ~isa(sys1,'ss') 
   error('The input system SYS1 must be an SS object')
end    
[p,m] = size(sys1);

Ts = sys1.Ts;
discr = (Ts ~= 0); 

if nargin < 2
   sys2 = 0; 
end

if nargin < 3
   options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS',3)
end    

% decode options and set default values

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

% sub-optimality degree
if isfield(options,'gamma')
   gamma = options.gamma;
   if ~isempty(gamma)
      validateattributes(gamma, {'double'}, {'real','scalar', '>=', 0,'finite'},'','OPTIONS.gamma') 
   end
else
   gamma = [];  % optimal solution sought
end
subopt = ~isempty(gamma);

% relative tolerance for convergence of gamma-iteration
if isfield(options,'reltol')
   reltol = options.reltol;
   validateattributes(reltol, {'double'},{'real','scalar','>=',0,'<=',1},'','OPTIONS.reltol') 
else
   reltol = 1.e-4;
end


if isa(sys2,'double') 
    if ~isscalar(sys2)
        error('SYS2 must be a scalar')
    end
    m2 = sys2;
    validateattributes(m2, {'double'},{'integer','scalar','>=',0,'<=',m},'','OPTIONS.m2') 
    sys = sys1;
    m1 = m-m2;
else
   if ~isa(sys2,'ss')
      error('The input system SYS2 must be an SS object')
   end
   if (discr && (sys2.Ts <= 0)) || (~discr && (sys2.Ts > 0))
      error('The systems SYS1 and SYS2 must have the same type')
   end
   if discr && (Ts ~= sys2.Ts) 
      error('The systems SYS1 and SYS2 must have the same sampling period')
   end
   m1 = size(sys1,2);
   m2 = size(sys2,2);
   m = m1+m2;
   sys = gir([sys1,sys2],tol);
end    
    
% address constant case
if order(sys) == 0 || m1 == 0
   % solve a 2-block Parrott problem
   sysx = sys(:,1:m1);
   mindist = norm(sys(:,m1+1:end),inf);
   return
end


[a,b,c,~,e] = dssdata(sys);
if tol == 0
   % set tolerance 
   tol = size(a,1)*eps*max([norm(e,1) norm(a,1) norm(b,1),norm(c,inf)]);
end

gl = norm(sys(:,m1+1:end),inf);
if subopt && gl > gamma
    error(['OPTIONS.gamma must be chosen greater than ',num2str(gl)])
end

if m2 == 0 || gl <= sqrt(eps)
  % solve one block LDP
  if subopt 
     % solve sub-optimal Nehari problem
     [sysx,mindist] = gnehari(sys(:,1:m1),gamma,tol);
     if nargout > 1
        mindist = norm(sys(:,1:m1)-sysx,inf);
     end
  else
     % solve optimal Nehari problem
     [sysx,mindist] = gnehari(gir(sys(:,1:m1),tol),gamma,tol);
  end
else
  % solve two block LDP
  opt_gsdec = struct('tol',tol,'job','unstable');
  if subopt 
     % solve sub-optimal LDP
     gam = max(gamma,gl*1.01);
     opt_glsfg = struct('stab',true,'tol',tol);
     V = glsfg(sys(:,m1+1:end),gam,opt_glsfg);
     syse = V\sys(:,1:m1);
     hlu = ghanorm(gsdec(syse,opt_gsdec)');
     % solve the sub-optimal Nehari problem
     %[sysx,mindist] = gnehari(syse,hlu);
     if hlu < 1
        sysx =  gbalmr(V*gnehari(syse,hlu,tol),tol);
     else
        subopt = false;
     end
     if nargout > 1
        mindist = norm(sys-sysx*eye(m1,m),inf);
     end
  end
  if ~subopt
     % solve optimal LDP using gamma-iteration
     % initialize the gamma-iteration
     % compute upper bound and inital gap
     hlu = ghanorm(gsdec(sys(:,1:m1),opt_gsdec)');
     gu = norm([hlu;gl]); gl = max(gl,hlu);
     gap = max(1,gu-gl);
     %gu = norm(sys,inf); gap = max(1,gu-gl);
     gam = (gl+gu)/2; 

     % stabilize G2 only once
     opt_grcfid = struct('stab',false,'tol',tol);
     G2s = grcfid(sys(:,m1+1:end),opt_grcfid); 

     % perform the gamma-iteration
     opt_glsfg = struct('mininf',true,'tol',tol);
     iter = 0; g0 = gam; hlu = 1;
     while gu-gl > reltol*gap
        iter = iter+1; 
        % try to catch failure of spectral factorization if iter > 1
        % by reusing the last admissible value
        try 
           lastwarn(''); 
           V = glsfg(G2s,gam,opt_glsfg); 
           if isempty(lastwarn)
              g0 = gam;
              fail = false;
           else
              if iter > 1 
                 gam = g0;   % redo factorization with last value
                 V = glsfg(G2s,gam,opt_glsfg); 
                 fail = true;
                 warning('Desired accuracy of gamma-iteration not achievable')
              else
                 error('Failure of the minimum-phase spectral factorization')
              end
           end
        catch
           if iter > 1 
              gam = g0;   % redo factorization with last value
              V = glsfg(G2s,gam,opt_glsfg); 
              fail = true;
              warning('Desired accuracy of gamma-iteration not achievable')
           else
              error('Failure of the minimum-phase spectral factorization')
           end
        end
        % compute inv(V)*G1
        syse = V\sys(:,1:m1);
        hlu = ghanorm(gsdec(syse,opt_gsdec)');
        if fail, break, end
        if hlu < 1
           gu = gam;           % gamma > gamma_opt
        else
           gl = gam;           % gamma <= gamma_opt
        end   
        gam = (gl+gu)/2;
     end
     
     if iter 
        if hlu >= 1
           % recompute spectral factorization if gamma < gamma_opt
           V = glsfg(G2s,gu,opt_glsfg);  
           % compute inv(V)*G1
           syse = V\sys(:,1:m1);
           gam = gu;
        end
     else
        % compute spectral factorization for increased gamma
        gam = min([gu gam*(1+reltol)]);
        V = glsfg(G2s,gam,opt_glsfg);  
        % compute inv(V)*G1
        syse = V\sys(:,1:m1);    
     end

     % solve the Nehari problem and compute solution
     sysx =  gbalmr(V*gnehari(syse),tol); 
     mindist = gam;
  end
end

% end of GLINFLDP
end
