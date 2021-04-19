function [sysx,info] = glasol(sysg,sysf,options)
%GLASOL Approximate solution of the linear rational matrix equation X*G = F.
%   [SYSX,INFO] = GLASOL(SYSG,SYSF,OPTIONS) returns a state-space  
%   realization SYSX of the approximate solution X(lambda) of the  
%   linear rational equation
%
%      X(lambda)*G(lambda) = F(lambda)                       (1)
%
%   which minimizes the error norm ||X(lambda)*G(lambda) - F(lambda)||.
%   G(lambda) and F(lambda) are the transfer function matrices of 
%   the state-space realizations in SYSG and SYSF, respectively, and SYSF 
%   is assumed stable. The resulting X(lambda) has all poles stable or  
%   lying on the boundary of the stability domain.
%
%   The OPTIONS structure specifies various user options, as follows:
%   OPTIONS.tol     - specifies the relative tolerance for rank 
%                     determinations. (Default: internally computed)
%   OPTIONS.offset  - stability boundary offset OFFSET, to specify the 
%                     marginally stable finite zeros of SYSG as follows: 
%                     in the continuous-time case are the finite zeros 
%                     having real parts in the interval [-OFFSET, OFFSET], 
%                     while in the discrete-time case having moduli in the 
%                     interval [1-OFFSET,1+OFFSET] 
%                     (Default: OFFSET = sqrt(eps)).
%   OPTIONS.sdeg    - specifies a prescribed stability degree for the free
%                     eigenvalues of SYSX 
%                     (Default: [], in which case:
%                                   -0.05 in the continuous-time case;
%                                    0.95 in the discrete-time case)
%   OPTIONS.poles   - specifies a complex conjugated set of desired poles
%                     to be assigned as free eigenvalues of SYSX 
%                     (Default: []) 
%   OPTIONS.mindeg  - minimum degree solution 
%                     true  - determine, if possible, a  minimum order 
%                             solution with all poles stable or lying on 
%                             the boundary of the stability domain 
%                     false - determine a particular solution which has 
%                             possibly non-minimal order and all its poles 
%                             stable or lying on the boundary of the 
%                             stability domain (default) 
%   OPTIONS.H2sol   - option to compute a H2-norm optimal solution 
%                     true  - perform a H2-norm optimal solution 
%                     false - perform a Hinf-norm optimal solution (default) 
%   OPTIONS.gamma   - specifies gamma > 0, the desired sub-optimality 
%                     degree to solve the gamma-suboptimal model-matching
%                     problem 
%                        ||X(lambda)*G(lambda)-F(lambda)||_inf <= gamma
%                     (Default: [] - the optimal error minimization problem 
%                                    is solved)
%   OPTIONS.reltol  - specifies the relative tolerance RELTOL  
%                     for the desired accuracy of the gamma-iteration. 
%                     The iterations are performed until the current 
%                     estimations of maximum gu and minimum gl of  
%                     the optimal value go, gl =< go =< gu, satisfies
%                     gu-gl < RELTOL*gapini, where gapini is an initial
%                     estimation of the distance gap. 
%                     (Default: RELTOL = 1.e-4)  
%
%   The resulting INFO structure contains additional information:
%   INFO.nrank is the normal rank of G(lambda); 
%   INFO.mindist is the achieved approximation error norm; 
%   INFO.nonstandard is true for a non-standard problem, with G(lambda) 
%       having zeros on the boundary of the stability domain, and false for
%       a standard problem, when G(lambda) has no zeros on the boundary 
%       of the stability domain.  
% 
%   [SYSX,INFO] = GLASOL(SYSGF,MF,OPTIONS) uses the compound realization 
%   SYSGF of [G(lambda); F(lambda)], where F(lambda) has MF rows.
%
%   See also GRASOL.

%  Copyright 2018 A. Varga
%  Author:       A. Varga, 24.07.2018.
%  Revision(s):  A. Varga, 21.08.2018, 18.09.2018.
%
%  Method: An extension of the approach of [1] to descriptor systems 
%  is implemented.
%
%  References:
%  [1]  B. A. Francis. A Course in H-infinity Theory, 
%       Vol. 88, of Lecture Notes in Control and Information Sciences. 
%       Springer-Verlag, New York, 1987.

narginchk(2,3)
nargoutchk(0,2)

% decode options
if nargin < 3
   options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

if ~isa(sysg,'ss')
   error('The input system SYSG must be an SS object')
end  
Ts = sysg.Ts;
[p,m] = size(sysg);
discr = (Ts ~= 0); 
if discr
    smax = 1-sqrt(eps); smin = 0;
else
    smax = -sqrt(eps);  smin = -inf;
end

if ~isa(sysf,'ss')
   if ~isa(sysf,'double') 
      error('SYSF must be an SS object or a positive integer')
   else
      validateattributes(sysf,{'double'},{'integer','scalar','>=',0,'<=',size(sysg,1)},'','MF')
      mf = sysf; 
      p = p-mf;
      sysf = gir(sysg(p+1:p+mf,:),tol);
      sysg = gir(sysg(1:p,:),tol);
   end
else
   if (discr && (sysf.Ts == 0)) || (~discr && (sysf.Ts ~= 0))
      error('The systems SYSG and SYSF must have the same type')
   end
   if discr && (Ts ~= sysf.Ts) 
      error('The systems SYSG and SYSF must have the same sampling period')
   end
   if size(sysg,2) ~= size(sysf,2) 
      error('The systems SYSG and SYSF must have the same number of inputs')
   end  
end  

% decode the rest of options

% offset for the stability region boundary
if isfield(options,'offset')
   offset = options.offset;
   validateattributes(offset, {'double'},{'real','scalar','>',0,'<',1},'','OPTIONS.offset') 
else
   offset = sqrt(eps);
end

if ~isempty(sysf)
   % check stability of SYSF
   [~,infop] = gpole(sysf,tol,offset);
   if ~infop.stable
      %error('SYSF must be a stable system') 
   end
end

% desired stability degree
if isfield(options,'sdeg')
   sdeg = options.sdeg;
   if ~isempty(sdeg)
      validateattributes(sdeg, {'double'},{'real','scalar','<=',smax,'>=',smin},'','OPTIONS.sdeg') 
   end
else
   sdeg = [];
end
if isempty(sdeg)
   if discr 
      sdeg = 0.95;
   else
      sdeg = -0.05;
   end
end

% desired poles
if isfield(options,'poles')
   poles = options.poles;
   if ~isempty(poles)
      validateattributes(poles, {'double'},{'vector'},'','OPTIONS.poles')
   end
else
   poles = [];
end
if ~isempty(poles)
    t = poles;
    if (min(size(t)) ~= 1 || ~isequal(sort(t(imag(t)>0)),sort(conj(t(imag(t)<0)))) )
        error('OPTIONS.poles must be a self-conjugated complex vector')
    end 
end    

% option to perform a H2-norm optimal solution
if isfield(options,'H2sol')
   H2sol = options.H2sol; 
   validateattributes(H2sol, {'logical'},{'binary'},'','OPTIONS.H2sol') 
else
   H2sol = false;
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

% relative tolerance for convergence of gamma-iteration
if isfield(options,'reltol')
   reltol = options.reltol;
   validateattributes(reltol, {'double'},{'real','scalar','>=',0,'<=',1},'','OPTIONS.reltol') 
else
   reltol = 1.e-4;
end

% minimum degree solution option
if isfield(options,'mindeg')
   mindeg = options.mindeg;
   validateattributes(mindeg, {'logical'},{'binary'},'','OPTIONS.mindeg') 
else
   mindeg = false;
end

% compute the extended quasi-co-outer-co-inner factorization
[Gi,Go,info1] = goifac(sysg,struct('tol',tol,'offset',offset)); 
ro = info1.nrank; 

% detect nonstandard problem 
%nonstandard = order(gcrange(Go,struct('tol',tol,'zeros','s-unstable'))) > 0;
nonstandard = (info1.nfuz + info1.niuz > 0);

% define and solve the LDP
F = sysf*Gi';
if H2sol
   % set tolerance for feedthrough matrix
   if tol
      told = tol;
   else
      told = 1.e4*eps;
   end
   if ~discr && norm(F.d(:,ro+1:end)) < told
      F.d(:,ro+1:end) = zeros(size(sysf,1),m-ro);
   end
   % solve the H2-LDP min||[ F1-Xt; F2 ] ||_2
   [Xt,Xtu] = gsdec(F(:,1:ro),struct('tol',tol,'job','stable'));
   gopt = norm(glcfid(gir([Xtu F(:,ro+1:end)],tol)),2); 
else
   % solve the H_inf-LDP min ||[ F1-Xt; F2 ] ||_inf
   opts = struct('tol',tol,'reltol',reltol,'gamma',gamma);
   [Xt,gopt] = glinfldp(F,m-ro,opts); 
end

if ro == p
   sysx = gir(Xt/Go,tol);
else
   opts_glsol = struct('tol',tol,'sdeg',sdeg,'poles',poles,'mindeg',mindeg);
   sysx = glsol(Go,Xt,opts_glsol); 
   if mindeg && ~isempty(sysx)
      ev = gpole(sysx,tol);
      if ~discr, ev = ev(isfinite(ev)); end 
      if ~isempty(ev) 
         if (~discr && max(real(ev)) > sqrt(eps)) || ...
            (discr && max(abs(ev)) > 1+sqrt(eps))
            warning('Least order solution has strictly unstable poles: a non-minimal solution computed instead') 
            opts_glsol.mindeg = false;
            sysx = glsol(Go,Xt,opts_glsol); 
         end
      end
   end
end

info = struct('nrank',ro,'mindist',gopt,'nonstandard',nonstandard);

% end GLASOL
end
