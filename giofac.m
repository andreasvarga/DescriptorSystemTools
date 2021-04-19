function [sysi,syso,info] = giofac(sys,options)
%GIOFAC Inner-outer/QR-like factorization.
% [SYSI,SYSO,INFO] = GIOFAC(SYS,OPTIONS) calculates, for a descriptor 
% system SYS = (A-lambda*E,B,C,D) with the transfer-function matrix 
% G(lambda), the square inner factor SYSI = (Ai-lambda*Ei,Bi,Ci,Di) with 
% the transfer function matrix Gi(lambda) and the minimum-phase quasi-outer
% or the full row rank factor SYSO = (A-lambda*E,B,Co,Do) with the 
% transfer function matrix Go(lambda) such that 
%
%                 G(lambda) = Gi(:,1:r)(lambda)*Go(lambda) ,  (*)
%
% where r is the normal rank of G(lambda). The resulting proper and stable 
% inner factor satisfies Gi'(lambda)*Gi(lambda) = I. The resulting factor 
% Go(lambda) has full row rank r. Depending on the selected factorization,
% Go(lambda) is either minimum phase, excepting possibly zeros on the 
% boundary of the appropriate stability domain, or contains all zeros 
% of G(lambda), in which case (*) is the extended QR-like factorization of
% G(lambda).
% Note: SYSO may generally contain free inner factors which can be 
%       eliminated by removing the finite unobservable eigenvalues. 
%
% The OPTIONS structure allows to specify various user options, 
% as follows:
% OPTIONS.tol      - tolerance for rank determinations and observability 
%                    tests (Default: internally computed)
% OPTIONS.minphase - option to compute a minimum-phase quasi-outer factor:
%                    true  - compute a minimum phase quasi-outer factor, 
%                            with all zeros stable, excepting possibly 
%                            zeros on the boundary of the appropriate 
%                            stability domain (Default);
%                    false - compute a full row rank factor, which 
%                            includes all zeros of G(lambda).  
% OPTIONS.offset   - stability boundary offset OFFSET, to specify the 
%                    marginally stable finite zeros as follows: 
%                    in the continuous-time case are the finite zeros 
%                    having real parts in the interval [-OFFSET, OFFSET], 
%                    while in the discrete-time case having moduli in the 
%                    interval [1-OFFSET,1+OFFSET] 
%                    (Default: OFFSET = sqrt(eps)).
%                    
% OPTIONS.balance  - specifies the balancing option for the Riccati 
%                    equation solvers: 
%                    true  - use balancing (Default);
%                    false - disable balancing.  
%   
% INFO is a structure containing additional information, as follows: 
% INFO.nrank  - normal rank r of the transfer function matrix G(lambda);
% INFO.nfuz   - number of finite unstable zeros of SYSO; these are the 
%               marginally stable finite zeros of SYS lying on the boundary
%               of the stability region;
% INFO.niuz   - number of infinite zeros of SYSO; these are the infinite  
%               zeros of SYS in the continuous-time case and 0 in the 
%               discrete-time case
% INFO.ricrez - diagnosis flag, as provided provided by the generalized 
%               Riccati equation solvers care and dare; if non-negative, 
%               this value represents the Frobenius norm of relative 
%               residual of the Riccati equation, while a negative value 
%               indicates failure of solving the Riccati equation. 
%
% See also GOIFAC.

% Copyright 2015-2018 A. Varga
% Author:      A. Varga, 24.12.2015.
% Revision(s): A. Varga, 09.11.2016, 26.07.2017, 21.08.2018, 09.09.2018,
%                        09.06.2019.
%
% References:
% [1] C. Oara and A. Varga.
%     Computation of the general inner-outer and spectral factorizations.
%     IEEE Trans. Autom. Control, vol. 45, pp. 2307--2325, 2000.
% [2] C. Oara.
%     Constructive solutions to spectral and inner–outer factorizations 
%     respect to the disk. Automatica,  41, pp. 1855–1866, 2005. 

narginchk(1,2)
nargoutchk(0,3)

if ~isa(sys,'ss')
   error('The input system SYS must be a state space system')
end

if nargin < 2
   options = struct('tol',0,'minphase',true);
else
   if isa(options,'double')
      % support for old syntax (undocummented)
      validateattributes(options, {'double'},{'real','scalar','>=',0},'','TOL') 
      options = struct('tol',options,'minphase',true);
   else    
      validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
   end
end

% decode options

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

% minimum phase option
if isfield(options,'minphase')
   minphase = options.minphase;
   validateattributes(minphase, {'logical'},{'binary'},'','OPTIONS.minphase')
   if minphase
      zerosopt = 'unstable';
   else
      zerosopt = 'none';
   end
else
   zerosopt = 'unstable';
end

% offset for the stability region boundary
if isfield(options,'offset')
   offset = options.offset;
   validateattributes(offset, {'double'},{'real','scalar','>',0,'<',1},'','OPTIONS.offset') 
else
   offset = sqrt(eps);
end

% balancing option for Riccati equation solvers
if isfield(options,'balance')
   balopt = options.balance;
   validateattributes(balopt, {'logical'},{'binary'},'','OPTIONS.balance')
   if balopt
      balance = 'balance';
   else
      balance = 'nobalance';
   end    
else
   balance = 'balance';
end

%  Reduce the system matrix pencil to the special Kronecker-like form
%
%                 ( Arg-lambda*Erg     *            *   *   )
%                 (   0            Abl-lambda*Ebl  Bbl  *   )
%  At-lambda*Et = (   0                0            0   Bn  )
%                 (-----------------------------------------)
%                 (   0               Cbl          Dbl  *   )
%
% where the subpencil
%                           ( Abl-lambda*Ebl  Bbl )
%                           (      Cbl        Ddl )
%
% is full column rank, the pair (Abl-lambda*Ebl,Bbl) is stabilizable, and
% Abl-lambda*Ebl contains either no zeros of SYS if 
% OPTIONS.minphase = false, or contains all unstable zeros of SYS if 
% OPTIONS.minphase = true. 


% enforce controllability
sys = gir(sys,tol,'contr'); 
[At,Et,dimsc,~,Z] = gsklf(sys,tol,zerosopt,offset,'noQ');

%  Select matrices for the extracted full column rank factor
%
[a,b,~,d,e,Ts] = dssdata(sys);
discr = (Ts ~= 0);
[n,m] = size(b); p = size(d,1);
mc    = dimsc(1); 
nric  = dimsc(2);
mric  = dimsc(3); 
nsinf = dimsc(4);
nm    = n+m;
nc    = n-nric-nsinf; 

ia = nc+1:nc+nric; ja = mc+1:mc+nric; jb = mc+nric+1:mc+nric+mric;
ic = n+1:n+p; 
A = At(ia,ja); E = Et(ia,ja); B = At(ia,jb); C = At(ic,ja); D = At(ic,jb);

if discr
   %[X,F,ricrez] = gdare1(A,B,C.'*C,D.'*D,C.'*D,E,true); 
   [X,~,F,ricrez] = dare(A,B,C.'*C,D.'*D,C.'*D,E,balance); F = -F;
   if ricrez < 0 || isnan(ricrez)
      if ricrez == -1
         error('Solution of the DARE failed: Symplectic matrix has eigenvalues on the unit circle')
      else
         error('Solution of the DARE failed: no finite stabilizing solution exists')
      end
   end
else
   %[X1,F1,ricrez] = gcare1(A,B,C.'*C,D.'*D,C.'*D,E,true);
   [X,~,F,ricrez] = care(A,B,C.'*C,D.'*D,C.'*D,E,balance); F = -F;
   if ricrez < 0 || isnan(ricrez)
      if ricrez == -1 
         error('Solution of the CARE failed: Hamiltonian has jw-axis eigenvalues')
      else
         error('Solution of the CARE failed: no finite stabilizing solution exists')
      end
   end
end

% assemble the inner and outer factors
if discr
   mm1 = p-mric;
   H = chol(D'*D+B'*X*B); 
   if mm1
      [~,s,V]=svd([A'*X C'; B'*X D'],0);
      if nric+mric > 1, s = diag(s);
      elseif nric+mric == 1, s = s(1);
      else, s = 0;
      end
      tolr = max(mric,nric) * max(s) * eps;
      r = sum(s > tolr); 
      Y = V(1:nric,r+1:nric+p);
      W = V(nric+1:nric+p,r+1:nric+p);
      T = chol(W'*W+Y'*X*Y);
      Y = Y/T; W = W/T; 
   else
      Y = []; W = [];
   end
else
   [Q,H]=qr(D); H = H(1:mric,:); 
   W = Q(:,mric+1:p); 
   if rank(X) < nric
      Y = -pinv(X)*((E')\C'*W); 
   else
      Y = -(E'*X)\C'*W;
   end
end
 
agi = A+B*F; bgi = [B/H Y]; cgi = C+D*F; dgi = [D/H W]; egi = E; 

CDt = [zeros(mric,mc) -H*F H zeros(mric,nsinf)]*Z'; 

% convert Gi to standard state-space if the original system 
% was in state-space form
if isequal(e,eye(n))
   [agi,egi,bgi,cgi,dgi] = sl_gminr(2,agi,egi,bgi,cgi,dgi,tol,1);
end

sysi = dss(agi,bgi,cgi,dgi,egi,Ts);
syso = dss(a,b,CDt(:,1:n),CDt(:,n+1:nm),e,Ts);
info = struct('nrank',mric,'nfuz',dimsc(5),'niuz',dimsc(6),'ricrez',ricrez);

% end GIOFAC
end
