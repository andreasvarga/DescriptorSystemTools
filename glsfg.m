function sysf = glsfg(sys,gamma,options)
%GLSFG Left spectral factorization of gamma^2*I-G*G'.
%  SYSF = GLSFG(SYS,gamma) computes for the LTI descriptor system
%  SYS = (A-lambda*E,B,C,D) with the transfer-function matrix G(lambda), 
%  the minimum-phase right spectral factor SYSF with the transfer-function 
%  matrix F(lambda), such that
%           F*F' = gamma^2*I-G*G' ,
%  where gamma > norm(G,inf).
%
%  SYSF = GLSFG(SYS,gamma,OPTIONS) uses OPTIONS structure to specify
%  various user options, as follows:
%  OPTIONS.tol       - specifies the tolerance for singular values
%                      based rank determination of E
%                      (Default: tol = max(size(E))*eps(norm(E,1)))
%  OPTIONS.tolmin    - specifies the tolerance for the singular values
%                      based controllability tests
%                      (Default: tolmin = max(size(B))*eps(norm(B,1)))
%  OPTIONS.stabilize - stabilzation option
%                      true - replace SYS by SYSN = SYS*SYSM, where
%                             SYSM is an inner denominator factor of a
%                             stable RCF SYS = SYSN*INV(SYSM) (default)
%                      false - no stabilization performed
%
%  Note: The system SYS must not have poles on the imaginary-axis in the 
%  continuous-time case or on the unit circle in the discrete-time case.
%  If the stabilization is not performed, SYS must be stable.
%
%  See also GRSFG.

%  Copyright 2016-2018 A. Varga
%  Author:    A. Varga 22.01.2016.
%  Revisions: A. Varga 21.08.2016, 04.12.2016, 08.06.2017.
%
%  Method:  Extensions of the factorization approaches of [1] are used.
%
%  References:
%  [1] K. Zhou, J. C. Doyle, and K. Glover. 
%      Robust and Optimal Control. Prentice Hall, 1996.

narginchk(2,3)
nargoutchk(0,1)

% check for state-space system
if ~isa(sys,'ss')
   error('The input system SYS must be an SS object')
end

if nargin < 3
   options = struct('tol',0);
end

% decode options and set default values

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

% tolerance for controllability tests
if isfield(options,'tolmin')
   tolmin = options.tolmin;
   validateattributes(tolmin, {'double'},{'real','scalar','>=',0},'','OPTIONS.tolmin') 
else
   tolmin = 0;
end

% stabilization option 
if isfield(options,'stabilize')
   stabilize = options.stabilize;
   validateattributes(stabilize, {'logical'},{'binary'},'','OPTIONS.stabilize') 
else
   stabilize = true;
end
 
% Perform right coprime factorization with inner denominator for
% unstable systems
if stabilize 
   % Compute the RCF with inner denominator
   options_grcfid = struct('tol',tol,'tolmin',tolmin,'mininf',true);
   sys = grcfid(sys,options_grcfid); 
end

[a,b,c,d,e,Ts] = dssdata(sys);
discr = (Ts ~= 0);

if ~isequal(e,eye(size(e))) && rcond(e) < sqrt(eps)
   sys = gss2ss(sys,tol); 
   [a,b,c,d,e,Ts] = dssdata(sys);
   if rcond(e) < sqrt(eps)
      error('Improper input system SYS')
   end
end

% Compute the stabilizing solution of the corresponding Riccati equation
[p,n] = size(c);
r = gamma*gamma*eye(p)-d*d.'; 

if n 
   if discr
      %[xric,kric] = gdare1(a.',c.',b*b.',-r,b*d.',e.',nochecks); kric = -kric;
      [xric,~,kric] = dare(a.',c.',b*b.',-r,b*d.',e.'); 
      [v,rdiag] = schur(xric);
      rdiag = sqrt(max(diag(rdiag),0));
      [u,s,~] = svd([d bsxfun(@times,c*v,rdiag.')]);
      % r = r - c*xric*c.'; 
   else
      %[~,kric] = gcare1(a.',c.',b*b.',-r,b*d.',e.',nochecks); kric = -kric;
      [~,~,kric] = care(a.',c.',b*b.',-r,b*d.',e.'); 
      [u,s,~] = svd(d);
   end
else
   [u,s,~] = svd(d);
   if ~isempty(s) && s(1,1) > abs(gamma)
      error('The condition gamma > norm(SYS) is not fulfilled')
   end
   kric = zeros(p,n);
end

% compute square-root factor
if min(size(s)) == 1
   s = s(1,1);
else
   s = diag(s); 
end
rsq = [s;zeros(p-length(s),1)].';
rsqrt = bsxfun(@times,u,sqrt(max(gamma*gamma-rsq.^2,0))); 
%rsqrt = sqrtm(r);

sysf = dss(a,kric.'*rsqrt,c,rsqrt,e,Ts);

% end GLSFG
end