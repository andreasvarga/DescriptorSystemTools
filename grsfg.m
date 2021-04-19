function sysf = grsfg(sys,gamma,options)
%GRSFG Right spectral factorization of gamma^2*I-G'*G.
%  SYSF = GRSFG(SYS,gamma) computes for the LTI descriptor system
%  SYS = (A-lambda*E,B,C,D) with the transfer-function matrix G(lambda), 
%  the minimum-phase right spectral factor SYSF with the transfer-function 
%  matrix F(lambda), such that
%           F'*F = gamma^2*I-G'*G ,
%  where gamma > norm(G,inf).
%
%  SYSF = GRSFG(SYS,gamma,OPTIONS) uses OPTIONS structure to specify
%  various user options, as follows:
%  OPTIONS.tol       - specifies the tolerance for singular values
%                      based rank determination of E
%                      (Default: tol = max(size(E))*eps(norm(E,1)))
%  OPTIONS.tolmin    - specifies the tolerance for the singular values
%                      based observability tests
%                      (Default: tolmin = max(size(C))*eps(norm(C,inf)))
%  OPTIONS.stabilize - stabilization option
%                      true - replace SYS by SYSN = SYSM*SYS, where
%                             SYSM is an inner denominator factor of a
%                             stable LCF SYS = INV(SYSM)*SYSN (default)
%                      false - no stabilization performed
%
%  Note: The system SYS must not have poles on the imaginary-axis in the 
%  continuous-time case or on the unit circle in the discrete-time case.
%  If the stabilization is not performed, SYS must be stable.
%
%  See also GLSFG.

%  Copyright 2016-2018 A. Varga
%  Author:    A. Varga 22.01.2016.
%  Revisions: A. Varga 14.09.2016, 07.06.2017.
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


% Perform left coprime factorization with inner denominator for
% unstable systems
if stabilize 
   % Compute the LCF with inner denominator
   options_glcfid = struct('tol',tol,'tolmin',tolmin,'mininf',true);
   sys = glcfid(sys,options_glcfid); 
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
[n,m] = size(b);
r = gamma*gamma*eye(m)-d.'*d; 

if n 
   if discr
      %[xric,fric] = gdare1(a,b,c.'*c,-r,c.'*d,e,true); fric = -fric;
      [xric,~,fric] = dare(a,b,c.'*c,-r,c.'*d,e); 
      %r = r-b.'*xric*b;
      [u,rdiag] = schur(xric);
      rdiag = sqrt(max(diag(rdiag),0));
      [~,s,v] = svd([d;bsxfun(@times,u.'*b,rdiag)]);
   else
      %[~,fric] = gcare1(a,b,c.'*c,-r,c.'*d,e,true); fric = -fric;
      [~,~,fric] = care(a,b,c.'*c,-r,c.'*d,e); 
      [~,s,v] = svd(d);
   end
else
   [~,s,v] = svd(d);
   if ~isempty(s) && s(1,1) > abs(gamma)
      error('The condition gamma > norm(SYS) is not fulfilled')
   end
   fric = zeros(m,n);
end

% compute square-root factor
if min(size(s)) == 1
   s = s(1,1);
else
   s = diag(s); 
end
%rsq = [diag(s);zeros(max(size(d))-size(d,1),1)];
rsq = [s;zeros(m-length(s),1)];
rsqrt = bsxfun(@times,v.',sqrt(max(gamma*gamma-rsq.^2,0))); 
%rsqrt = sqrtm(r);

% assemble the spectral factor
sysf = dss(a,b,rsqrt*fric,rsqrt,e,Ts);

% end GRSFG
end