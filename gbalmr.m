function [sysr,hs] = gbalmr(sys,tol,balance)
%GBALMR Balancing-based model reduction of a LTI descriptor system.
%       [SYSR,HS] = GBALMR(SYS)  calculates for a stable LTI
%       descriptor system SYS = (A-lambda*E,B,C,D) a minimal realization 
%       SYSR = (Ar-lambda*Er,Br,Cr,D). 
%       HS contains the decreasingly ordered Hankel singular values of SYS. 
%
%       [SYSR,HS] = GBALMR(SYS,TOL) uses the tolerance TOL for minimal
%       realization or order reduction. The order of the realization of 
%       SYSR is the number of Hankel singular values greater than 
%       TOL*max(1,HS(1)). (Default: TOL = sqrt(eps)). 
%
%       [SYSR,HS] = GBALMR(SYS,TOL,'balance') calculates for the 
%       stable LTI descriptor system SYS = (A-lambda*E,B,C,D) a standard 
%       LTI balanced minimal realization SYSR = (Ar-lambda*I,Br,Cr,D). 
%

%  Author:      A. Varga, 17.01.2016.
%  Revision(s): A. Varga, 07.11.2016, 23.07.2017, 26.11.2017.
%
%  Method:  For the order reduction of a standard system, the 
%  balancing-free method of [1] or the balancing-based method of [2] are 
%  used. For a descriptor system the balancing related order reduction 
%  methods of [3] are used.
%
%  References:
%  [1] A. Varga
%      Efficient minimal realization procedure based on balancing.
%      In A. El Moudni, P. Borne, and S.G. Tzafestas (Eds.), 
%      Prepr. of the IMACS Symp. on Modelling and Control of Technological 
%      Systems, Lille, France, vol. 2, pp.42-47, 1991.
%  [2] M. S. Tombs and I. Postlethwaite. 
%      Truncated balanced realization of a stable non-minimal state-space 
%      system. Int. J. Control, vol. 46, pp. 1319–1330, 1987.
%  [3] T. Stykel. 
%      Gramian based model reduction for descriptor systems. 
%      Mathematics of Control, Signals, and Systems, 16:297–319, 2004.

narginchk(1,3)
nargoutchk(0,2)

if ~isa(sys,'ss')
   error('The input system SYS must be an SS object')
end

if nargin < 2 || isempty(tol) 
   tol = sqrt(eps); 
else
   validateattributes(tol,{'double'},{'real','scalar', '>=', 0, '<=',1},'','TOL') 
end
if tol == 0, tol = sqrt(eps); end


if nargin < 3
   bal = false;
else
   validateattributes(balance,{'char'},{'nonempty'},'','BALANCE') 
   bal = strcmp(balance,'balance');
   if ~bal
      error('Improper balancing option')
   end
end
 
[a,b,c,d,e,Ts] = dssdata(sys);
n = size(a,1);
standsys = isequal(e,eye(n));  
discr = (Ts ~= 0);

if size(a,1) == 0,
    sysr = sys; 
    hs = [];
    return
end


if  standsys
    % reduce the system to real Schur coordinate form
    [Q,as] = schur(a,'real');
    ev = ordeig(as); 
    % check stability
    if (discr && max(abs(ev)) >= 1-sqrt(eps)) || ...
       (~discr && max(real(ev)) >= -sqrt(eps))
          error('The system SYS is unstable')
    end
    bs = Q'*b; cs = c*Q; es = e;
 else
    % reduce the system to generalized real Schur coordinate form
    [as,es,Q,Z] = qz(a,e,'real');
    if rcond(es) < eps
       error('The system SYS = (A-lambda*E,B,C,D) has singular E matrix')
    end
    ev = ordeig(as,es); 
    if (discr && max(abs(ev)) >= 1-sqrt(eps)) || ...
       (~discr && max(real(ev)) >= -sqrt(eps))
          error('The system SYS is unstable')
    end
    bs = Q*b; cs = c*Z; 
end

if discr
   flag = [ 1 1]; trans = 1;
   S = sl_glme(4,as,es,bs,flag,trans);
   trans = 0;
   R = sl_glme(4,as,es,cs,flag,trans);
else
   flag = [ 0 1]; trans = 1;
   S = sl_glme(4,as,es,bs,flag,trans);
   trans = 0;
   R = sl_glme(4,as,es,cs,flag,trans);
end

[u,sig,v]= svd(R*es*S); hs = diag(sig); 

% Form the truncated transformation matrices by selecting the
% Hankel singular values greater than TOL*hs(1)
if isempty(hs)
    nr = 0;
else
    nr = length(find(hs>tol*max(1,hs(1))));
end

if nr == 0
   sysr = ss(d);
   sysr.Ts = Ts;
   return
end

ind = 1:nr;

if bal
   % apply balancing formulas
   hsi2 = 1./sqrt(hs(ind)); 
   % efficient versions of L = diag(hsi2)*(u(:,ind)'*R) and
   % T = S*v(:,ind)*diag(hsi2)
   L = bsxfun(@times,u(:,ind)'*R,hsi2); 
   T = bsxfun(@times,S*v(:,ind),hsi2');
   % simpler versions using implicit expansion introduced in MATLAB R2016b 
   % L = (u(:,ind)'*R).*hsi2';  
   % T = S*v(:,ind).*hsi2; 
   % build the minimal balanced system
   sysr = ss(L*as*T,L*bs,cs*T,d,Ts);   
else
   if nr < n 
      % apply balancing-free formulas
      [L,~] = qr(R'*u(:,ind),0); 
      [T,~] = qr(S*v(:,ind),0);
      L = L';
      % build the minimal system
      sysr = dss(L*as*T,L*bs,cs*T,d,L*es*T,Ts);  
   else
      sysr = sys;  % keep original system if order is preserved
   end
end

% end GBALMR
end
