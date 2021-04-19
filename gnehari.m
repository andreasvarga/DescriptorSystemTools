function [sysx,s1] = gnehari(sys,gamma,tol)
%GNEHARI Generalized Nehari approximation.
%        [SYSX,S1] = GNEHARI(SYS) computes for the descriptor system SYS
%        the optimal stable Nehari solution SYSX using
%        zero order antistable Hankel-norm optimal approximation.
%        S1 is the L-infinity norm of the optimal error system SYS-SYSX.
%
%        [SYSX,S1] = GNEHARI(SYS,GAMMA) computes the stable
%        Nehari suboptimal solution using suboptimal zero order Hankel norm
%        approximation. GAMMA must be choosen greater than, S1, the largest
%        Hankel singular value of the conjugate of the antistable part.
%
%        [SYSX,S1] = GNEHARI(SYS,GAMMA,TOL) uses tolerance TOL for rank 
%        computations. GAMMA = [] must be used to compute an optimal Nehari
%        approximation. 
%
%        Note: SYS must not have eigenvalues on the boundary of the 
%        stability domain (i.e., on the imaginary axis, for a 
%        continuous-time system or on the unit circle, for a 
%        discrete-time system.

%  Author:      A. Varga, 18.01.2016.
%  Revision(s): A. Varga, 09.06.2017, 26.11.2017, 16.02.2018, 25.07.2018,
%                         09.04.2021. 

%  Method:  The Hankel-norm approximation methods of [1] (see also [2]), 
%  with extensions for descriptor systems, is used for the approximation
%  of the unstable part.   
%
%  References:
%  [1] K. Glover. 
%      All optimal Hankel-norm approximations of linear
%      multivariable systems and their L-infinity error bounds,
%      Int. J. Control, vol. 39, pp. 1115-1193, 1984.
%
%  [2] M. G. Safonov, R. Y. Chiang, and D. J. N. Limebeer. 
%      Optimal Hankel model reduction for nonminimal systems. 
%      IEEE Trans. Automat. Control, vol. 35, pp. 496–502, 1990.

narginchk(1,3)
nargoutchk(0,2)

if ~isa(sys,'ss')
   error('The input system SYS must be a state-space system')
end

% set default inputs
if nargin >= 2
   if isempty(gamma) 
      subopt_flg = false;
   else
      validateattributes(gamma, {'double'}, {'real','scalar', '>', 0,'finite'},'','GAMMA',2) 
      subopt_flg = true;
   end
else
   subopt_flg = false;
end

if nargin == 3
   validateattributes(tol, {'double'},{'real','scalar','>=',0,'<',1},'','TOL',3) 
else
   tol = 0;
end

discr = (sys.Ts ~= 0);

% set stability margin to separate stable and unstable parts
s2eps = sqrt(eps);
if discr 
   smarg = 1-eps^0.33;
else
   smarg = -s2eps;
end

% Stable/antistable separation; possible eigenvalues on the boundary of 
% the stability domain are included in the unstable part 
opt_gsdec = struct('tol',tol,'smarg',smarg,'job','stable');
[sysx,sysu] = gsdec(sys,opt_gsdec); 

% only proper systems are handled in the continuous-time case
if ~discr && rcond(sysu.e) < s2eps
   sysu = gss2ss(sysu,s2eps);
   if rcond(sysu.e) < s2eps
       error('Improper continuous-time system SYS')
   else
       sysx.d = sysx.d+sysu.d;
       sysu.d = zeros(size(sysu.d));
   end
end

% exit if the system has only stable part
if  size(sysu,'order') == 0
    s1 = 0;
    return
end

[a,b,c,d,e] = dssdata(sysu);
ev = eig(a,e); 
if discr
   if min(abs(ev)) <= 1+s2eps
       error('The system SYS has possibly poles on the unit circle')
   end
else
    if min(real(ev)) <= s2eps
       error('The system SYS has possibly poles on the imaginary axis')
    end
end

if discr
   % use the bilinear transformation to compute an equivalent unstable, 
   % possibly non-minimal, continuous-time system  
   sqrt2 = sqrt(2);
   esave = e;
   e = a+e; a = a-esave; 
   c = c/e; 
   d = d-c*b; 
   b = sqrt2*b; c = sqrt2*c*esave;
end

% compute a balanced minimal realization of antistable part and the
% Hankel singular values
[sysr,hs] = gbalmr(dss(-a,-b,c,d,e),s2eps,'balance');
[a,b,c,d] = ssdata(sysr);


na = size(a,1); 
[p,m] = size(d);

s1 = hs(1);  % expected L-inf norm of the optimal approximation error 

% Determine the type of required approximation
if subopt_flg 
  if abs(gamma-s1) <= s1*s2eps
     subopt_flg = false; % enforce optimal approximation
  elseif gamma < s1
    error('Suboptimal antistable Hankel approximation cannot be computed')
  end
else
  gamma = hs(1);
end

if subopt_flg
   sk = gamma; r = 0; 
else
   sk = s1; ep = s1*s2eps;
   if na
      r = 1;
      for i = 2:na
          if abs(sk-hs(i)) > ep, break, end
          r = r+1;
      end
   else
      r = 0;
   end
end
%          ^
%  Compute F(s)
ns = na-r; sk2 = sk*sk; i1 = r+1:na; 

% make system square by padding with zeros rows or column 
if m < p, 
    b = [b zeros(na,p-m)]; 
elseif m > p
    c = [c; zeros(m-p,na)]; 
end
if subopt_flg
   u = zeros(max(p,m));
   b1 = b; c1 = c;
else
   i2 = 1:r;
   % compute orthogonal u such that u*b(i2,:)'+c(:,i2) = 0
   [q1,~,~] = svd(b(i2,:)');
   [q2,~,~] = svd(c(:,i2));
   u = -q2*q1'; 
   b1 = b(i1,:); c1 = c(:,i1);
end
jm = 1:m; jp = 1:p;
if ns
   % optimal solution
   e = diag(hs(i1).*hs(i1)-sk2);
   a = -((sk2*a(i1,i1)'...
            +diag(hs(i1))*a(i1,i1)*diag(hs(i1))-sk*c1'*u*b1'));
   b = -(diag(hs(i1))*b1(:,jm)+sk*c1'*u(:,jm));
   c = c1(jp,:)*diag(hs(i1))+sk*u(jp,:)*b1';
else
   a = []; e = []; b = zeros(0,m); c = zeros(p,0); 
end
if na 
   d = d-sk*u(jp,jm);
end

if discr 
   % compute the equivalent discrete-time solution  
   esave = e;
   e = e-a; a = a+esave; 
   c = c/e; 
   d = d+c*b; 
   b = sqrt2*b; c = sqrt2*c*esave;
end

sysx = sysx + dss(a,b,c,d,e,sys.Ts);

% end GNEHARI
end

