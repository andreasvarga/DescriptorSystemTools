function [sys1,sys2,q,z] = gsdec(sys,options)
%GSDEC  Generalized additive spectral decompositions.
%  [SYS1,SYS2] = GSDEC(SYS,OPTIONS) calculates for the LTI descriptor
%  system SYS = (A-lambda*E,B,C,D) with the transfer-function matrix
%  G(lambda), an additive decomposition
%
%          G(lambda) = G1(lambda) + G2(lambda) ,
%
%  such that G1(lambda), the transfer-function matrix of
%  SYS1 = (A1-lambda*E1,B1,C1,D), has only poles in a certain domain
%  of interest Cg of the complex plane and G2(lambda), the
%  transfer-function matrix of SYS2 = (A2-lambda*E2,B2,C2,0), has
%  only poles outside of Cg. The domain of interest Cg is defined via the
%  options structure OPTIONS, as follows:
%  OPTIONS.tol   - specifies the relative tolerance TOL < 1 used for
%                  rank determinations. (Default: internally computed)
%  OPTIONS.smarg - specifies SMARG, the stability margin for the
%                  stable eigenvalues of SYS, such that,
%                  in the continuous-time case, the stable eigenvalues
%                  have real parts less than or equal to SMARG, and
%                  in the discrete-time case, the stable eigenvalues
%                  have moduli less than or equal to SMARG.
%                  (Default: SMARG = -sqrt(eps) for a continuous-time system
%                  and SMARG = 1-sqrt(eps) for a discrete-time system)
%  OPTIONS.job   - option for specific spectral separation tasks
%                 'finite'   - SYS1 has only finite poles and
%                              SYS2 has only infinite poles (default)
%                 'infinite' - SYS1 has only infinite poles and
%                              SYS2 has only finite poles
%                 'stable'   - SYS1 has only stable poles and
%                              SYS2 has only unstable and infinite poles
%                 'unstable' - SYS1 has only unstable and infinite poles
%                              and SYS2 has only stable poles
%
%  [SYS1,SYS2,Q,Z] = GSDEC(SYS,OPTIONS) provides additionally the
%  employed left and right transformation matrices, respectively, used
%  to achieve the block diagonalization of the matrices A and E as
%
%               ( A1   0  )               ( E1   0  )
%      Q*A*Z  = (         ) ,     Q*E*Z = (         )
%               ( 0    A2 )               ( 0    E2 )
%
%  and to obtain the transformed matrices
%
%        Q*B = ( B1 )    and  C*Z = ( C1 C2 ).
%              ( B2 )
%
%  See also SL_GSEP.

%  Author:      A. Varga, 17.01.2016.
%  Revision(s): A. Varga, 07.10.2016, 08.05.2017.
%
%  Method:  The spectral separation based approach of [1] is used.
%
%  References:
%  [1] B. Kagstrom and P. Van Dooren. 
%      Additive decomposition of a transfer function with respect
%      to a specified region. In Proc. MTNS Symp., Brussels, 1989.

narginchk(1,2)
nargoutchk(0,4)

if ~isa(sys,'ss')
   error('The input system SYS must be a state space system')
end

if nargin == 1
   options = struct('tol',0,'job','finite'); 
end

[a,b,c,d,e,Ts] = dssdata(sys);
standsys = isequal(e,eye(size(a,1)));  
if Ts == 0
   discr = 0;
   smax = -sqrt(eps); smin = -inf;
else
   discr = 1;
   smax = 1-sqrt(eps); smin = 0;
end

% decode options and set default values

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0,'<',1},'','OPTIONS.tol') 
else
   tol = 0;
end

% stability margin
if isfield(options,'smarg')
   smarg = options.smarg;
%   validateattributes(smarg, {'double'},{'real','scalar','<',smax,'>=',smin},'','OPTIONS.smarg') 
else
   smarg = smax;
end

if isfield(options,'job')
   jobopt = options.job;
   validateattributes(jobopt,{'char'},{'nonempty'},'','OPTIONS.job') 
else
   jobopt ='finite';
end

% set method and job options for SL_GSEP
switch jobopt
    case 'finite'
       job = 0; meth = 4;
    case 'infinite'
       job = 0; meth = 4;
    case 'stable'
       job = 1; meth = 6;
    case 'unstable'
       job = 1; meth = 6;
    otherwise
       error('No such OPTIONS.job option')
end

% perform spectral separation of eigenvalues
if nargout <= 2
   [a,e,b,c,dims] = sl_gsep(meth,a,e,b,c,job,tol,discr,smarg); 
else
   [a,e,b,c,dims,~,q,z] = sl_gsep(meth,a,e,b,c,job,tol,discr,smarg); 
end

n = size(a,1); n1 = dims(1); 
i1 = 1:n1; i2 = n1+1:n;

if strcmp(jobopt,'infinite') || strcmp(jobopt,'stable')
   if standsys
      sys1 = ss(a(i1,i1),b(i1,:),c(:,i1),d,Ts);
      sys2 = ss(a(i2,i2),b(i2,:),c(:,i2),zeros(size(d)),Ts);
   else
      sys1 = dss(a(i1,i1),b(i1,:),c(:,i1),d,e(i1,i1),Ts);
      sys2 = dss(a(i2,i2),b(i2,:),c(:,i2),zeros(size(d)),e(i2,i2),Ts);
   end
else 
   if standsys
      sys1 = ss(a(i2,i2),b(i2,:),c(:,i2),d,Ts);
      sys2 = ss(a(i1,i1),b(i1,:),c(:,i1),zeros(size(d)),Ts);
   else
      sys1 = dss(a(i2,i2),b(i2,:),c(:,i2),d,e(i2,i2),Ts);
      sys2 = dss(a(i1,i1),b(i1,:),c(:,i1),zeros(size(d)),e(i1,i1),Ts);
   end
end

% end of GSDEC
end


