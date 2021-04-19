function [syst,sysi1] = gbilin(sys,sys1,options)
%GBILIN  Generalized bilinear transformation. 
%    [SYST,SYSI1] = GBILIN(SYS,SYS1,OPTIONS) computes for a given 
%    descriptor system SYS, with the transfer function matrix G(lambda), 
%    and a first degree real rational transfer function SYS1, with the 
%    transfer function g(delta), the descriptor system realization of  
%    G(g(delta)) corresponding to the bilinear transformation 
%                 lambda = g(delta) = (a*delta+b)/(c*delta+d). 
%    For a continuous-time system SYS1, delta = s, the complex variable in 
%    the Laplace transform, while for a discrete-time system SYS1,  
%    delta = z, the complex variable in the Z-transform. 
%    SYST inherits the sampling-time of SYS1. 
%    SYSI1 is the transfer function  
%                       ginv(lambda) = (d*lambda-b)/(-c*lambda+a)
%    representing the inverse of the bilinear transformation g(delta) 
%    (i.e., g(ginv(lambda)) = lambda).
%    The OPTIONS structure (optional) specifies various user 
%    options, as follows:
%    OPTIONS.tol     - specifies the relative tolerance for rank 
%                      determinations. (Default: internally computed)
%    OPTIONS.compact - option to compute a compact descriptor realization
%                      without non-dynamic modes
%                      true  - determine, a compact descriptor realization 
%                              of SYST without nondynamic modes (default) 
%                      false - disable ellimination of non-dynamic modes  
%    OPTIONS.minimal - option to compute minimal descriptor realization 
%                      true  - determine, a  minimal descriptor
%                              realization of SYST (i.e., irreducible and 
%                              without nondynamic modes) 
%                      false - no minimal realization computed (default) 
%    OPTIONS.ss      - option to compute a standard state-space (if possible)
%                      realizations of SYST:
%                      true  - standard state-space realization, 
%                      false - descriptor system realization (default).  
%
%    See also  GBILIN1.

%  Copyright 2018 A. Varga 
%  Author:    A. Varga, 22.08.2018.
%  Revisions: 

narginchk(2,3)
nargoutchk(0,2)

if ~isa(sys,'ss')
    error('The input system SYS must be an SS object')
end

if ~isa(sys1,'tf') && any(size(sys1) ~= 1)
    error('The input system SYS1 must be a SISO TF object')
end

if nargin < 3
    options = struct('tol',0);
end

% decode options

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

% compact realization option, without non-dynamic modes
if isfield(options,'compact')
   compact = options.compact;
   validateattributes(compact, {'logical'},{'binary'},'','OPTIONS.compact')
else
   compact = true;
end
  
% minimal realization option
if isfield(options,'minimal')
   minimal = options.minimal;
   validateattributes(minimal, {'logical'},{'binary'},'','OPTIONS.minimal')
else
   minimal = false;
end

% standard state-space realization option
if isfield(options,'ss')
   standard = options.ss;
   validateattributes(standard, {'logical'},{'binary'},'','OPTIONS.ss')
else
   standard = false;
end

if compact
   if standard
       opt_ss = 'ident';
   else
       opt_ss = 'triu';
   end
end

sys1 = minreal(sys1); 
num = sys1.num{1}; degn = length(num)-1; 
den = sys1.den{1}; degd = length(den)-1;

validateattributes(num, {'double'}, {'real','vector'},'','numerator',2)
validateattributes(den, {'double'}, {'real','vector'},'','denominator',2)

if degn > 1 || degd > 1 || max(degn,degd) == 0
    error('The McMillan degree of SYS1 must be one')
end

if max(abs(num)) == 0
    error('SYS1 must be nonzero')
end

Ts = sys.Ts;
Ts1 = sys1.Ts;   
if Ts && Ts1 && Ts ~= Ts1
    error('SYS and SYS1 must have the same sampling periods')
end

% assume g(delta) = (a*delta+b)/(c*delta+d)
if degn
   a = num(1); b = num(2);
else
   a = 0; b = num(1);
end
if degd
   c = den(1); d = den(2);
else
   a = a/den(2); b = b/den(2); c = 0; d = 1;
end

[A,B,C,D,E] = dssdata(sys);
[n,m] = size(B); p = size(C,1);

if degd 
   % rational case
   At = [-b*E+d*A d*B; zeros(m,n) -eye(m)];
   Et = [a*E-c*A -c*B; zeros(m,n+m)];
   Bt = [zeros(n,m); eye(m)];
   Ct = [C D]; Dt = zeros(p,m);
   syst = dss(At,Bt,Ct,Dt,Et,Ts1);
   if minimal 
      if standard
         syst = gss2ss(gir(syst,tol),tol);
      else
         syst = gminreal(syst,tol);
      end
   else
      syst = gss2ss(syst,tol,opt_ss);
   end
else
   % polynomial case
   if isequal(E,eye(size(E,1)))
      % preserve standard system form
      syst = ss((A-b*E)/a,B/a,C,D,Ts1);
   else
      syst = dss(A-b*E,B,C,D,a*E,Ts1);
      if standard
          syst = gss2ss(syst,tol);
      end
   end
end

if nargout > 1
   if Ts == 0 && Ts1 == 0
      Tsi = 0;
   elseif Ts == 0 && Ts1 ~= 0
      Tsi = Ts1;
   elseif Ts ~= 0 && Ts1 == 0
      Tsi = 0;
   else %Ts ~= 0 && Ts1 ~= 0
      if Ts1 < 0
          Tsi = Ts;
      else
          Tsi = Ts1;
      end
   end
   sysi1 = tf([d -b],[-c a],Tsi);
end

% end GBILIN
