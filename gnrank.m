function nr = gnrank(sys,tol,freq)
%GNRANK  Normal rank of the transfer function matrix of a LTI system.
%  NR = GNRANK(SYS) returns NR, the normal rank of the transfer-function
%  matrix of the LTI system SYS. 
%
%  NR = GNRANK(SYS,TOL) uses tolerance TOL for rank determinations.
% 
%  NR = GNRANK(SYS,TOL,FREQ) uses the complex frequency values in the 
%  vector FREQ to estimate the normal rank NRS of the system matrix
%  pencil S(lambda) = [lambda*E-A B; C D] as the maximum rank of S(lambda)
%  for lambda taking all frequency values contained in FREQ. 
%  (A-lambda*E,B,C,D) is a descriptor system realization of order N of SYS. 
%  The normal rank of the transfer function matrix of SYS is NR = NRS - N.
% 
%  See also GZERO.

%  Author:       A. Varga, 15.12.2016.
%  Revision(s):  A. Varga, 25.04.2017, 19.02.2018, 14.09.2018, 29.12.2018. 
%
%  References:
%  [1]  Misra P., Van Dooren, P., Varga, A.:
%       Computation of structural invariants of generalized state space systems. 
%       Automatica, vol. 30, pp. 1921-1936, 1994.

narginchk(1,3)
nargoutchk(0,1)

if ~isa(sys,'ss') && ~isa(sys,'tf') && ~isa(sys,'zpk')
   error('The input system SYS must be an LTI system object')
end

if nargin == 1
   tol = 0;
else
   if isempty(tol)
       tol = 0;
   else
       validateattributes(tol, {'double'}, {'real','scalar', '>=', 0},'','TOL',2)
   end
end

if nargin < 3
   freq = [];
else
   validateattributes(freq, {'double'}, {'vector'},'','FREQ',3) 
end

[a,b,c,d,e] = dssdata(sys);

if isempty(a)
   if tol
      nr = rank(d,tol);
   else
      nr = rank(d);
   end
   return
end

n = size(a,1);
nr = 0; 
if isempty(freq)
   [~,ni] = sl_gzero(a,e,b,c,d,tol); 
   nr = max(nr,ni(2)-n);
else
   for i = 1:length(freq)
       if tol
          nr = max(nr,rank([a-freq(i)*e b; c d],tol)-n);
       else
          nr = max(nr,rank([a-freq(i)*e b; c d])-n);
       end
   end
end

% end GNRANK

