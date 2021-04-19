function nr = nrank(sys,tol)
%NRANK  Normal rank of the transfer function matrix of a LTI system.
%  NR = NRANK(SYS) returns NR, the normal rank of the transfer-function
%  matrix of the LTI system SYS.
%
%  NR = NRANK(SYS,TOL) uses tolerance TOL for rank determinations.
% 
%  Note: This version is obsolete and will be replaced in future releases 
%        by GNRANK. 

%  Author:       A. Varga, 15.12.2016.
%  Revision(s):  A. Varga, 25.04.2017, 19.02.2018, 10.09.2018. 
%
%  References:
%  [1]  Misra P., Van Dooren, P., Varga, A.:
%       Computation of structural invariants of generalized state space systems. 
%       Automatica, vol. 30, pp. 1921-1936, 1994.

narginchk(1,2)
nargoutchk(0,1)

if ~isa(sys,'lti')
   error('The input system SYS must be an LTI system object')
end

if nargin == 1
   tol = 0;
else
   validateattributes(tol, {'double'}, {'real','scalar', '>=', 0},'','TOL') 
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

[~,ni] = sl_gzero(a,e,b,c,d,tol); 
nr = ni(2)-size(a,1);

% end NRANK

