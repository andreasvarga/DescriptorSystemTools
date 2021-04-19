function [sysr,rankE] = gss2ss(sys,tol,Eshape)
%GSS2SS Conversions to SVD-like coordinate forms without non-dynamic modes. 
%       [SYSR,rankE] = GSS2SS(SYS) calculates for a given descriptor system 
%       SYS = (A-lambda E,B,C,D) an input-output equivalent descriptor 
%       system SYSR = (Ar-lambda Er,Br,Cr,Dr) having the same transfer 
%       function matrix, where Er is in the form Er = diag(E1,0), with E1 
%       a rankE-by-rankE identity matrix. If Er is non-singular, then 
%       the resulting system SYSR is a standard state-space system.
%
%       [SYSR,rankE] = GSS2SS(SYS,TOL) uses the tolerance TOL for the 
%       determination of rank of E. If TOL = 0, an internally computed
%       default value is used. 
%
%       [SYSR,rankE] = GSS2SS(SYS,TOL,ESHAPE) specifies the option for
%       the resulting shape of submatrix E1 as follows:
%       ESHAPE = 'diag'     - E1 diagonal; the diagonal elements are the 
%                             decreasingly ordered nonzero singular values 
%                             of E
%                'triu'     - E1 upper triangular 
%                'ident'    - E1 identity (default) 
%
%       See also SL_GMINR and SL_GSTRA.

%  Author: A. Varga, 27-11-2015.
%  Revision(s): 

if ~isa(sys,'ss')
   error('The input system SYS must be a state space system object')
end

[a,b,c,d,e]=dssdata(sys);

% finisch for a standard state space system
if isequal(e,eye(size(e,1)))
    sysr = sys;
    rankE = size(e,1);
    return
end

% set default inputs
if nargin < 2
   tol = 0;
end

if nargin < 3
   Eshape = 'ident';
end

% set job option for SL_GMINR
if strcmp(Eshape,'triu')
   job = 0;
elseif strcmp(Eshape,'diag')
   job = 0;
elseif strcmp(Eshape,'ident')
   job = 1;
else
    error('Improper shape option')
end

% try to exploit the generalized Hessenberg form of (A,E)
n = size(e,1);
if job == 1 && istriu(e) && isbanded(a,1,n-1) ...
            && rcond(e) > eps && max(max(abs(e))) > tol
   indneg = (diag(e) < 0); 
   if any(indneg)
      e(indneg,:) = -e(indneg,:);
      a(indneg,:) = -a(indneg,:);
      b(indneg,:) = -b(indneg,:); 
   end
   e2 = sqrtm(e); 
   sysr = ss((e2\a)/e2,e2\b,c/e2,d,sys);
   rankE = n;
   return
end


meth = 2;
[a,e,b,c,d,info] = sl_gminr(meth,a,e,b,c,d,tol,job);
if strcmp(Eshape,'diag')
  meth = 4;  
  [a,e,b,c,~,~,~]=sl_gstra(meth,a,e,b,c,tol);
end

rankE = info(8);   
if rankE == size(e,1) && job == 1
   sysr = ss(a,b,c,d,sys);
else
   sysr = dss(a,b,c,d,e,sys);
end

end

% end GSS2SS
