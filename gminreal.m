function [sysmin,info] = gminreal(sys,tol,ndmonlyflag)
%GMINREAL   Minimal realization of a LTI descriptor system.
%  [SYSMIN,INFO] = GMINREAL(SYS) calculates for a given descriptor system 
%  SYS = (A-lambda*E,B,C,D), with transfer-function matrix G(lambda),
%  an input-output equivalent minimal descriptor system 
%  SYSMIN = (Am-lambda*Em,Bm,Cm,Dm), such that 
%                              -1                       -1                
%     G(lambda) =  C(lambda*E-A) B + D = Cm(lambda*Em-Am) Bm + D.
%  INFO is a three-dimensional vector which contains information on the 
%  number of removed eigenvalues, as follows:
%  INFO(1) - the number of removed uncontrollable eigenvalues;
%  INFO(2) - the number of removed unobservable eigenvalues;
%  INFO(3) - the number of removed non-dynamic (infinite) eigenvalues.
%
%  [SYSMIN,INFO] = GMINREAL(SYS,TOL) uses tolerance TOL for rank 
%  determinations.
%
%  [SYSMIN,INFO] = GMINREAL(SYS,TOL,'ndmonly') removes only the 
%  non-dynamic modes (eigenvalues).  
%         
%  See also GIR.

%  Author:      A. Varga, 02.11.2016.
%  Revision(s): A. Varga, 29.04.2017.  
%
%  References:
%  [1] A. Varga, "Computation of irreducible generalized state-space 
%      realizations", Kybernetika, vol. 26, pp. 89-106, 1989.

narginchk(1,3)
nargoutchk(0,2)

if nargin < 2
   tol = 0;
else
   validateattributes(tol, {'double'}, {'real','scalar','>=',0,'<',1},'','TOL') 
end

if nargin < 3
   ndmonly = false;
else
   validateattributes(ndmonlyflag,{'char'},{'nonempty'},'','NDMONLYFLAG') 
   ndmonly = strcmp(ndmonlyflag,'ndmonly');
   if ~ndmonly
      error('Improper option: use instead ''ndmonly'' to remove non-dynamic modes only')
   end
end
 
% set task and job options for SL_GMINR
if ndmonly
    task = 2; job = 1;
else
    task = 3; job = 0;
end

[a,b,c,d,e]=dssdata(sys);
[aa,ee,bb,cc,dd,info1] = sl_gminr(task,a,e,b,c,d,tol,job);
s = warning('query','all');
info(1) = max([info1(1),0])+max([info1(2),0]);
info(2) = max([info1(3),0])+max([info1(4),0]);
info(3) = max([info1(5),0]);
if ~strcmp(s(1).state,'off') && nargout < 2
   if info(1), disp([num2str(info(1)) ' uncontrollable eigenvalue(s) removed']), end
   if info(2), disp([num2str(info(2)) ' unobservable eigenvalue(s) removed']), end
   if info(3), disp([num2str(info(3)) ' non-dynamic eigenvalue(s) removed']), end
end

if isequal(ee,eye(size(aa,1)))
   sysmin = ss(aa,bb,cc,dd,sys);
else
   sysmin = dss(aa,bb,cc,dd,ee,sys);
end

% end GMINREAL
end


