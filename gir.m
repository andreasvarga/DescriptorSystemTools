function sysr = gir(sys,tol,jobopt)
%GIR   Reduced order realizations of LTI descriptor systems.
%   SYSR = GIR(SYS)  calculates for a given descriptor system 
%   SYS = (A-lambda*E,B,C,D) an input-output equivalent irreducible
%   descriptor system SYSR = (Ar-lambda*Er,Br,Cr,D) having the
%   same transfer-function matrix. 
%
%   SYSR = GIR(SYS,TOL) uses TOL as tolerance for rank determinations. 
%   (Default: internally computed if TOL = 0)
%
%   SYSR = GIR(SYS,TOL,JOBOPT) determines various reduced order 
%   realizations, as specified by the option parameter JOBOPT, as follows:
%   JOBOPT = 'irreducible'    - compute an irreducible realization (Default)
%            'finite'         - compute a finite controllable and
%                               finite observable realization
%            'infinite'       - compute an infinite controllable and
%                               infinite observable realization
%            'contr'          - compute a controllable realization
%            'obs'            - compute an observable realization
%            'finite_contr'   - compute a finite controllable realization
%            'infinite_contr' - compute an infinite controllable realization
%            'finite_obs'     - compute a finite observable realization
%            'infinite_obs'   - compute an infinite observable realization
%
%
%   See also GMINREAL.

%  Author:      A. Varga, 27.12.2015.
%  Revision(s): A. Varga, 03.05.2017, 01.06.2018, 15.09.2018. 
%
% References
% 
% [1] A. Varga
%     Computation of Irreducible Generalized State-Space Realizations.
%     Kybernetika, vol. 26, pp. 89-106, 1990.
% 

narginchk(1,3)
nargoutchk(0,1)

if ~isa(sys,'ss')
   error('The input system SYS must be a state-space system object')
end


% set default inputs
if nargin < 2
   tol = 0;
else
   validateattributes(tol, {'double'}, {'real','scalar', '>=', 0},'','TOL',2) 
end

if nargin < 3
   jobopt = 'irreducible';
else
   validateattributes(jobopt,{'char'},{'nonempty'},'','JOBOPT',3) 
end

% set job and systype options for SL_GMINR
task = 1;
switch jobopt
    case 'irreducible'
       job = 0; systype = 0;
    case 'finite'
       job = 0; systype = 1;
    case 'infinite'
       job = 0; systype = 2;
    case 'contr'
       job = 1; systype = 0;
    case 'finite_contr'
       job = 1; systype = 1;
    case 'infinite_contr'
       job = 1; systype = 2;
    case 'obs'
       job = 2; systype = 0;
    case 'finite_obs'
       job = 2; systype = 1;
    case 'infinite_obs'   
       job = 2; systype = 2;
    otherwise
       error('No such JOBOPT option')
end

sysr = sys;

[~,~,M,N] = size(sys); 
s = warning('query','all');
warnon = ~strcmp(s(1).state,'off');

for i = 1:M
    for j = 1:N
        [a,b,c,d,e] = dssdata(sys(:,:,i,j));
        n = size(a,1);
        standsys = isequal(e,eye(n)); 
        % no computations necessary for standard systems and handling
        % impulsive modes
        if ~standsys || systype ~= 2
           [aa,ee,bb,cc,dd,info] = sl_gminr(task,a,e,b,c,d,tol,job,systype);
           info(info(1:4) < 0) = 0; 
           nc = sum(info(1:2)); no = sum(info(3:4)); 
           if warnon
              if N > 1 || M > 1
                 mesi = ['System(',num2str(i),',',num2str(j),'): ']; 
                 if nc, disp([mesi,num2str(nc) ' uncontrollable eigenvalue(s) removed']), end
                 if no, disp([mesi,num2str(no) ' unobservable eigenvalue(s) removed']), end
              else
                 if nc, disp([num2str(nc) ' uncontrollable eigenvalue(s) removed']), end
                 if no, disp([num2str(no) ' unobservable eigenvalue(s) removed']), end
              end
           end
           % save again only if order reduction took place
           if size(aa,1) < n
              if standsys
                 sysr(:,:,i,j) = ss(aa,bb,cc,dd,sys);
              else
                 sysr(:,:,i,j) = dss(aa,bb,cc,dd,ee,sys);
              end
           end
        end
    end
end

% end GIR
end
