function [sysr,rankE] = gss2ss(sys,tol,Eshape)
%GSS2SS Conversions to SVD-like forms without non-dynamic modes.
%       [SYSR,rankE] = GSS2SS(SYS) calculates for a given descriptor system 
%       SYS = (A-lambda*E,B,C,D) an input-output equivalent descriptor 
%       system SYSR = (Ar-lambda*Er,Br,Cr,Dr) having the same 
%       transfer-function matrix, where Er is in the form Er = diag(E1,0), 
%       with E1 a rankE-by-rankE identity matrix. If Er is non-singular, 
%       then the resulting system SYSR is a standard state-space system.
%
%       [SYSR,rankE] = GSS2SS(SYS,TOL) uses the tolerance TOL for rank 
%       determinations. If TOL = 0, internally computed default values
%       are used. 
%
%       [SYSR,rankE] = GSS2SS(SYS,TOL,ESHAPE) specifies the option for
%       the resulting shape of the submatrix Er1 as follows:
%       ESHAPE = 'diag'     - E1 diagonal; the diagonal elements are the 
%                             decreasingly ordered nonzero singular values 
%                             of E
%                'triu'     - E1 upper triangular 
%                'ident'    - E1 identity (default) 
%

%  Author:      A. Varga, 10.11.2016.
%  Revision(s): A. Varga, 03.05.2017, 29.08.2017.

narginchk(1,3)
nargoutchk(0,2)

if ~isa(sys,'ss')
   error('The input system SYS must be a state space system object')
end

% set default inputs
if nargin < 2
   tol = 0;
else
   validateattributes(tol, {'double'}, {'real','scalar', '>=', 0},'','TOL') 
end

if nargin < 3
   Eshape = 'ident';
else
   validateattributes(Eshape,{'char'},{'nonempty'},'','ESHAPE') 
end

[a,b,c,d,e]=dssdata(sys);
n = size(a,1);

% finish for a standard state space system
if isequal(e,eye(n))
    sysr = sys;
    rankE = n;
    return
end


if tol == 0 
   % set absolute tolerance for determination of rank of E 
   tole = n*n*eps(norm(e,1));
   tola = n*n*eps(norm(a,1));
%   tole = max(tola,tole);
else
    tole = tol; 
    tola = tol; 
end


% try to exploit the generalized Hessenberg form of (A,E)
maxE = max(max(abs(e)));
if strcmp(Eshape,'ident') && istriu(e) && isbanded(a,1,n-1) ...
            && rcond(e) > eps && maxE > tole
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


Ediag = strcmp(Eshape,'diag') || strcmp(Eshape,'ident');

% Using the non-orthogonal transformation matrices Q and Z, reduce the 
% matrices A, E, B and C to the forms
%
%              [ At11  At12 At13 ]                 [ Et11  0  0 ]   
% At = Q*A*Z = [ At21   I    0   ] ,  Et = Q*E*Z = [  0    0  0 ] , 
%              [ At31   0    0   ]                 [  0    0  0 ]
%
%            [ Bt1 ] 
% Bt = Q*B = [ Bt2 ] ,  Ct = C*Z = [ Ct1  Ct2  Ct3 ]
%            [ Bt3 ]
%
if Ediag
   % reduce (A,E) to a SVD-coordinate form using orthogonal transformations
   if maxE <= tole
      % skip reduction of E if E is negligible
      rankE = 0;
      [et,at,bt,ct,ra22] = sl_gstra(4,e,a,b,c,tola);
   else
      [at,et,bt,ct,ranks] = sl_gstra(5,a,e,b,c,tol);
      rankE = ranks(1); ra22 = ranks(2); 
   end
   ninf = n-rankE-ra22;
   i1 = 1:rankE; i2 = rankE+1:rankE+ra22; i3 = n-ninf+1:n; 
   % make A22 = diag(I,0)
   if ra22
      tid = 1./sqrt(diag(at(i2,i2))); at(i2,i2) = eye(ra22);
      at(i2,i1) = bsxfun(@times,at(i2,i1),tid); 
      at(i1,i2) = bsxfun(@times,at(i1,i2),tid'); 
      bt(i2,:)  = bsxfun(@times,bt(i2,:),tid); 
      ct(:,i2)  = bsxfun(@times,ct(:,i2),tid'); 
   end
elseif strcmp(Eshape,'triu')
   % reduce (A,E) to a SVD-like coordinate form with E upper-triangular
   if maxE <= tole
      % skip reduction of E if E is negligible
      rankE = 0; i1 = [];
      [et,at,bt,ct,ra22] = sl_gstra(4,e,a,b,c,tola);
   else
      [at,et,bt,ct,rankE] = sl_gstra(6,a,e,b,c,tole);
      i1 = 1:rankE; i2 = rankE+1:n; 
      if max(max(abs(at(i2,i2)))) <= tola
          ra22 = 0;
      else
         [et(i2,i2),at(i2,i2),bt(i2,:),ct(:,i2),ra22,Q2,Z2] = ...
             sl_gstra(4,et(i2,i2),at(i2,i2),bt(i2,:),ct(:,i2),tola);
         at(i2,i1) = Q2'*at(i2,i1); at(i1,i2) = at(i1,i2)*Z2;
      end
   end 
   ninf = n-rankE-ra22;
   i2 = rankE+1:rankE+ra22; i3 = n-ninf+1:n; 
   % make A22 = diag(I,0)
   if ra22
      tid = 1./sqrt(diag(at(i2,i2))); at(i2,i2) = eye(ra22);
      at(i2,i1) = bsxfun(@times,at(i2,i1),tid); 
      at(i1,i2) = bsxfun(@times,at(i1,i2),tid'); 
      bt(i2,:)  = bsxfun(@times,bt(i2,:),tid); 
      ct(:,i2)  = bsxfun(@times,ct(:,i2),tid'); 
   end
else
    error('Improper shape option for the E matrix')
end
% 
if rankE < n
   % apply residualization formulas to eliminate non-dynamic modes
   % Dr = D-Ct2*Bt2
   d = d - ct(:,i2)*bt(i2,:);
   % Br = [ B1-At12*Bt2 ]
   %      [    Bt3      ]
   bt = [ bt(i1,:)-at(i1,i2)*bt(i2,:); bt(i3,:)];   
   % Cr = [ C1-Ct2*At21  Ct3 ]
   ct = [ ct(:,i1)-ct(:,i2)*at(i2,i1) ct(:,i3) ];
   % Ar = [ A11-At12*At21  At13 ]
   %      [    At31          0  ]
   at = [ at(i1,i1)-at(i1,i2)*at(i2,i1) at(i1,i3);
          at(i3,i1) zeros(ninf) ];
   % Er = [ E11  0 ]
   %      [  0   0 ]
   et = blkdiag(et(i1,i1),zeros(ninf));    
end
   
if ~isempty(et) && strcmp(Eshape,'ident') 
   % make E11 identity  
   tid = 1./sqrt(diag(et(i1,i1))); 
   at(i1,:) = bsxfun(@times,at(i1,:),tid); 
   at(:,i1) = bsxfun(@times,at(:,i1),tid'); 
   bt(i1,:) = bsxfun(@times,bt(i1,:),tid); 
   ct(:,i1) = bsxfun(@times,ct(:,i1),tid'); 
   if size(at,1) == rankE
       et = [];
   else
       et(i1,i1) = eye(rankE);
   end
end

if isempty(et)
   sysr = ss(at,bt,ct,d,sys);
else
   sysr = dss(at,bt,ct,d,et,sys);
end

% end GSS2SS
end


