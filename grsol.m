function [sysx,info,sysgen] = grsol(sysg,sysf,options)
%GRSOL  Solution of the linear rational matrix equation G*X = F.
%   [SYSX,INFO] = GRSOL(SYSG,SYSF,OPTIONS) returns a state-space  
%   realization SYSX of the solution X(lambda) of the linear rational 
%   equation
%
%      G(lambda)*X(lambda)=F(lambda)                       (1)
%
%   given the state-space realizations SYSG and SYSF of G(lambda) and 
%   F(lambda). The OPTIONS structure (optional) specifies various user 
%   options, as follows:
%   OPTIONS.tol     - specifies the relative tolerance for rank 
%                     determinations. (Default: internally computed)
%   OPTIONS.sdeg    - specifies a prescribed stability degree for the free
%                     eigenvalues of SYSX 
%                     (Default: [], i.e., no stabilization performed)
%   OPTIONS.poles   - specifies a complex conjugated set of desired poles
%                     to be assigned as free eigenvalues of SYSX 
%                     (Default: []) 
%   OPTIONS.mindeg  - minimum degree solution 
%                     true  - determine minimum order solution
%                     false - determine a particular solution which has 
%                             possibly non-minimal order (default) 
%   The resulting INFO structure contains additional information:
%   INFO.nrank is the normal rank of G(lambda); 
%   INFO.rdeg is a vector which contains the relative column degrees of 
%       X(lambda) (i.e., the numbers of integrators/delays needed to make 
%       each column of X(lambda) proper); 
%   INFO.tcond is the maximum of the Frobenius-norm condition numbers of 
%        the employed non-orthogonal transformation matrices.
%   INFO.fnorm is the Frobenius-norm norm of the employed 
%       state-feedback/feedforward used for dynamic cover computation if 
%       OPTIONS.mindeg = true, or for stabilization of the free eigenvalues
%       if OPTION.sdeg is non-empty;  
%   INFO.nr is the number of freely assignable poles of the solution X(lambda). 
%   Note: Large values of INFO.tcond and/or INFO.fnorm indicate possible
%         loss of numerical stability. 
% 
%   [SYSX,INFO] = GRSOL(SYSGF,MF,OPTIONS) uses the compound realization 
%   SYSGF of [G(lambda) F(lambda)], where F(lambda) has MF columns.
%
%   [SYSX,INFO,SYSGEN] = GRSOL(SYSG,SYSF,OPTIONS) determines additionally
%   a generator of all solution SYSGEN = [ SYS0 SYSN ], where SYS0 is a 
%   particular solution of (1) and SYSN is a proper right nullspace basis
%   of G(lambda) satisfying G(lambda)*SYSN(lambda) = 0. All solutions 
%   of (1) can be generated as
%               SYSX = SYS0 + SYSN*SYSY, 
%   where SYSY is an arbitrary system with appropriate dimensions. 
%   The descriptor system realization of SYSGEN is usually not minimal 
%   (uncontrollable and/or non-dynamic modes present) and has the form 
%   (Ag-lambda*Eg,[B0 BN],Cg,[D0 DN]), where
%
%                     ( Ar-lambda*Er     *              *      )  
%      Ag-lambda*Eg = (     0        Af-lambda*Ef       *      ), 
%                     (     0            0        Ai-lambda*Ei ) 
%
%                   ( B1 | Br )
%      [B0 | BN ] = ( B2 | 0  ),  Cg  =   ( Cr   *    *  ) 
%                   ( B3 | 0  )
%   
%   with Er, Ef and Ai invertible and upper triangular, Ei nillpotent
%   and upper triangular, and DN full row rank. A minimal order descriptor 
%   system realization of the proper basis SYSN is (Ar-lambda*Er,Br,Cr,DN).
%   Br and DN have M-INFO.nrank columns, where M is the number of 
%   inputs of SYSG.
%
%   The INFO structure contains further information on the dimensions of
%   the square diagonal blocks of the pencil Ag-lambda*Eg, as follows:
%   INFO.nr is the dimension of Ar-lambda*Er (also the row dimension of Br);
%   INFO.nf is the dimension of Af-lambda*Ef (the number of finite eigenvalues);
%   INFO.ninf is the dimension of Ai-lambda*Ei (the number of infinite
%      eigenvalues).
%
%   See also GLSOL.

%  Copyright 2016-2018 A. Varga 
%  Author:       A. Varga, 14.01.2016.
%  Revision(s):  A. Varga, 04.10.2016, 28.10.2016, 16.05.2017, 21.08.2018,
%                          26.07.2019.
%
%  Method:  The method of [1] to solve rational systems is used.
%
%  References:
%  [1] A. Varga, "Computation of least order solutions of linear 
%      rational equations", Proc. MTNS'04, Leuven, Belgium, 2004.

narginchk(2,3)
nargoutchk(0,3)

if ~isa(sysg,'ss')  
   error('The input system SYSG must be an SS object')
end  

[p,m] = size(sysg);
if ~isa(sysf,'ss')
   if ~isa(sysf,'double') 
      error('SYSF must be an SS object or a positive integer')
   else
      validateattributes(sysf,{'double'},{'integer','scalar','>=',0,'<=',m},'','MF')
   end
end  

discr = sysg.Ts ~= 0; 
if discr
    smax = 1-sqrt(eps); smin = 0;
else
    smax = -sqrt(eps);  smin = -inf;
end

% decode options
if nargin < 3
   options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

% desired stability degree
if isfield(options,'sdeg')
   sdeg = options.sdeg;
   if ~isempty(sdeg)
      validateattributes(sdeg, {'double'},{'real','scalar','<=',smax,'>=',smin},'','OPTIONS.sdeg') 
   end
else
   sdeg = [];
end

% desired poles
if isfield(options,'poles')
   poles = options.poles;
   if ~isempty(poles)
      validateattributes(poles, {'double'},{'vector'},'','OPTIONS.poles')
   end
else
   poles = [];
end
if ~isempty(poles)
    t = poles;
    if (min(size(t)) ~= 1 || ~isequal(sort(t(imag(t)>0)),sort(conj(t(imag(t)<0)))) )
        error('OPTIONS.poles must be a self-conjugated complex vector')
    end 
end    

if isempty(sdeg) && isempty(poles)
   stabilize = false;
else
   stabilize = true;
end

% minimum degree solution option
if isfield(options,'mindeg')
   mindeg = options.mindeg;
   validateattributes(mindeg, {'logical'},{'binary'},'','OPTIONS.mindeg') 
else
   mindeg = false;
end

Ts = sysg.Ts;
if isa(sysf,'double') 
    mf = sysf; 
    [a,b,c,d,e] = dssdata(sysg);
    m = m-mf;
else
   if ~isa(sysf,'ss')
      error('The input system SYSF must be an SS object')
   end
   if (discr && (sysf.Ts <= 0)) || (~discr && (sysf.Ts > 0))
      error('The systems SYSG and SYSF must have the same type')
   end
   if discr && (Ts ~= sysf.Ts) 
      error('The systems SYSG and SYSF must have the same sampling period')
   end
   if p ~= size(sysf,1) 
      error('The systems SYSG and SYSF must have the same number of outputs')
   end
   mf = size(sysf,2);
   [a,b,c,d,e] = dssdata(gir([sysg,sysf],tol));
end    
    
n = size(a,1);
bf = b(:,m+1:end); df = d(:,m+1:end);
b = b(:,1:m); d = d(:,1:m);
tola = tol; 
if tol == 0
   % set absolute tolerance for Kronecker-form computation
   tola = n*eps*max([norm(e,1) norm(a,1) norm(b,1),norm(c,inf)]);
end


% compute the Kronecker-like form of the system matrix of G to obtain
%                 ( Br  Ar - s Er      *             *           *     )
%   Q'*(A-sE)*Z = (  0     0      Afin - s Efin      *           *     )
%                 (  0     0           0        Ainf - s Einf    *     )
%                 (  0     0           0            0         Al - s El)
 
A = [a b; c d]; E = [e zeros(n,m);zeros(p,n+m)]; 
[A,E,info1,Q,Z] = gklf(A,E,tola,'reverse'); 

%                        ( B1 ) 
% compute Q'*[-Bf;-Df] = ( B2 )
%                        ( B3 ) 
%                        ( B4 )
f = Q'*[-bf;-df]; 

% determine the orders of Ar, Af and Ai
ninf = sum(info1.minf); nf = info1.mf; nr = sum(info1.mr); 
nl = sum(info1.nl);
nreg = nf+ninf; ninv = nr+nreg+nl; mr = n+m-ninv; 
info = struct('nrank',ninv-n,'rdeg',[],'tcond',1,'fnorm',0,'nr',nr, ...
                 'nf',nf,'ninf',ninf);

% check compatibility condition, i.e., B4 = 0.
if sum(info1.ml)
   if norm(f(nr+nreg+1:end,:),1) >= tola
       error('System not compatible')
   end
end

% form X0 = (A0-sE0,B0,C0,D0) and XN = (Ar-sEr,Br,Cr,DN), where
%
%                ( Ar - s*Er     *           *        )        ( B1 )
%       A0-sE0 = (    0       Af - s*Ef      *        ) , B0 = ( B2 )
%                (    0          0        Ai - s*Einf )        ( B3 )
%      
%           C0 = ( Cr  Cf Ci );  D0 = 0, with
%
%             [ 0  Im ]*Z := [ DN Cr  Cf Ci Cl ];
n0 = nr+nreg;
A0 = A(1:n0,mr+1:mr+n0); E0 = E(1:n0,mr+1:mr+n0); 
B0 = f(1:n0,:);    C0 = Z(n+1:n+m,mr+1:mr+n0); D0 = zeros(m,mf); 
Br = A(1:nr,1:mr); BN = [Br;zeros(nreg,mr)]; DN = Z(n+1:n+m,1:mr); % CN = C0


i2 = nr+nf+1:n0;
rdeg = zeros(1,mf);
%if length(minf) > 1
if length(info1.minf) > 1
   % compute relative column degrees as the number of controllable 
   % infinite poles (in excess with 1 with respect to infinite eigenvalues)
   for i=1:mf
       [~,~,~,~,dims] = sl_gstra(8,A0(i2,i2),E0(i2,i2),B0(i2,i),C0(:,i2),2,tola);
       ni = dims(1)-1;
       if ni > 0
           rdeg(i)=ni;
       end
   end
end  
info.rdeg = rdeg;
if nargout > 2 && mindeg
   % form a pair of generators 
   sysgen = dss(A0,[B0 BN],C0,[D0 DN],E0,Ts);
end
    
  
if mindeg 
   % block-diagonalize E0 to allow working on the proper part
   i1 = 1:nr; i2 = nr+1:n0;
   Y = -E0(i1,i1)\E0(i1,i2);
   A0(i1,i2) = A0(i1,i2)+A0(i1,i1)*Y;
   C0(:,i2) = C0(:,i2)+C0(:,i1)*Y;
   
   % form a pair of generators for the proper part 
   sysr = dss(A0(i1,i1),[ B0(i1,:) A0(i1,i2) Br],C0(:,i1),[D0 C0(:,i2) DN],E0(i1,i1),Ts);
   
   % compute minimum order for the proper part
   [sysr,info2] = grmcover2(sysr,mf+nreg,tol); 
   [aa,bb,cc,dd,ee] = dssdata(sysr); na = size(aa,1);
   i3 = mf+1:mf+nreg;
   A0 = [aa bb(:,i3); zeros(nreg,na) A0(i2,i2)]; E0 = blkdiag(ee,E0(i2,i2));
   B0 = [ bb(:,1:mf); B0(i2,:)]; C0 = [cc dd(:,i3)]; D0 = dd(:,1:mf); 
   info.tcond = info2.tcond;
   info.fnorm = max(info2.fnorm,info2.gnorm);
else
   F = [];
   if stabilize  
      % make spurious poles stable
      if nr
         i1 = 1:nr; 
         opt_poles = struct('tol',tol);
         F = gsfstab(A0(i1,i1),E0(i1,i1),Br,poles,sdeg,opt_poles);
         A0(i1,i1) = A0(i1,i1)+Br*F; C0(:,i1) = C0(:,i1)+DN*F;
      end
   end
   info.fnorm = norm(F,'fro');
   if nargout > 2
      % form a pair of generators 
      sysgen = dss(A0,[B0 BN],C0,[D0 DN],E0,Ts);
   end
end

% eliminate possible uncontrollable eigenvalues
task = 1; job = 1; systype = 0;
[A0,E0,B0,C0,D0]=sl_gminr(task,A0,E0,B0,C0,D0,tol,job,systype);
% eliminate possible simple infinite eigenvalues
% task = 2; job = 0;
% [A0,E0,B0,C0,D0] = sl_gminr(task,A0,E0,B0,C0,D0,tola,job);
sysx = dss(A0,B0,C0,D0,E0,Ts);
sysx = gss2ss(sysx,tol,'triu');

% end GRSOL
end
