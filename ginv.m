function [sysinv,info] = ginv(sys,options)
%GINV  Generalized inverses of LTI systems.
%   [SYSINV,NR] = GINV(SYS,OPTIONS) computes a generalized inverse system 
%   SYSINV such that the transfer function matrices G(lambda) and Gi(lambda)
%   of SYS and SYSINV satisfy two or more of the Moore-Penrose conditions 
%         (1) G(lambda)*Gi(lambda)*G(lambda) = G(lambda);        
%         (2) Gi(lambda)*G(lambda)*Gi(lambda) = Gi(lambda);
%         (3) G(lambda)*Gi(lambda) = (G(lambda)*Gi(lambda))'
%         (4) Gi(lambda)*G(lambda) = (Gi(lambda)*G(lambda))'.
%   The OPTIONS structure (optional) specifies various user 
%   options, as follows:
%   OPTIONS.tol     - specifies the relative tolerance for rank 
%                     determinations. (Default: internally computed)
%   OPTIONS.sdeg    - specifies a prescribed stability degree for the free
%                     eigenvalues of SYSINV 
%                     (Default: [], i.e., no stabilization performed)
%   OPTIONS.poles   - specifies a complex conjugated set of desired poles
%                     to be assigned as free eigenvalues of SYSINV 
%                     (Default: []) 
%   OPTIONS.mindeg  - minimum degree solution 
%                     true  - determine a minimum order generalized inverse
%                     false - determine a particular generalized inverse  
%                             which has possibly non-minimal order (default) 
%   OPTIONS.type    - specifies the type of the desired generalized inverse. 
%                     The allowed values are:
%                     '1-2'     - for a generalized inverse which satisfies 
%                                 conditions (1) and (2) (default)
%                     '1-2-3'   - for a generalized inverse which satisfies 
%                                 conditions (1), (2) and (3)
%                     '1-2-4'   - for a generalized inverse which satisfies 
%                                 conditions (1), (2) and (4)
%                     '1-2-3-4' - for the Moore-Penrose pseudoinverse,
%                                 which satisfies all conditions (1)-(4)
%   OPTIONS.freq    - specifies a vector of complex frequency values  
%                     to evaluate the normal rank of G(lambda)  
%                     (Default: [], i.e., a random real frequency value used)
% 
%   The resulting INFO structure contains additional information:
%   INFO.nrank is the normal rank of G(lambda); 
%   INFO.tcond is the maximum of the Frobenius-norm condition numbers of 
%       the employed transformation matrices.
%   INFO.fnorm is the Frobenius norm of the employed 
%       state-feedback/feedforward used for dynamic cover computation if 
%       OPTIONS.mindeg = true, or for stabilization of free eigenvalues if 
%       OPTION.sdeg is non-empty;  
%   INFO.nfp is the number of freely assignable poles of the generalized
%       inverse Gi(lambda). 
%
%   See also GNRANK.

%  Author: A. Varga, 30-07-2019.
%  Revision(s): 
%
%  References:
%   [1] Varga, A.
%       A note on computing range space bases of rational matrices. 
%       arXiv:1707.0048, https://arxiv.org/abs/1707.00489, 2017.
%   [2] A. Varga.
%       Computing generalized inverse systems using matrix pencil methods.
%       Int. J. of Applied Mathematics and Computer Science, vol. 11, 
%       pp. 1055-1068, 2001.

narginchk(1,2)
nargoutchk(0,2)

if ~isa(sys,'ss')
   error('The input system SYS must be a state-space system object')
end


if nargin < 2
   options = struct('tol',0,'type','1-2');
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

% inverse selection option
if isfield(options,'type')
   type = options.type;
   validateattributes(type, {'char'},{'nonempty'},'','OPTIONS.type') 
else
   type = '1-2';
end

switch type
    case '1-2'
       job = 0; 
    case '1-2-3'
       job = 1;
    case '1-2-4'
       job = 2;
    case '1-2-3-4'
       job = 3; 
    otherwise
       error('No such inverse selection option')
end

% test frequencies for rank determination
if isfield(options,'freq')
   freq = options.freq;
   validateattributes(freq, {'double'},{'vector'},'','OPTIONS.freq') 
else
   freq = rand;
end

% desired stability degree
if sys.Ts ~= 0
    smax = 1-sqrt(eps); smin = 0;
else
    smax = -sqrt(eps);  smin = -inf;
end

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

[p,m]=size(sys);
nrank = gnrank(sys,tol,freq);
nfp = 0;
tcond = 1;
fnorm = 0;

if p == m && m == nrank
   % compute a standard inverse for an invertible system 
   sysinv = inv(sys);
   info.nrank = nrank;
   info.nfp = nfp;
   info.tcond = 1;
   info.fnorm = 0;
   return
end

if job == 0 
   % compute an (1,2)-inverse
   if min(p,m) == nrank
      if m == nrank
         % compute a left inverse
         [sysinv,info1] = glsol([sys; eye(m)],m,options);
         if info1.nrank ~= nrank
            warning('Inconsistent rank evaluations: check tolerances')
         end
         if ~mindeg
            nfp = info1.nl;
         end
      else
         % compute a right inverse
         [sysinv,info1] = grsol([sys eye(p)],p,options);
         if info1.nrank ~= nrank
            warning('Inconsistent rank evaluations: check tolerances')
         end
         if ~mindeg
            nfp = info1.nr;
         end
      end
      tcond = info1.tcond;
      fnorm = info1.fnorm;
   else
      if mindeg
         % compute a full rank factorization G = U*V
         [U,V] = grange(sys,struct('tol',tol));
         if nrank ~= size(U,2)
            warning('Inconsistent rank evaluations: check tolerances')
         end
         % compute a left inverse of least order of U
         [UL,info1] = glsol([U; eye(nrank)],nrank,options);
         tcond = info1.tcond;
         fnorm = info1.fnorm;
         % compute a right inverse of least order of V
         [VR,info1] = grsol([V eye(nrank)],nrank,options);
         tcond = max(tcond,info1.tcond);
         fnorm = max(fnorm,info1.fnorm);
         sysinv = gminreal(VR*UL,tol);
      else
         % compute a generalized 1-2 inverse using the Kronecker form 
         % of the system pencil for a system with no full rank TFM
         [a,b,c,d,e,Ts] = dssdata(sys);
         n = size(a,1);
         tola = tol; 
         if tol == 0
            % set absolute tolerance for Kronecker-form computation
            tola = eps((n+p)*(n+m))*max([norm(a,1) norm(e,1) norm(b,1), norm(c,inf)]);
         end
         [A,E,info2,Q,Z] = gklf([a b; c d],[e zeros(n,m);zeros(p,n+m)],tola);
         B = Q(n+1:end,:)'; C =  -Z(n+1:end,:); D = zeros(m,p);
         nreg = info2.mf+sum(info2.minf); nr1 = sum(info2.mr); nl = sum(info2.nl); ninv = nr1+nreg+nl;
         mr = n+m-ninv; % pr = n+p-ninv;
         if stabilize
            % make spurious zeros stable 
            if nr1
               i1 = 1:nr1; j1 = 1:mr; j2 = mr+1:mr+nr1; 
               F = gsfstab(A(i1,j2),E(i1,j2),A(i1,j1),poles,sdeg,struct('tol',tola)); 
               A(i1,j2) = A(i1,j2)+A(i1,j1)*F; C(:,j2) = C(:,j2)+C(:,j1)*F;
               fnorm = norm(F,'fro');
            end
            if nl
               i1 = nr1+nreg+1:nr1+nreg+nl; j1=n+m-nl+1:n+m; i2 = ninv+1:n+p;
               K = gsfstab(A(i1,j1)',E(i1,j1)',A(i2,j1)',poles,sdeg,struct('tol',tola)); K = K';
               A(i1,j1) = A(i1,j1)+K*A(i2,j1); B(i1,:) = B(i1,:)+K*B(i2,:); 
               fnorm = max(fnorm,norm(K,'fro'));
            end
         end
         nfp = nr1+nl;
         % select a square system
         ic = mr+1:n+m; ir = 1:ninv; 
         sysinv = gminreal(dss(A(ir,ic),B(ir,:),C(:,ic),D,E(ir,ic),Ts),tol);
      end
   end
elseif job == 1
   % compute an (1,2,3)-inverse
   if p == nrank
      % for full row rank, any (1,2)-inverse is an (1,2,3)-inverse
      [sysinv,info1] = grsol([ sys eye(p)],p,options);
      if ~mindeg
         nfp = info1.nr;
      end
   else
      opt = struct('tol',tol,'inner',true);
      % compute the full-rank factorization G = U*G1 with U inner
      [U,G1] = grange(sys,opt);
      if size(G1,1) ~= nrank
         warning('Inconsistent rank evaluations: check tolerances')
      end
      [G2,info1] = grsol([ G1 eye(nrank)],nrank,options);
      if ~mindeg
         nfp = info1.nr;
      end
      sysinv = gminreal(G2*U',tol);
   end
   tcond = info1.tcond;
   fnorm = info1.fnorm;
elseif job == 2
   % compute an (1,2,4)-inverse
   if m == nrank
      % for full column rank, any (1,2)-inverse is an (1,2,4)-inverse
      [sysinv,info1] = glsol([ sys; eye(m)],m,options);
      if ~mindeg
         nfp = info1.nl;
      end
   else
      opt = struct('tol',tol,'coinner',true);
      % compute the full-rank factorization G = G1*V with V coinner
      [V,G1] = gcrange(sys,opt);
      if size(G1,2) ~= nrank
         warning('Inconsistent rank evaluations: check tolerances')
      end
      [G2,info1] = glsol([ G1; eye(nrank)],nrank,options);
      if ~mindeg
         nfp = info1.nl;
      end
      sysinv = gminreal(V'*G2,tol);
   end
   tcond = info1.tcond;
   fnorm = info1.fnorm;
else
   % compute an (1,2,3,4)-inverse
   opt = struct('tol',tol,'inner',true,'coinner',true);
   if m == nrank
      % for full column rank, any (1,2,3)-inverse is an (1,2,3,4)-inverse
      % compute the full-rank factorization G = U*G2 with U inner and G2
      % square and invertible
      [U,G2] = grange(sys,opt);
      sysinv = gminreal(G2\U',tol);  
   elseif p == nrank
      % for full row rank, any (1,2,4)-inverse is an (1,2,3,4)-inverse
      % compute the full-rank factorization G = G2*V with V coinner and G2
      % square and invertible
      [V,G2] = gcrange(sys,opt);
      sysinv = gminreal(V'/G2,tol);  
   else
      % compute the full-rank factorization G = U*G1 with U inner
      [U,G1] = grange(sys,opt);
      if size(G1,1) ~= nrank
         warning('Inconsistent rank evaluations: check tolerances')
      end
      % compute the full-rank factorization G1 = G2*V with V coinner
      [V,G2] = gcrange(G1,opt);
      % compute the Moore-Penrose pseudo-inverse 
      sysinv = gminreal(V'*(G2\(U')),tol);   
   end
end
info.nrank = nrank;
info.tcond = tcond;
info.fnorm = fnorm;
info.nfp = nfp;

% end GINV

