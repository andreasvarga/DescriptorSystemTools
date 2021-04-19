function [sysr,sysx,info] = grange(sys,options)
%GRANGE Range space basis of a transfer function matrix. 
% [SYSR,SYSX,INFO] = GRANGE(SYS,TOL,OPTIONS) calculates, for a descriptor 
% system SYS = (A-lambda*E,B,C,D) with the transfer function matrix  
% G(lambda), the proper descriptor system SYSR = (Ar-lambda*Er,Br,Cr,Dr), 
% with the full column rank proper transfer function matrix R(lambda) 
% such that 
%
%          Range G(lambda) = Range R(lambda) , 
% 
% and the descriptor system SYSX = (A-lambda*E,B,Cx,Dx) with the 
% full row rank transfer function matrix X(lambda), which satisfies
%
%        G(lambda) = R(lambda)*X(lambda) ,
%
% representing a full rank factorization of G(lambda). 
% The number of columns of R(lambda) is the normal rank of G(lambda). 
% The columns of R(lambda) form a rational basis of the range (or image)
% space of the rational matrix G(lambda). 
% A selected set of zeros of G(lambda) are included as zeros of R(lambda). 
%
% The OPTIONS structure allows to specify various user options, 
% as follows:
% OPTIONS.tol     - tolerance for rank determinations 
%                   (Default: internally computed)
% OPTIONS.zeros   - option for the selection of zeros to be included in
%                   the computed range basis:
%                   'none'       - include no zeros (Default) 
%                   'all'        - include all zeros of SYS
%                   'unstable'   - include all unstable zeros of SYS
%                   's-unstable' - include all strictly unstable zeros 
%                                  of SYS, both finite and infinite
%                   'stable'     - include all stable zeros of SYS
%                   'finite'     - include all finite zeros of SYS
%                   'infinite'   - include all infinite zeros of SYS
% OPTIONS.offset  - stability boundary offset OFFSET, to be used to assess 
%                   the finite zeros which belong to the boundary of the 
%                   stability domain as follows: 
%                   in the continuous-time case, these are the finite zeros
%                      having real parts in the interval [-OFFSET, OFFSET]; 
%                   in the discrete-time case, these are the finite zeros 
%                      having moduli in the interval [1-OFFSET,1+OFFSET]; 
%                   (Default: OFFSET = sqrt(eps)).
% OPTIONS.inner   - specifies the option to compute an inner basis 
%                   (can be used only in conjunction with 
%                   OPTIONS.zeros = 'none' or OPTIONS.zeros = 'unstable'):
%                   true  - compute an inner basis
%                   false - no inner basis is computed (default) 
% OPTIONS.balance - specifies the balancing option for the Riccati 
%                   equation solvers: 
%                   true  - apply balancing (Default)
%                   false - disable balancing  
%  
% INFO is a structure containing additional information, as follows: 
% INFO.nrank  - normal rank r of the transfer function matrix G(lambda);
% INFO.nfuz   - number of finite unstable zeros of SYS lying
%               on the boundary of the stability region;
% INFO.niuz   - number of infinite zeros of SYS in the continuous-time case
%               and 0 in the discrete-time case
% INFO.ricrez - diagnosis flag, as provided provided by the generalized 
%               Riccati equation solvers care and dare; if non-negative, 
%               this value represents the Frobenius norm of relative 
%               residual of the Riccati equation, while a negative value 
%               indicates failure of solving the Riccati equation. 
%
% See also GCRANGE.

% Copyright 2017-2018 A. Varga
% Author:      A. Varga, 04.07.2017.
% Revision(s): A. Varga, 10.11.2017, 21.08.2018, 15.09.2018, 07.07.2019. 
%
% Method: The range computation method is described in [1] and is based on 
% the reduction algorithm of [2], which has been adapted to deal with 
% several zero selection options. The computation of the involved 
% Kronecker-like form is based on the algorithm of [3]. 
%
% References
% [1] Varga, A.
%     A note on computing the range of rational matrices. 
%     arXiv:1707.0048, https://arxiv.org/abs/1707.0048, 2017.
% [2] Oara, C. 
%     Constructive solutions to spectral and inner–outer factorizations 
%     with respect to the disk. Automatica, 41:1855–1866, 2005.
% [3] Beelen, Th. and Van Dooren, P.
%     An improved algorithm for the computation of Kronecker's
%     canonical form of a singular pencil.
%     Linear Algebra and Applications, vol. 105, pp. 9-65, 1988.



narginchk(1,2)
nargoutchk(0,3)

if ~isa(sys,'ss')
   error('The input system SYS must be a state-space system object')
end

if nargin < 2
   options = struct('zeros','none','inner',false);
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

% zeros selection option
if isfield(options,'zeros')
   zerosopt = options.zeros;
   validateattributes(zerosopt, {'char'},{'nonempty'},'','OPTIONS.zeros') 
else
   zerosopt = 'none';
end

switch zerosopt
    case 'none'
       job = 0; 
    case 'unstable'
       job = 1;
    case 'all'
       job = 2;
    case 'infinite'
       job = 3; 
    case 'finite'
       job = 4; 
    case 'stable'
       job = 5;
    case 's-unstable'
       job = 6;
    otherwise
       error('No such zero selection option')
end

% offset for the stability region boundary
if isfield(options,'offset')
   offset = options.offset;
   validateattributes(offset, {'double'},{'real','scalar','>',0,'<',1},'','OPTIONS.offset') 
else
   offset = sqrt(eps);
end

% inner option
if isfield(options,'inner')
   inner = options.inner;
   validateattributes(inner, {'logical'},{'binary'},'','OPTIONS.inner') 
else
   inner = false;
end

if inner && job > 1
   warning('No inner range computation possible for the selected zeros')
   inner = false;
end

% balancing option for Riccati equation solvers
if isfield(options,'balance')
   balopt = options.balance;
   validateattributes(balopt, {'logical'},{'binary'},'','OPTIONS.balance')
   if balopt
      balance = 'balance';
   else
      balance = 'nobalance';
   end    
else
   balance = 'balance';
end

%  Reduce the system matrix pencil to the special Kronecker-like form
%
%                 ( Arg-lambda*Erg     *            *   *   )
%                 (   0            Abl-lambda*Ebl  Bbl  *   )
%  At-lambda*Et = (   0                0            0   Bn  )
%                 (-----------------------------------------)
%                 (   0               Cbl          Dbl  *   )

% enforce impulse controllability
sys = gir(sys,tol,'infinite_contr'); 
[At,Et,dimsc,~,Z] = gsklf(sys,tol,zerosopt,offset,'noQ');

%  Select matrices for the extracted full column rank factor
%
[a,b,~,d,e,Ts] = dssdata(sys);
discr = (Ts ~= 0);
[n,m] = size(b); p = size(d,1);
mc    = dimsc(1); 
nbl   = dimsc(2);
nrank = dimsc(3); 
nsinf = dimsc(4);
nfuz   = dimsc(5);
niuz   = dimsc(6);
nm    = n+m;
nc    = n-nbl-nsinf; 


ia = nc+1:nc+nbl; ja = mc+1:mc+nbl; jb = mc+nbl+1:mc+nbl+nrank;
ic = n+1:n+p; 
A = At(ia,ja); E = Et(ia,ja); B = At(ia,jb); C = At(ic,ja); D = At(ic,jb);
if inner 
   % enforce finite stabilizability
   [A,E,B,C,D] = sl_gminr(1,A,E,B,C,D,tol,1,1);
   if discr
      [X,~,F,ricrez] = dare(A,B,C.'*C,D.'*D,C.'*D,E,balance); F = -F;
      if ricrez < 0 || isnan(ricrez)
         if ricrez == -1 
            error('Solution of the DARE failed: Symplectic matrix has eigenvalues on the unit circle')
         else
            error('Solution of the DARE failed: no finite stabilizing solution exists')
         end
      end
      H = chol(D'*D+B'*X*B); 
   else
      [~,~,F,ricrez] = care(A,B,C.'*C,D.'*D,C.'*D,E,balance); F = -F;
      if ricrez < 0 || isnan(ricrez)
         if ricrez == -1 
            error('Solution of the CARE failed: Hamiltonian has jw-axis eigenvalues')
         else
            error('Solution of the CARE failed: no finite stabilizing solution exists')
         end
      end
      [~,H] = qr(D); H = H(1:nrank,:); 
   end
   sysr = dss(A+B*F,B/H,C+D*F,D/H,E,Ts); 
   if nargout >= 2
      CDt = [zeros(nrank,mc) -H*F H zeros(nrank,nsinf)]*Z'; 
      sysx = dss(a,b,CDt(:,1:n),CDt(:,n+1:nm),e,Ts);
   end
else
   ricrez = 0;
   sysr = dss(A,B,C,D,E,Ts);
   if nargout >= 2
      sysx = dss(a,b,Z(1:n,jb)',Z(n+1:nm,jb)',e,Ts);
   end
end

if nargout == 3
   info = struct('nrank',mric,'nfuz',nfuz,'niuz',niuz,'ricrez',ricrez);
end

% end GRANGE
end

