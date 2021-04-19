function [sysr,sysx,info] = gcrange(sys,options)
%GCRANGE Coimage space basis of a transfer function matrix. 
% [SYSR,SYSX,INFO] = GCRANGE(SYS,TOL,OPTIONS) calculates, for a descriptor 
% system SYS = (A-lambda*E,B,C,D) with the transfer function matrix 
% G(lambda), the proper descriptor system SYSR = (Ar-lambda*Er,Br,Cr,Dr),
% with the full row rank proper transfer function matrix R(lambda) 
% such that 
%
%          Coimage G(lambda) = Coimage R(lambda) , 
% 
% and the descriptor system SYSX = (A-lambda*E,Bx,C,Dx) with the 
% full column rank transfer function matrix X(lambda), which satisfies
%
%        G(lambda) = X(lambda)*R(lambda) ,
%
% representing a full rank factorization of G(lambda). 
% The number of rows of R(lambda) is the normal rank of G(lambda). 
% The rows of R(lambda) form a rational basis of the coimage 
% space (row space) of the rational matrix G(lambda). 
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
% OPTIONS.coinner - specifies the option to compute co-inner basis 
%                   (can be used only in conjunction with 
%                   OPTIONS.zeros = 'none' or OPTIONS.zeros = 'unstable'):
%                   true  - compute an co-inner basis
%                   false - no co-inner basis is computed (default) 
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
% See also GRANGE.

% Author:      A. Varga, 10.11.2017.
% Revision(s): A. Varga, 21.08.2018, 15.09.2018. 
%
% Method: The range computation method, described in [1], is applied to the
% dual descriptor system realization corresponding to the transpose of the 
% rational matrix G(lambda). The underlying pencil reduction algorithm of 
% [2], has been adapted to deal with several zero selection options. 
% The computation of the involved Kronecker-like form is based on 
% the algorithm of [3]. 
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
   % coinner option
   if isfield(options,'coinner')
      inner = options.coinner;
      validateattributes(inner, {'logical'},{'binary'},'','OPTIONS.coinner') 
   else
      inner = false;
   end
   options.inner = inner;
end

if nargout < 2
   sysr = grange(sys.',options);
   sysr = sysr.';    
elseif nargout == 2
  [sysr,sysx] = grange(sys.',options);
  sysr = sysr.'; sysx = sysx.';
else
  [sysr,sysx,info] = grange(sys.',options);
  sysr = sysr.'; sysx = sysx.';
end

% end GCRANGE
end

