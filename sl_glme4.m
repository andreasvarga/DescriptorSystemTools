%sl_glme  MEX-function for solving generalized positive Lyapunov equations.
%
% X = sl_glme(4,A,E,C,FLAG,TRANS) computes for FLAG(1) = 0 the 
% upper triangular solution X of the generalized positive continuous-time
% Lyapunov equation
%
%  op(A)'*op(X)'*op(X)*op(E) + op(E)'*op(X)'*op(X)*op(A) =  -op(C)'*op(C), 
%
% and for FLAG(1) = 1 the upper triangular solution X of the generalized  
% positive discrete-time Lyapunov (or Stein) equation 
%
%  op(A)'*op(X)'*op(X)*op(A) - op(E)'*op(X)'*op(X)*op(E) = -op(C)'*op(C),
%
% where op(M) = M, if TRANS = 0, and op(M) = M', if TRANS = 1.
% For FLAG(1) = 0, the pair (A,E) must have all its generalized eigenvalues 
% with negative real parts, while for FLAG(1) = 0, the pair (A,E) must have
% all its generalized eigenvalues with moduli strictly less than 1. 
% 
% FLAG is a two-dimensional vector, whose second component can be set 
% as follows:
%   FLAG(2) = 0 : if the pair (A,E) is in a general form;
%   FLAG(2) = 1 : if the pair (A,E) is in a generalized real Schur form.
%
% See also sl_glme.  

% Authors: A. Varga and V. Sima, 11-01-2016.
% Revision(s): 
%
% References:
% [1] Hammarling, S.J.
%     Numerical solution of the stable, non-negative definite
%     Lyapunov equation.
%     IMA J. Num. Anal., 2, pp. 303-323, 1982.
% 
% [2] Penzl, T.
%     Numerical solution of generalized Lyapunov equations.
%     Advances in Comp. Math., vol. 8, pp. 33-48, 1998.


