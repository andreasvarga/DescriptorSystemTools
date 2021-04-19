%sl_glme  MEX-function for solving generalized Lyapunov equations.
%
% X = sl_glme(3,A,E,C,FLAG,TRANS) computes for FLAG(1) = 0 the solution X
% of the generalized continuous-time Lyapunov equation
%
%        op(A)'*X*op(E) + op(E)'*X*op(A) = C,   
%
% and for FLAG(1) = 1, the solution X of the generalized discrete-time  
% Lyapunov (or Stein) equation 
%
%       op(A)'*X*op(A) - op(E)'*X*op(E) = C,
%
% where op(M) = M, if TRANS = 0, and op(M) = M', if TRANS = 1.
% 
% FLAG is a two-dimensional vector, whose second component can be set 
% as follows:
%   FLAG(2) = 0 : if the pair (A,E) is in a general form;
%   FLAG(2) = 1 : if the pair (A,E) is in a generalized real Schur form.
%
% [X,SEP] = sl_glme(3,A,E,C,FLAG,TRANS) computes additionally SEP, 
% an estimate of the separation of the Lyapunov operator. 
%
% See also sl_glme.  

% Authors: A. Varga and V. Sima, 11-01-2016.
% Revision(s): 
%
% References:
% [1] Penzl, T.
%     Numerical solution of generalized Lyapunov equations.
%     Advances in Comp. Math., vol. 8, pp. 33-48, 1998.