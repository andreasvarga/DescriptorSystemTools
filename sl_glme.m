%SL_GLME  MEX-function for solving generalized linear matrix equations.
%
% [X,...] = sl_glme(TASK,A,...) solves, according to the selected value of 
% TASK, a class of linear matrix equations. 
% The posible choices of TASK are:
%     TASK = 1 : solve generalized Sylvester equations;   
%     TASK = 3 : solve generalized Lyapunov equations;
%     TASK = 4 : solve generalized positive Lyapunov equations.
%
% See also sl_glme1, sl_glme3, sl_glme4 corresponding to the choices 
% TASK = 1, 3, 4, respectively. 

%  Authors: A. Varga and V. Sima, 11-01-2016.
%  Revision(s): 
%
%
%SL_GLME  MEX-function for solving generalized Sylvester equations.
%
% [X,Y] = sl_glme(1,A,D,B,E,C,F,FLAG,TRANS) computes for TRANS = 0 the  
% solution (X,Y) of the generalized Sylvester equation
%
%       A*X - Y*B = C,
%       D*X - Y*E = F,
%
% and for TRANS = 1 the solution (X,Y) of the generalized Sylvester 
% equation 
%
%       A'*X + D'*Y = C,
%       X*B' + Y*E' = -F.
% 
% FLAG is an optional two-dimensional vector, whose components can be set 
% as follows:
%   FLAG(1) = 0 : if the pair (A,D) is in a general form;
%   FLAG(1) = 1 : if the pair (A,D) is in a generalized real Schur form;
%   FLAG(2) = 0 : if the pair (B,E) is in a general form;
%   FLAG(2) = 1 : if the pair (B,E) is in a generalized real Schur form.
%
% [X,Y,DIF] = sl_glme(1,A,D,B,E,C,F,,FLAG,TRANS) computes additionally DIF, 
% an estimate of the separation between the matrix pairs (A,D) and (B,F). 
%
% References:
% 
% [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
%     for Solving the Generalized Sylvester Equation and Estimating the
%     Separation between Regular Matrix Pairs, Report UMINF - 93.23,
%     Department of Computing Science, Umea University, S-901 87 Umea,
%     Sweden, December 1993, Revised April 1994, Also as LAPACK Working
%     Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,
%     No 1, 1996.
% 
% [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester
%     Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal.
%     Appl., 15(4):1045-1060, 1994
% 
% [3] B. Kagstrom and L. Westin, Generalized Schur Methods with
%     Condition Estimators for Solving the Generalized Sylvester
%     Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7,
%     July 1989, pp 745-751.
%
%
%SL_GLME  MEX-function for solving generalized Lyapunov equations.
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
% References:
% [1] Penzl, T.
%     Numerical solution of generalized Lyapunov equations.
%     Advances in Comp. Math., vol. 8, pp. 33-48, 1998.
%
%
%SL_GLME  MEX-function for solving generalized positive Lyapunov equations.
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
% References:
% [1] Hammarling, S.J.
%     Numerical solution of the stable, non-negative definite
%     Lyapunov equation.
%     IMA J. Num. Anal., 2, pp. 303-323, 1982.
% 
% [2] Penzl, T.
%     Numerical solution of generalized Lyapunov equations.
%     Advances in Comp. Math., vol. 8, pp. 33-48, 1998.





