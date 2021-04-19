%sl_glme  MEX-function for solving generalized Sylvester equations.
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
% See also sl_glme.  

% Authors: A. Varga and V. Sima, 11-01-2016.
% Revision(s): 
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

