%sl_gsep MEX-function for Schur finite-infinite spectral separation.
%
% [At,Et,Bt,Ct,DIMS,TAU,Q,Z] = sl_gsep(3,A,E,B,C,JOB,TOL) computes, for
% the descriptor system triple (A-lambda*E,B,C), the orthogonal 
% transformation matrices Q and Z, such that the transformed descriptor 
% system triple
%        (At-lambda*Et,Bt,Ct) = (Q*A*Z-lambda*Q*E*Z,Q*B,C*Z) 
% has the pair (At,Et) in a generalized real Schur form with At and Et 
% having the block upper triangular forms 
%
%          ( A11  A12 )          ( E11  E12 )    
%     At = (          ) ,   Et = (          ) , 
%          (  0   A22 )          (  0   E22 )    
%    
% where the pairs (A11,E11) and (A22,E22) have distinct sets of
% generalized eigenvalues, as specified by the input parameter JOB:
%   JOB  = 0: (A11,E11) contains the infinite spectrum, with A11 invertible
%             and upper triangular and E11 nilpotent and upper triangular;
%             (A22,E22) contains the finite spectrum, with A22 upper
%             quasi-triangular and E22 invertible and upper triangular; 
%   JOB  = 1: (A11,E11) contains the finite spectrum, with A11 upper 
%             quasi-triangular and E11 invertible and upper triangular, and
%             (A22,E22) contains the infinite spectrum, with A22 
%             invertible and upper triangular and E11 nilpotent and 
%             upper triangular. 
% TOL specifies an optional relative tolerance for rank computations 
%   (default: internally computed). 
%
% The output parameters DIMS and TAU provide additional information on the 
% structure of the resulting pair (At,Et), as follows: 
% DIMS is a vector [n1,0,n2,rankE], where n1 is the order of the pencil 
%   A11-lambda*E11, n2 is the order of the pencil A22-lambda*E22, and 
%   rankE is the rank of E;
% TAU contains the dimensions of the diagonal blocks in the subpencil 
%   which contains the infinite eigenvalues. The number of 
%   simple (nondynamic) infinite generalized eigenvalues of the pair (A,E)
%   is contained in TAU(END) for JOB = 0 and TAU(1) for JOB = 1. 
%
% See also sl_gsep. 

%  Author:      A. Varga, 13-01-2016.
%  Revision(s):  

