%SL_GSEP MEX-function for descriptor system spectral separations.
%
% [At,Et,Bt,Ct,...] = SL_GSEP(TASK,A,E,B,C,...) computes, according to 
% the selected value TASK = 1-6, a specific similarity transformation 
% to separate the generalized eigenvalues of the pair (A,E), as follows: 
%   TASK = 1: finite-infinite separation    
%   TASK = 2: reduction to generalized Hessenberg form with 
%             finite-infinite separation    
%   TASK = 3: reduction to generalized real Schur form with 
%             finite-infinite separation    
%   TASK = 4: reduction to block diagonal generalized real Schur form
%             with finite-infinite separation    
%   TASK = 5: reduction to generalized Schur form with ordered 
%             finite-infinite separation    
%   TASK = 6: reduction to block diagonal generalized real Schur form
%             with ordered finite-infinite separation       
% 
% See also sl_gsep1, sl_gsep2, sl_gsep3, sl_gsep4, sl_gsep5, sl_gsep6,  
% corresponding to the choices TASK = 1, 2, ..., 6, respectively. 

%  Author:      A. Varga, 13-01-2016.
%  Revision(s):  
%
%SL_GSEP MEX-function for finite-infinite spectral separation.
%
% [At,Et,Bt,Ct,DIMS,TAU,Q,Z] = sl_gsep(1,A,E,B,C,JOB,TOL) computes, for
% the descriptor system triple (A-lambda*E,B,C), the orthogonal 
% transformation matrices Q and Z, such that the transformed descriptor 
% system triple
%        (At-lambda*Et,Bt,Ct) = (Q*A*Z-lambda*Q*E*Z,Q*B,C*Z) 
% has the matrices At and Et in the block upper triangular form 
%
%          ( A11  A12 )          ( E11  E12 )    
%     At = (          ) ,   Et = (          ) , 
%          (  0   A22 )          (  0   E22 )    
%    
% where the pairs (A11,E11) and (A22,E22) have distinct sets of
% generalized eigenvalues, as specified by the input parameter JOB:
%   JOB  = 0: (A11,E11) contains the infinite spectrum, with A11 invertible
%             and upper triangular and E11 nilpotent and upper triangular;
%             (A22,E22) contains the finite spectrum, with 
%             E22 invertible (default); 
%   JOB  = 1: (A11,E11) contains the finite spectrum,  with E11 invertible,
%             and (A22,E22) contains the infinite spectrum, with A22 
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
%
%SL_GSEP MEX-function for Hessenberg finite-infinite spectral separation.
%
% [At,Et,Bt,Ct,DIMS,TAU,Q,Z] = sl_gsep(2,A,E,B,C,JOB,TOL) computes, for
% the descriptor system triple (A-lambda*E,B,C), the orthogonal 
% transformation matrices Q and Z, such that the transformed descriptor 
% system triple
%        (At-lambda*Et,Bt,Ct) = (Q*A*Z-lambda*Q*E*Z,Q*B,C*Z) 
% has the pair (At,Et) in a generalized Hessenberg form with At and Et 
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
%             (A22,E22) contains the finite spectrum, with A22 in upper
%             Hessenberg form and E22 invertible and upper triangular; 
%   JOB  = 1: (A11,E11) contains the finite spectrum, with A11 in upper 
%             Hessenberg form and E11 invertible and upper triangular, and
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
%
%SL_GSEP MEX-function for Schur finite-infinite spectral separation.
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
%
%SL_GSEP MEX-function for block-diagonal finite-infinite spectral separation.
%
% [At,Et,Bt,Ct,DIMS,TAU,Q,Z] = sl_gsep(4,A,E,B,C,JOB,TOL) computes, for
% the descriptor system triple (A-lambda*E,B,C), the transformation 
% matrices Q and Z, such that the transformed descriptor system triple
%        (At-lambda*Et,Bt,Ct) = (Q*A*Z-lambda*Q*E*Z,Q*B,C*Z) 
% has the pair (At,Et) in a generalized real Schur form with At and Et 
% having the block-diagonal forms 
%
%          ( A11   0  )          ( E11   0  )    
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
%
%SL_GSEP MEX-function for ordered Schur finite-infinite spectral separation.
%
% [At,Et,Bt,Ct,DIMS,TAU,Q,Z] = sl_gsep(5,A,E,B,C,JOB,TOL,DISC,SMARG) 
% computes, for the descriptor system triple (A-lambda*E,B,C), the  
% orthogonal transformation matrices Q and Z, such that the transformed  
% descriptor system triple
%        (At-lambda*Et,Bt,Ct) = (Q*A*Z-lambda*Q*E*Z,Q*B,C*Z) 
% has the pair (At,Et) in a generalized real Schur form with At and Et 
% having the block upper-triangular forms 
%
%          ( A11  A12  A13 )          ( E11  E12  E13 )    
%     At = (  0   A22  A23 ) ,   Et = (  0   E22  E23 ) , 
%          (  0    0   A33 )          (  0    0   E33 )    
%    
% where the pairs (A11,E11), (A22,E22) and (A33,A33) have distinct sets of
% generalized eigenvalues, as specified by the input parameter JOB:
%   JOB  = 0: (A11,E11) contains the infinite spectrum, with A11 invertible
%             and upper triangular and E11 nilpotent and upper triangular;
%             (A22,E22) contains the unstable finite spectrum, with A22 
%             upper quasi-triangular and E22 invertible and 
%             upper triangular;  and
%             (A33,E33) contains the stable finite spectrum, with A33 
%             upper quasi-triangular and E33 invertible and 
%             upper triangular;  
%   JOB  = 1: (A11,E11) contains the stable finite spectrum, with A11 
%             upper quasi-triangular and E11 invertible and 
%             upper triangular;  
%             (A22,E22) contains the unstable finite spectrum, with A22 
%             upper quasi-triangular and E22 invertible;  and
%             (A33,E33) contains the infinite spectrum, with A33 invertible
%             and upper triangular and E33 nilpotent and upper triangular;
%   JOB  = 2: (A11,E11) contains the infinite spectrum, with A11 invertible
%             and upper triangular and E11 nilpotent and upper triangular;
%             (A22,E22) contains the stable finite spectrum, with A22 
%             upper quasi-triangular and E22 invertible and 
%             upper triangular;  and 
%             (A33,E33) contains the unstable finite spectrum with A33 
%             upper quasi-triangular and E33 invertible and 
%             upper triangular;  
%   JOB  = 3: (A11,E11) contains the unstable finite spectrum, with A11 
%             upper quasi-triangular and E11 invertible and 
%             upper triangular;   
%             (A22,E22) contains the stable finite spectrum, with A22 
%             upper quasi-triangular and E22 invertible and 
%             upper triangular;  and
%             (A33,E33) contains the infinite spectrum, with A33 invertible
%             and upper triangular and E33 nilpotent and upper triangular.
% TOL specifies an optional relative tolerance for rank computations 
%   (default: internally computed). 
% DISC and SMARG specifies the stability domain for the generalized 
%   eigenvalues of the pair (A,E), as follows:
%   DISC = 0: complex numbers with real parts at most SMARG (default);
%   DISC = 1: complex numbers with moduli at most SMARG. 
%   The default values of SMARG are -sqrt(eps) for DISC = 0 and  
%   1-sqrt(eps) for DISC = 1.
%
% The output parameters DIMS and TAU provide additional information on the 
% structure of the resulting pair (At,Et), as follows: 
% DIMS is a vector [n1,n2,n3,rankE], where n1 is the order of the pencil 
%   A11-lambda*E11, n2 is the order of the pencil A22-lambda*E22, n3 is the   
%   order of the pencil A33-lambda*E33, and rankE is the rank of E;
% TAU contains the dimensions of the diagonal blocks in the subpencil 
%   which contains the infinite eigenvalues. The number of 
%   simple (nondynamic) infinite generalized eigenvalues of the pair (A,E)
%   is contained in TAU(END) for JOB = 0 or 2 and TAU(1) for JOB = 1 or 3. 
%
%
%SL_GSEP MEX-function for ordered block-diagonal spectral separation.
%
% [At,Et,Bt,Ct,DIMS,TAU,Q,Z] = sl_gsep(6,A,E,B,C,JOB,TOL,DISC,SMARG) 
% computes, for the descriptor system triple (A-lambda*E,B,C), the  
% transformation matrices Q and Z, such that the transformed descriptor 
% system triple
%        (At-lambda*Et,Bt,Ct) = (Q*A*Z-lambda*Q*E*Z,Q*B,C*Z) 
% has the pair (At,Et) in a generalized real Schur form with At and Et 
% having a block-diagonal form 
%
%          ( A11  A12  A13 )          ( E11  E12  E13 )    
%     At = (  0   A22  A23 ) ,   Et = (  0   E22  E23 ) , 
%          (  0    0   A33 )          (  0    0   E33 )    
%  
% with either A13, A23, E13, E23 zero or A12, A13, E12, E13 zero, and
% where the pairs (A11,E11), (A22,E22) and (A33,A33) have distinct sets of
% generalized eigenvalues, as specified by the input parameter JOB:
%   JOB  = 0: A13, A23, E13, E23 are zero, and 
%             (A11,E11) contains the infinite spectrum, with A11 invertible
%             and upper triangular and E11 nilpotent and upper triangular;
%             (A22,E22) contains the unstable finite spectrum, with A22 
%             upper quasi-triangular and E22 invertible and 
%             upper triangular;  and
%             (A33,E33) contains the stable finite spectrum, with A33 
%             upper quasi-triangular and E33 invertible and 
%             upper triangular;
%   JOB  = 1: A12, A13, E12, E13 are zero and 
%             (A11,E11) contains the stable finite spectrum, with A11 
%             upper quasi-triangular and E11 invertible and 
%             upper triangular;  
%             (A22,E22) contains the unstable finite spectrum, with A22 
%             upper quasi-triangular and E22 invertible;  and
%             (A33,E33) contains the infinite spectrum, with A33 invertible
%             and upper triangular and E33 nilpotent and upper triangular;
%   JOB  = 2: A13, A23, E13, E23 are zero, and 
%             (A11,E11) contains the infinite spectrum, with A11 invertible
%             and upper triangular and E11 nilpotent and upper triangular;
%             (A22,E22) contains the stable finite spectrum, with A22 
%             upper quasi-triangular and E22 invertible and 
%             upper triangular;  and 
%             (A33,E33) contains the unstable finite spectrum with A33 
%             upper quasi-triangular and E33 invertible and 
%             upper triangular;  
%   JOB  = 3: A12, A13, E12, E13 are zero and 
%             (A11,E11) contains the unstable finite spectrum, with A11 
%             upper quasi-triangular and E11 invertible and 
%             upper triangular;   
%             (A22,E22) contains the stable finite spectrum, with A22 
%             upper quasi-triangular and E22 invertible and 
%             upper triangular;  and
%             (A33,E33) contains the infinite spectrum, with A33 invertible
%             and upper triangular and E33 nilpotent and upper triangular.
% TOL specifies an optional relative tolerance for rank computations 
%   (default: internally computed). 
% DISC and SMARG specifies the stability domain for the generalized 
%   eigenvalues of the pair (A,E), as follows:
%   DISC = 0: complex numbers with real parts at most SMARG (default);
%   DISC = 1: complex numbers with moduli at most SMARG. 
%   The default values of SMARG are -sqrt(eps) for DISC = 0 and  
%   1-sqrt(eps) for DISC = 1.
%
% The output parameters DIMS and TAU provide additional information on the 
% structure of the resulting pair (At,Et), as follows: 
% DIMS is a vector [n1,n2,n3,rankE], where n1 is the order of the pencil 
%   A11-lambda*E11, n2 is the order of the pencil A22-lambda*E22, n3 is the   
%   order of the pencil A33-lambda*E33, and rankE is the rank of E;
% TAU contains the dimensions of the diagonal blocks in the subpencil 
%   which contains the infinite eigenvalues. The number of 
%   simple (nondynamic) infinite generalized eigenvalues of the pair (A,E)
%   is contained in TAU(END) for JOB = 0 or 2 and TAU(1) for JOB = 1 or 3. 
%


