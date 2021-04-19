%SL_GSTRA  MEX-function for descriptor system coordinate transformations.
%
% [At,Et,...] = sl_gstra(TASK,A,E,...) computes, according to the selected
% value TASK = 1-11, a specific similarity transformation on the
% matrices A, E, B, C of a descriptor system triple (A-lambda*E, B, C). 
% The posible choices of TASK are:
%     TASK = 1 : descriptor system scaling    
%     TASK = 2 : reduction to QR-coordinate form    
%     TASK = 3 : reduction to RQ-coordinate form  
%     TASK = 4 : reduction to SVD-coordinate form
%     TASK = 5 : reduction to a more detailed SVD-coordinate form
%     TASK = 6 : reduction to SVD-like coordinate form  
%     TASK = 7 : reduction to a more detailed SVD-like coordinate form  
%     TASK = 8 : reduction to controllability forms
%     TASK = 9 : reduction to observability forms
%     TASK = 10: reduction to a special controllability form
%                of the pair (A,B), where B = [B1,B2]
%     TASK = 11: reduction to Hessenberg coordinate form
% 
% See also sl_gstra1, sl_gstra2, sl_gstra3, sl_gstra4, sl_gstra5, 
% sl_gstra6, sl_gstra7, sl_gstra8, sl_gstra9, sl_gstra10, sl_gstra11  
% corresponding to the choices TASK = 1, 2, ..., 11, respectively. 

%  Authors: A. Varga and V. Sima, 11-01-2016.
%  Revision(s): 
%
%
%SL_GSTRA  MEX-function for descriptor system scaling.
%
% [At,Et,Bt,Ct,Dl,Dr] = sl_gstra(1,A,E,B,C,JOB,THRESH) scales, for a 
% descriptor system (A-lambda*E,B,C), the matrices of the system pencil
% 
%           S =  ( A  B ) - lambda * ( E  0 )                (1)
%                ( C  0 )            ( 0  0 )
% 
%  by applying diagonal similarity transformations Dl and Dr to determine
%  (At-lambda*Et,Bt,Ct) = (Dl*A*Dr - lambda*Dl*E*Dr, Dl*B, C*Dr) such that
%  the rows and columns of the matrices of the transformed system pencil
% 
%                diag(Dl,I) * S * diag(Dr,I)
% 
%  are as close in norm as possible. Scaling may reduce the 1-norms
%  of the matrices A, E, B, and C.
% 
%  JOB specifies the option for performing the balancing on the following
%  particular system pencils: 
%  JOB = 0: scale all matrices A, E, B, C using S in (1) (default);
%  JOB = 1: scale only A, E, B using S = ( A-lambda*E  B );
%  JOB = 2: scale only A, E, C using  S = ( A-lambda*E );
%                                         (     C      )
%  JOB = 3: scale using only A, E using S = A-lambda*E.
% 
% THRESH is an optional threshold on the size of the elements of Dl and Dr.
%
%
%SL_GSTRA  MEX-function for reduction to QR-coordinate form.
%
% [At,Et,Bt,Q] = sl_gstra(2,A,E,B) computes, for a descriptor pair 
% (A-lambda*E,B), an orthogonal transformation matrix Q such that 
% the transformed pair 
%       (At-lambda*Et, Bt) := (Q'*A-lambda*Q'*E, Q'*B) 
% has Et upper triangular. 
%
%
%SL_GSTRA  MEX-function for reduction to RQ-coordinate form.
%
% [At,Et,Ct,Z] = sl_gstra(3,A,E,C) computes, for a descriptor pair 
% (A-lambda*E,C), an orthogonal transformation matrix Z such that 
% the transformed pair 
%         (At-lambda*Et, Ct) := (A*Z-lambda*E*Z, C*Z) 
% has Et upper triangular. 
%
%
%SL_GSTRA  MEX-function for reduction to SVD-coordinate form.
%
% [At,Et,Bt,Ct,rankE,Q,Z] = sl_gstra(4,A,E,B,C,TOL) computes, for a 
% descriptor triple (A-lambda*E,B,C), orthogonal transformation matrices
% Q and Z such that the transformed triple
%
%       (At-lambda*Et,Bt,Ct) = (Q'*A*Z-lambda*Q'*E*Z,Q'*B,C*Z) 
%
% is in an SVD (singular value decomposition) coordinate form.  
% The resulting At and Et are in the form 
% 
%            ( A11  A12 )           ( Er  0 )
%       At = (          ) ,    Et = (       ) ,
%            ( A21  A22 )           (  0  0 )
% 
% where Er is a rankE x rankE invertible diagonal matrix having on its 
% diagonal the decreasingly ordered nonzero singular values of E.
%
% TOL is the relative tolerance to be used in determining the rank of E. 
% For TOL > 0, 1/TOL is a lower bound for the reciprocal condition numbers 
% of the diagonal matrices from the SVD decomposition of E. 
% (Default: prod(size(E))*eps). 
%
%
%SL_GSTRA  MEX-function for reduction to a detailed SVD-coordinate form.
%
% [At,Et,Bt,Ct,RANKS,Q,Z] = sl_gstra(5,A,E,B,C,TOL) computes, for a 
% descriptor triple (A-lambda*E,B,C), orthogonal transformation matrices
% Q and Z such that the transformed triple
%
%       (At-lambda*Et,Bt,Ct) = (Q'*A*Z-lambda*Q'*E*Z,Q'*B,C*Z) 
%
% is in a more detailed SVD (singular value decomposition) coordinate form. 
% The resulting At and Et are in the form 
% 
%         ( A11  A12 )   ( A11  *   *  )         ( Er  0  0 )
%    At = (          ) = (  *   Ar  0  ) ,  Et = ( 0   0  0 ) ,
%         ( A21  A22 )   (  *   0   0  )         ( 0   0  0 )
% 
% where Er is a rankE x rankE invertible diagonal matrix having on its 
% diagonal the decreasingly ordered nonzero singular values of E, A11 is 
% a rankE x rankE matrix, and Ar is a rankA22-by-rankA22 invertible 
% diagonal matrix, with decreasingly ordered positive diagonal elements. 
% RANKS is the two-dimensional vector RANKS = [ rankE rankA22 ].
%
% TOL is the relative tolerance to be used in determining the ranks of 
% E and A22. For TOL > 0, 1/TOL is a lower bound for the reciprocal  
% condition numbers of the diagonal matrices from the SVD decompositions  
% of E and A22. (Default: prod(size(E))*eps). 
%
%
%SL_GSTRA  MEX-function for reduction to SVD-like-coordinate form.
%
% [At,Et,Bt,Ct,rankE,Q,Z] = sl_gstra(6,A,E,B,C,TOL) computes, for a 
% descriptor triple (A-lambda*E,B,C), orthogonal transformation matrices
% Q and Z such that the transformed triple
%
%       (At-lambda*Et,Bt,Ct) = (Q'*A*Z-lambda*Q'*E*Z,Q'*B,C*Z) 
%
% is in an SVD (singular value decomposition)-like coordinate form.  
% The resulting At and Et are in the form 
% 
%            ( A11  A12 )           ( Er  0 )
%       At = (          ) ,    Et = (       ) ,
%            ( A21  A22 )           (  0  0 )
% 
% where Er is a rankE x rankE upper triangular invertible matrix.
%
% TOL is the relative tolerance to be used in determining the rank of E. 
% For TOL > 0, 1/TOL is a lower bound for the reciprocal condition numbers 
% of the leading upper triangular matrices in the QR decomposition of E. 
% (Default: prod(size(E))*eps). 
%
%
%SL_GSTRA  MEX-function for reduction to a detailed SVD-like coordinate form.
%
% [At,Et,Bt,Ct,RANKS,Q,Z] = sl_gstra(7,A,E,B,C,TOL) computes, for a 
% descriptor triple (A-lambda*E,B,C), orthogonal transformation matrices
% Q and Z such that the transformed triple
%
%       (At-lambda*Et,Bt,Ct) = (Q'*A*Z-lambda*Q'*E*Z,Q'*B,C*Z) 
%
% is in a more detailed SVD(singular value decomposition)-like coordinate  
% form. The resulting  At and Et are in the form 
% 
%         ( A11  A12 )   ( A11  *   *  )         ( Er  0  0 )
%    At = (          ) = (  *   Ar  0  ) ,  Et = ( 0   0  0 ) ,
%         ( A21  A22 )   (  *   0   0  )         ( 0   0  0 )
% 
% where Er is a rankE x rankE invertible upper triangular matrix, A11 is 
% a rankE x rankE matrix, and Ar is a rankA22-by-rankA22 invertible 
% upper triangular matrix.
% RANKS is the two-dimensional vector RANKS = [ rankE rankA22 ].
%
% TOL is the relative tolerance to be used in determining the ranks of 
% E and A22. For TOL > 0, 1/TOL is a lower bound for the reciprocal  
% condition numbers of the leading upper triangular matrices in the QR 
% decompositions of E and A22. (Default: prod(size(E))*eps). 
%
%
%SL_GSTRA  MEX-function for reduction to controllability forms.
%
% [At,Et,Bt,Ct,NCONT,TAU,Q,Z] = sl_gstra(8,A,E,B,C,JOB,TOL) computes, for  
% a descriptor system triple (A-lambda*E,B,C), orthogonal transformation 
% matrices Q and Z such that the transformed triple
%
%       (At-lambda*Et,Bt,Ct) = (Q'*A*Z-lambda*Q'*E*Z,Q'*B,C*Z) 
%
% is in a controllability form with the system matrices having the forms
% 
%          ( Ac  *  )         ( Ec  *  )         ( Bc )
%     At = (        ) ,  Et = (        ) ,  Bt = (    ) ,
%          ( 0  Anc )         ( 0  Enc )         ( 0  )
% 
%     Ct = ( Cc Cnc ) ,
% 
% where the NCONT-th order descriptor system triple (Ac-lambda*Ec,Bc,Cc)
% is finite and/or infinite controllable. The pencil Anc - lambda*Enc 
% is regular and contains the uncontrollable finite and/or infinite 
% eigenvalues of the pencil A-lambda*E.
% 
% For JOB = 0 or JOB = 2, the pencil ( Bc Ec-lambda*Ac ) has full
% row rank NCONT for all finite lambda and is in a staircase form, having
% Ac block upper triangular with TAU(i)xTAU(i) upper triangulardiagonal 
% blocks, and (Bc,Ec) in a staircase form with TAU(i)-by-TAU(i-1) full row 
% rank diagonal blocks, where TAU(0) is the number of columns of B. 
% 
% For JOB = 1, the pencil ( Bc Ac-lambda*Ec ) has full row rank NCONT for 
% all finite lambda and is in a staircase form, having Ec block upper 
% triangular with TAU(i)xTAU(i) upper triangular diagonal blocks, and
% (Bc,Ac) in a staircase form with TAU(i)-by-TAU(i-1) full row rank 
% diagonal blocks, where TAU(0) is the number of columns of B.
% 
% TOL is a relative tolerance to be used for rank determinations when 
% transforming the pair (A-lambda*E, B). (Default: prod(size(E))*eps). 
%
%
%SL_GSTRA  MEX-function for reduction to controllability forms.
%
% [At,Et,Bt,Ct,NOBSV,TAU,Q,Z] = sl_gstra(9,A,E,B,C,JOB,TOL) computes, for  
% a descriptor system triple (A-lambda*E,B,C), orthogonal transformation 
% matrices Q and Z such that the transformed triple
%
%       (At-lambda*Et,Bt,Ct) = (Q'*A*Z-lambda*Q'*E*Z,Q'*B,C*Z) 
%
% is in an observability form with the system matrices having the forms
% 
%          ( Ano  * )         ( Eno  * )         ( Bno )
%     At = (        ) ,  Et = (        ) ,  Bt = (     ) ,
%          ( 0   Ao )         ( 0   Eo )         ( Bo  )
% 
%     Ct = ( 0   Co ) ,
% 
% where the NOBSV-th order descriptor system triple (Ao-lambda*Eo,Bo,Co) is 
% finite and/or infinite observable. The pencil Ano - lambda*Eno is regular
% and contains the unobservable finite and/or infinite eigenvalues of the
% pencil A-lambda*E.
% 
% For JOB = 0 or JOB = 2, the pencil ( Eo-lambda*Ao ) has full
%                                    (      Co      )
% column rank NOBSV for all finite lambda and is in a staircase form, 
% having Ao block upper triangular with TAU(NT-i+1)xTAU(NT-i+1) upper 
% triangular diagonal blocks, and ( Eo ) in a staircase form with 
%                                 ( Co ) 
% TAU(NT-i+1)-by-TAU(NT-i) full column rank diagonal blocks, where 
% TAU(0) is the number of rows of C and NT is the dimension of vector TAU.
% 
% For JOB = 1, the pencil ( Eo-lambda*Ao ) has full column rank NOBSV 
%                         (      Co      )
% for all finite lambda and is in a staircase form, having Eo block upper 
% triangular with TAU(NT-i+1)xTAU(NT-i+1) upper triangular diagonal blocks, 
% and ( Ao ) in staircase form with TAU(NT-i+1)-by-TAU(NT-i) full column 
%     ( Co ) 
% rank diagonal blocks, where TAU(0) is the number of rows of C and NT is 
% the dimension of vector TAU.
% 
% TOL is a relative tolerance to be used for rank determinations when 
% transforming the pair (A-lambda*E, C). (Default: prod(size(E))*eps). 
%
%
%SL_GSTRA  MEX-function for reduction to a special controllability form.
%
% [At,Et,Bt,Ct,NCONT,TAU,Q,Z] = sl_gstra(10,A,E,B,C,M1,TOL) computes, for  
% a descriptor system triple (A-lambda*E,[B1,B2],C), orthogonal 
% transformation matrices Q and Z such that the transformed triple
%
%    (At-lambda*Et,[Bt1,Bt2],Ct) = (Q'*A*Z-lambda*Q'*E*Z,Q'*[B1,B2],C*Z) 
%
% is in a special controllability form with the system matrices having 
% the forms
% 
%          ( Ac  *  )         ( Ec  *  )                ( Bc1 Bc2 )
%     At = (        ) ,  Et = (        ) ,  [Bt1 Bt2] = (         ) ,
%          ( 0  Anc )         ( 0  Enc )                (  0   0  )
% 
%     Ct = ( Cc Cnc ) ,
% 
% where the NCONT-th order descriptor system triple 
% (Ac-lambda*Ec,[Bc1,Bc2],Cc) is finite controllable. The pencil 
% Anc-lambda*Enc is regular and contains the uncontrollable finite and 
% possibly some uncontrollable infinite eigenvalues of the pencil 
% A-lambda*E. M1 specifies the column dimension of the submatrices B1 of 
% B = [B1 B2] and Bt1 of Bt = [ Bt1 Bt2 ].  
% 
% The pencil ( Bc1 Bc2 Ac-lambda*Ec ) has full row rank NCONT for all 
% finite lambda and is in a staircase form, having Ec block upper 
% triangular with TAU(i)xTAU(i) upper triangular diagonal blocks, and 
% [Bc1 Bc2 Ac] in a staircase form with TAU(i)-by-TAU(i-2) full row rank 
% diagonal blocks, where TAU(-1) is M1, the number of columns of B1, and
% TAU(0) is the number of columns of B2. 
% 
% TOL is a relative tolerance to be used for rank determinations when 
% transforming the pair (A-lambda*E, B). (Default: prod(size(E))*eps). 
%
%
%SL_GSTRA  MEX-function for reduction to Hessenberg coordinate form.
%
% [At,Et,Bt,Ct,Q,Z] = sl_gstra(11,A,E,B,C) computes, for a descriptor 
% system triple (A-lambda*E,B,C), orthogonal transformation matrices 
% Q and Z such that the transformed triple
%
%    (At-lambda*Et,Bt,Ct) = (Q'*A*Z-lambda*Q'*E*Z,Q'*B,C*Z) 
%
% is in a Hessenberg coordinate form, with At upper Hessenberg and 
% Et upper triangular.  
% 

