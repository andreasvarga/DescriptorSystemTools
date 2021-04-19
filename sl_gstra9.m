%sl_gstra  MEX-function for reduction to controllability forms.
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
% See also sl_gstra, sl_gstra8.



