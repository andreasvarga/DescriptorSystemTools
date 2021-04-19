%sl_gstra  MEX-function for reduction to controllability forms.
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
% See also sl_gstra, sl_gstra9, sl_gstra10.



