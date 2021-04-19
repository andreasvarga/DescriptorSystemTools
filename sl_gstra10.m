%sl_gstra  MEX-function for reduction to a special controllability form.
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
% See also sl_gstra, sl_gstra8.



