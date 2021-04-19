%sl_gstra  MEX-function for reduction to Hessenberg coordinate form.
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
% See also sl_gstra.



