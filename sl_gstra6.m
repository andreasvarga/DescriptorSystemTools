%sl_gstra  MEX-function for reduction to SVD-like-coordinate form.
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
% See also sl_gstra, sl_gstra4, sl_gstra7.



