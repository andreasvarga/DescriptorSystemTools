%sl_gstra MEX-function for reduction to a detailed SVD-like coordinate form.
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
% See also sl_gstra, sl_gstra6.
