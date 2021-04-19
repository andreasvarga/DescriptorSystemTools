%sl_gstra  MEX-function for reduction to a detailed SVD-coordinate form.
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
% See also sl_gstra, sl_gstra4, sl_gstra7.



