%SL_KLF  MEX-function for pencil reduction to Kronecker-like forms.
%
% [At,Et,Q,Z,IMUK,INUK,IMUKINF] = SL_KLF(A,E,TOL,MODE,QZ) computes for a 
% a linear pencil A - sE, the orthogonal transformations Q and Z 
% such that the transformed pencil 
%        At - sEt := Q'(A - sE)Z 
% is in an upper block triangular Kronecker-like form, in accordance with
% the option parameter MODE. TOL is an optional tolerance used for 
% rank determinations. QZ specifies which of the matrices Q and Z are 
% accumulated as follows: 
%     QZ = 0, if Q and Z are not accumulated and Q and Z are set [];
%     QZ = 1, if Q is accumulated and Z is set [];
%     QZ = 2, if Z is accumulated and Q is set [];
%     QZ = 3, if both Q and Z are accumulated (default).
% According to the value of MODE, the following reductions are performed:
%
%     If MODE = 0, then the matrices A and E are transformed into the
%     following Kronecker-like form by unitary transformations Q1
%     and Z1 :
%
%                    | A(r,inf)-sE(r,inf) |       X        |
%      Q1'(A-sE)Z1 = |--------------------|----------------|.   (1)
%                    |            O       | A(f,l)-sE(f,l) |
%
%     The pencil A(r,inf)-sE(r,inf) is in staircase form, and it
%     contains all Kronecker right indices and infinite elementary
%     divisors of the pencil A-sE. The pencil A(f,l)-sE(f,l) contains all
%     Kronecker left indices and finite elementary divisors of A-sE.
%     Note: X is a pencil.
%
%     If MODE = 1, then the submatrices having full row and column
%     rank in the pencil A(r,inf)-sE(r,inf) in (1) are
%     triangularized by applying unitary transformations Q2 and Z2 to
%     Q1'*(A-sE)*Z1.
%
%     If MODE = 2 (default), then the pencil A(r,inf)-sE(r,inf) in 
%     (1) is separated into A(r)-sE(r) and A(inf)-sE(inf) by applying
%     unitary transformations Q3 and Z3 to Q2'*Q1'*(A-sE)*Z1*Z2.
%
%     This gives
%
%                 | A(r)-sE(r) |        X       |       X        |
%                 |------------|----------------|----------------|
%                 |        O   | A(inf)-sE(inf) |       X        |
%     Q'(A-sE)Z = |=============================|================|, (2)
%                 |                             |                |
%                 |            O                | A(f,l)-sE(f,l) |
%              
%     where Q = Q1*Q2*Q3 and Z = Z1*Z2*Z3.
%
% The output parameters IMUK, INUK, and IMUKINF contain information on the 
% structure of the computed subpencils, as follows: 
%   IMUK    - contains the column dimensions of the submatrices having 
%             full column rank in the pencil A(x)-sE(x),
%             where  x = r,inf  if MODE = 0 or 1,
%             or     x = r      if MODE = 2.
%   INUK    - contains the row dimensions of the submatrices having 
%             full row rank in the pencil A(x)-sE(x),
%             where  x = r,inf  if MODE = 0 or 1,
%             or     x = r      if MODE = 2.
%   IMUKIN -  if MODE = 2, then the elements of this vector contain 
%             the dimensions of the square diagonal submatrices in the  
%             pencil A(inf)-sE(inf). 
%

% Authors: A. Varga and V. Sima, 11-01-2016.
% Revision(s): 
%
% References
% [1] Beelen, Th. and Van Dooren, P.
%     An improved algorithm for the computation of Kronecker's
%     canonical form of a singular pencil.
%     Linear Algebra and Applications, vol. 105, pp. 9-65, 1988.



