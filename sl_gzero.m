% SL_GZERO  MEX-function for the computation of system zeros.
% 
% [ZF,NI,MI,KRONR,INFE,KRONL] = SL_GZERO(A,E,B,C,D) computes for a LTI 
% descriptor system SYS = (A-lambda*E,B,C,D), the finite zeros of the 
% system pencil 
%                 ( A - lambda*E  B )
%                 (     C         D ) 
% and information on its infinite zeros and Kronecker structure. 
% Use E = [] or E = I for a standard state-space system. 
% Use D = [] for a system with no direct feedthrough matrix.
% For B = [], C = [] and D = [], the zeros (finite and infinite) of an 
% arbitrary rectangular pencil A - lambda*E are computed. 
% The computed output parameters are:                                
%   ZF      - finite invariant zeros  
%   NI      - two-dimensional vector, with NI(1) containing the number of 
%             infinite zeros and NI(2) containing the normal rank of the 
%             system pencil
%   MI      - multiplicities of the infinite zeros
%   KRONR   - right Kronecker indices
%   INFE    - dimensions of elementary infinite blocks 
%             (also the multiplicities of infinite eigenvalues)
%   KRONL   - left Kronecker indices
%
% [ZF,NI,MI,KRONR,INFE,KRONL] = SL_GZERO(A,E,B,C,D,TOL) uses the tolerance 
% TOL for rank determinations. 
%
% [ZF,NI,MI,KRONR,INFE,KRONL] = SL_GZERO(A,E,B,C,D,TOL,SCALE) uses the 
% optional parameter SCALE, to specify the scaling option as follows:
%    SCALE = 0, to perform preliminary scaling of the system pencil;
%    SCALE = 1, no scaling (default).
%

%  Authors:     A. Varga and V. Sima, 11-01-2016.
%  Revision(s): V. Sima, 01-10-2016.
%
%  References:
%
%  [1] Misra P., Van Dooren, P., Varga, A.:
%      Computation of structural invariants of generalized state-space 
%      systems. 
%      Automatica, vol. 30, pp. 1921-1936, 1994.
%
%  [2] Svaricek, F.
%      Computation of the structural invariants of linear
%      multivariable systems with an extended version of
%      the program ZEROS.
%      System & Control Letters, vol. 6, pp. 261-266, 1985.
%
%  [3] Emami-Naeini, A. and Van Dooren, P.
%      Computation of zeros of linear multivariable systems.
%      Automatica, vol. 18, pp. 415-430, 1982.


