%SL_GMINR  MEX-function for ellimination of non-dynamic modes.
%
% [Ar,Er,Br,Cr,Dr,REDINFO] = sl_gminr(2,A,E,B,C,D,TOL,JOB) computes for a 
% system SYS with a descriptor system representation (A-lambda*E,B,C,D) a 
% reduced order system SYSR with the descriptor system representation 
% (Ar-lambda*Er,Br,Cr,Dr), such that SYS and SYSR have the same
% transfer-function matrices and the pencil Ar-lambda*Er has no simple
% infinite eigenvalues (non-dynamic modes). For a standard system with
% E = [] or E = I, no computation is performed. The following options can 
% be specified via the input parameter JOB for the resulting shape of Er: 
%    JOB = 0 - Er is either unreduced or upper triangular;
%    JOB = 1 - Er is diagonal with the leading nonzero block identity.
% TOL is an optional tolerance to be used for rank determinations
% (Default: internally computed if TOL = 0)
%
% The order reduction is achieved by eliminating the simple infinite 
% eigenvalues (non-dynamic modes) of the pencil A-lambda*E.
% REDINFO is an eight-dimensional vector, which contains additional 
% informaton on the achieved order reduction, on the shape of the 
% resulting matrices Ar and Er, and on the rank of E:
% REDINFO(k) <= 0 - for k = 1..4
% REDINFO(5)  > 0 - number of removed non-dynamic modes 
% REDINFO(5) <= 0 - no order reduction performed 
% REDINFO(6) >= 0 - number of nonzero subdiagonals of Ar 
% REDINFO(7) >= 0 - number of nonzero subdiagonals of Er 
% REDINFO(8) >= 0 - rank of E 
% REDINFO(8)  < 0 - rank of E not computed
%
% See also sl_gminr.

% Author: A. Varga and V. Sima, 11-01-2016.
% Revision(s): 
%
