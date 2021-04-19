%SL_GMINR  MEX-function for irreducible realization of descriptor systems.
%
% [Ar,Er,Br,Cr,Dr,REDINFO] = sl_gminr(1,A,E,B,C,D,TOL,JOB,SYSTYPE)
% computes for a system SYS with the descriptor system representation 
% (A-lambda*E,B,C,D) a reduced order system SYSR with the descriptor system 
% representation (Ar-lambda*Er,Br,Cr,Dr), such that SYS and SYSR have the 
% same transfer-function matrices. For standard systems, E = [] can be used. 
% According to the options specified by the combination of the input 
% parameters JOB and SYSTYPE, the resulting system representation of SYSR 
% has the following properties:
%    JOB = 0, SYSTYPE = 0 - irreducible (controllable and observable)
%    JOB = 0, SYSTYPE = 1 - finite controllable and observable
%    JOB = 0, SYSTYPE = 2 - infinite controllable and observable
%    JOB = 1, SYSTYPE = 0 - controllable
%    JOB = 1, SYSTYPE = 1 - finite controllable
%    JOB = 1, SYSTYPE = 2 - infinite controllable 
%    JOB = 2, SYSTYPE = 0 - observable
%    JOB = 2, SYSTYPE = 1 - finite observable
%    JOB = 2, SYSTYPE = 2 - infinite observable
% By calling sl_gminr with JOB+10 instead JOB, a preliminary scaling of the 
% initial descriptor system realization (A-lambda*E,B,C,D) is performed.
% TOL is an optional tolerance to be used for controllability and 
% observability tests (Default: internally computed if TOL = 0)
%
% The order reduction is achieved in four stages, by successively 
% eliminating uncontrollable and/or unobservable eigenvalues of the 
% pencil A-lambda*E, as follows:
%    Stage 1 - eliminate uncontrollable finite eigenvalues 
%    Stage 2 - eliminate uncontrollable nonzero finite eigenvalues and
%              uncontrollable infinite eigenvalues
%    Stage 3 - eliminate unobservable finite eigenvalues
%    Stage 4 - eliminate unobservable nonzero finite eigenvalues and
%              unobservable infinite eigenvalues 
% REDINFO is an eight-dimensional vector, which contains additional 
% informaton on the achieved order reductions at each stage and on the 
% structures of the resulting matrices Ar and Er:
% REDINFO(k)  > 0  - for k = 1..4, the reduction achieved in Stage k 
% REDINFO(k) <= 0  - for k = 1..4, no reduction achieved in Stage k 
% REDINFO(6) >= 0  - number of nonzero subdiagonals of Ar 
% REDINFO(7) >= 0  - number of nonzero subdiagonals of Er 
% REDINFO(k)  < 0  - for k = 5, 8 (contains no relevant information)
%
% See also sl_gminr.

% Author: A. Varga and V. Sima, 11-01-2016.
% Revision(s): 
%
% References
% 
% [1] A. Varga
%     Computation of Irreducible Generalized State-Space Realizations.
%     Kybernetika, vol. 26, pp. 89-106, 1990.
% 


