%SL_GMINR  MEX-function for minimal realization of descriptor systems.
%
% [Ar,Er,Br,Cr,Dr,...] = sl_gminr(TASK,A,E,B,C,D,...) computes for a 
% system SYS with the descriptor representation (A-lambda*E,B,C,D) 
% a reduced order system SYSR with a descriptor representation
% (Ar-lambda*Er,Br,Cr,Dr), such that SYS and SYSR have the same
% transfer-function matrices. According to the selected value of TASK, the 
% resulting system representation of SYSR has the following properties:
%     TASK = 1 : irreducible or finite/infinite controllable or
%                finite/infinite observable   
%     TASK = 2 : without non-dynamic modes   
%     TASK = 3 : without non-dynamic modes  and either irreducible or 
%                finite/infinite controllable or finite/infinite observable 
% 
% See also sl_gminr1, sl_gminr2, sl_gminr3 corresponding to the choices 
% TASK = 1, 2, 3, respectively. 

%  Authors: A. Varga and V. Sima, 11-01-2016.
%  Revision(s): 
%
% References
% 
% [1] A. Varga
%     Computation of Irreducible Generalized State-Space Realizations.
%     Kybernetika, vol. 26, pp. 89-106, 1990.
% 
%
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
%
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
%
%SL_GMINR  MEX-function for minimal realization of descriptor systems.
%
% [Ar,Er,Br,Cr,Dr,REDINFO] = sl_gminr(3,A,E,B,C,D,TOL,JOB,SYSTYPE)
% computes for a system SYS with a descriptor system representation 
% (A-lambda*E,B,C,D) a reduced order system SYSR with the descriptor system 
% representation (Ar-lambda*Er,Br,Cr,Dr), such that SYS and SYSR have the 
% same transfer-function matrices and the pencil Ar-lambda*Er has no 
% simple infinite eigenvalues (non-dynamic modes).  
% For a standard system, E = [] can be used. 
% According to the options specified by the combination of the input 
% parameters JOB and SYSTYPE, the resulting system representation of SYSR 
% has the following additional properties:
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
% The order reduction is achieved in five stages, by successively 
% eliminating the uncontrollable, unobservable and/or simple infinite 
% eigenvalues (non-dynamic modes) of the pencil A-lambda*E, as follows:
%    Stage 1 - eliminate uncontrollable finite eigenvalues 
%    Stage 2 - eliminate uncontrollable nonzero finite eigenvalues and
%              uncontrollable infinite eigenvalues
%    Stage 3 - eliminate unobservable finite eigenvalues
%    Stage 4 - eliminate unobservable nonzero finite eigenvalues and
%              unobservable infinite eigenvalues 
%    Stage 5 - eliminate simple infinite eigenvalues
% REDINFO is an eight-dimensional vector, which contains additional 
% informaton on the achieved order reductions at each stage, on the 
% structures of the resulting matrices Ar and Er, and on the rank of E:
% REDINFO(k)  > 0 - for k = 1..5, the reduction achieved in Stage k 
% REDINFO(k) <= 0 - for k = 1..5, no reduction achieved in Stage k 
% REDINFO(6) >= 0 - number of nonzero subdiagonals of Ar 
% REDINFO(7) >= 0 - number of nonzero subdiagonals of Er 
% REDINFO(8) >= 0 - rank of E 
% REDINFO(8)  < 0 - rank of E not computed
%
%



