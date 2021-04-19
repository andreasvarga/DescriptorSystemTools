function [varargout] = glsol(sysg,sysf,options)
%GLSOL  Solution of the linear rational matrix equation X*G = F.
%   [SYSX,INFO] = GLSOL(SYSG,SYSF,OPTIONS) returns a state-space  
%   realization SYSX of the solution X(lambda) of the linear rational 
%   equation
%
%      X(lambda)*G(lambda) = F(lambda)                  (1)
%
%   given the state-space realizations SYSG and SYSF of G(lambda) and 
%   F(lambda). The OPTIONS structure (optional) specifies various user 
%   options, as follows:
%   OPTIONS.tol     - specifies the relative tolerance for rank 
%                     determinations. (Default: internally computed)
%   OPTIONS.sdeg    - specifies a prescribed stability degree for the free
%                     eigenvalues of SYSX 
%                     (Default: [], i.e., no stabilization performed)
%   OPTIONS.poles   - specifies a complex conjugated set of desired poles
%                     to be assigned as free eigenvalues of SYSX 
%                     (Default: []) 
%   OPTIONS.mindeg  - minimum degree solution 
%                     true  - determine minimum order solution
%                     false - determine a particular solution which has 
%                             possibly non-minimal order (default) 
%   The resulting INFO structure contains additional information:
%   INFO.nrank is the normal rank of G(lambda); 
%   INFO.rdeg is a vector which contains the relative row degrees of 
%       X(lambda) (i.e., the numbers of integrators/delays needed to make 
%       each row of X(lambda) proper). 
%   INFO.tcond is the maximum of the Frobenius-norm condition numbers 
%       of the employed transformation matrices.
%   INFO.fnorm is the Frobenius-norm of the employed 
%       state-feedback/feedforward used for dynamic cover computation if 
%       OPTIONS.mindeg = true, or for stabilization of the free eigenvalues
%       if OPTION.sdeg is non-empty;  
%   INFO.nl is the number of freely assignable poles of the solution X(lambda). 
%   Note: Large values of INFO.tcond and/or INFO.fnorm indicate possible
%         loss of numerical stability. 
%
%   [SYSX,INFO] = GLSOL(SYSGF,MF,OPTIONS) uses the compound realization 
%   SYSGF of [G(lambda); F(lambda)], where F(lambda) has MF rows.
%
%   [SYSX,INFO,SYSGEN] = GLSOL(SYSG,SYSF,OPTIONS) determines additionally
%   a generator of all solution SYSGEN = [ SYS0; SYSN ], where SYS0 is a 
%   particular solution of (1) and SYSN is a proper left nullspace basis
%   of G(lambda) satisfying SYSN(lambda)*G(lambda) = 0. All solutions 
%   of (1) can be generated as
%               SYSX = SYS0 + SYSY*SYSN, 
%   where SYSY is an arbitrary system with appropriate dimensions. 
%   The descriptor system realization of SYSGEN is usually not minimal 
%   (unobservable and/or non-dynamic modes present) and has the form 
%   (Ag-s*Eg,Bg,[C0;CN],[D0;DN]), where
%                   ( Ai-lambda*Ei     *           *       )       ( *  )
%    Ag-lambda*Eg = (    0        Af-lambda*Ef     *       ), Bg = ( *  ),
%                   (    0             0      Al-lambda*El )       ( Bl )
%
%      ( C0 )  = (   C1       C2      C3    ) 
%      ( CN )    (   0        0       Cl    )
%   with Ai, Ef and El invertible and upper triangular, Ei nillpotent and
%   upper triangular, and DN full column rank. A minimal order descriptor 
%   system realization of the proper basis SYSN is (Al-lambda*El,Bl,Cl,DN).
%   Cl and DN have P-INFO.nrank rows, where P is the number of 
%   outputs of SYSG.
%
%   The INFO structure contains further information on the dimensions of
%   the square diagonal blocks of the pencil Ag-lambda*Eg, as follows:
%   INFO.ninf is the dimension of Ai-lambda*Ei (the number of infinite
%      eigenvalues);
%   INFO.nf is the dimension of Af-lambda*Ef (the number of finite eigenvalues);
%   INFO.nl is the dimension of Al-lambda*El (also the column dimension 
%      of Cl).
%
%   See also GRSOL.

%  Copyright 2016-2018 A. Varga 
%  Author:      A. Varga, 14.01.2016.
%  Revision(s): A. Varga, 02.11.2016, 07.05.2017, 21.08.2018, 27.07.2019.
%
%  Method:  The method of [1] to solve linear rational systems is used.
%
%  References:
%  [1] A. Varga, "Computation of least order solutions of linear 
%      rational equations", Proc. MTNS'04, Leuven, Belgium, 2004.

narginchk(2,3)
nargoutchk(0,3)

if ~isa(sysg,'ss')  
   error('The input system SYSG must be an SS object')
end  

if ~isa(sysf,'ss')
   if ~isa(sysf,'double') 
      error('SYSF must be an SS object or a positive integer')
   else
      validateattributes(sysf,{'double'},{'integer','scalar','>=',0,'<=',size(sysg,1)},'','MF')
   end
else
   if size(sysg,2) ~= size(sysf,2) 
      error('The systems SYSG and SYSF must have the same number of inputs')
   end    
end  

if nargin < 3
   options = struct('tol',0);
end

if nargout > 2
   [sysx,info,sysgen] = grsol(sysg.',sysf.',options); 
   sysgen = xperm(sysgen,order(sysgen):-1:1).';
else
   [sysx,info] = grsol(sysg.',sysf.',options); 
end    
sysx = xperm(sysx,order(sysx):-1:1).';

% set INFO
info.nl = info.nr;
info = rmfield(info,'nr');
info = orderfields(info,[1 2 3 4 6 5 7]);

if nargout > 2
   varargout = {sysx,info,sysgen};
else
   varargout = {sysx,info}; 
end    

% end GLSOL
end
