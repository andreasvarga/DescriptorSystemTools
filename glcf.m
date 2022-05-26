function [sysn,sysm] = glcf(sys,options)
%GLCF   Left coprime factorization with proper and stable factors.
%       [SYSN,SYSM] = GLCF(SYS,OPTIONS) calculates for the LTI descriptor
%       system SYS = (A-lambda E,B,C,D) with the transfer function matrix 
%                                        -1
%               G(lambda) =  C(lambda*E-A) B + D
%
%       a left coprime factorization
%                                    -1          
%               G(lambda) = M(lambda)  *N(lambda)   ,
% 
%       where N(lambda) and M(lambda) are proper and stable transfer 
%       function matrices of two systems SYSN and SYSM, respectively, 
%       with their eigenvalues having a prescribed stability degree. 
%       The OPTIONS structure allows to specify various user options, 
%       as follows:
%       OPTIONS.tol    - specifies the tolerance for the singular values 
%                        based rank determination of E 
%                        (Default: tol = prod(size(E)) * eps(norm(E,1)))
%       OPTIONS.tolmin - specifies the tolerance for the singular values 
%                        based observability tests 
%                        (Default: tol = prod(size(C)) * eps(norm(C',1)))
%       OPTIONS.smarg  - sets the stability margin
%                        (Default: -sqrt(eps) for a continuous-time system 
%                         and 1-sqrt(eps) for a discrete-time system) 
%       OPTIONS.sdeg   - specifies a prescribed stability degree for the 
%                        eigenvalues of the factors
%                        (Default:  -0.05 for a continuous-time system and
%                                    0.95  for a discrete-time system) 
%       OPTIONS.poles  - specifies a complex conjugated set of desired 
%                        poles to be assigned for the factors (Default: []) 
%       OPTIONS.mindeg - minimum degree option for the denominator SYSM
%                        true  - determine minimum degree denominator
%                        false - both factors have the same order (default) 
%       OPTIONS.mininf - specifies option for removal of simple infinite
%                        eigenvalues of SYSN and SYSM
%                        true  - remove simple infinite eigenvalues
%                        false - keep simple infinite eigenvalues (default) 
%
%  See also GRCF.

%  Author:    A. Varga 23.11.2015.
%  Revisions: A. Varga 05.12.2016, 06.06.2017.
%
%  Method:  The Procedure GRCF from [1] is applied to the dual system. 
%  The underlying algorithm represents an extension of the recursive 
%  factorization approach of [2] to cope with infinite poles. 
%  All infinite poles are assigned to finite real values. 
%  If OPTIONS.poles is empty or does not contain a sufficient 
%  number of real values, then a part or all of infinite poles are 
%  assigned to the value specified by OPTIONS.sdeg.  
%
%  References:
%  [1] Varga A.
%      On recursive computation of coprime factorizations of rational 
%      matrices. arXiv:1703.07307, https://arxiv.org/abs/1703.07307, 2017.
%  [2] A. Varga. 
%      Computation of coprime factorizations of rational matrices.
%      Linear Algebra and Its Applications, vol. 271, pp.88-115, 1998.

narginchk(1,2)
nargoutchk(0,2)

% check for state-space system
if ~isa(sys,'ss')
   error('The input system SYS must be an SS object')
end

if nargin == 1
   options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

% build the dual system, by trying to preserve upper-triangular shapes
if nargout == 1
   sysn = grcf(xperm(sys,order(sys):-1:1).',options);
else
   [sysn,sysm] = grcf(xperm(sys,order(sys):-1:1).',options);
end
   
% build dual factors, by preserving upper-triangular shapes
sysn = xperm(sysn,order(sysn):-1:1).';
if nargout == 2
   sysm = xperm(sysm,order(sysm):-1:1).';
end

% end GLCF
end
