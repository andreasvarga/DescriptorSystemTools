function [sysn,sysm] = glcfid(sys,options)
%GLCFID Left coprime factorization with inner denominator.
%       [SYSN,SYSM] = GLCFID(SYS,OPTIONS) calculates for the LTI descriptor
%       system SYS = (A-lambda E,B,C,D) with the transfer function matrix 
%                                        -1
%               G(lambda) =  C(lambda*E-A) B + D
%
%       a left coprime factorization 
%                                              -1
%               G(lambda) = N(lambda)*M(lambda)
% 
%       where N(lambda) and M(lambda) are proper and stable transfer 
%       function matrices of two systems SYSN and SYSM, respectively, 
%       and M(lambda) is inner.  
%       The OPTIONS structure allows to specify various user options, 
%       as follows:
%       OPTIONS.tol    - specifies the tolerance for the singular values 
%                        based rank determination of E 
%                        (Default: tol = prod(size(E)) * eps(norm(E,1)))
%       OPTIONS.tolmin - specifies the tolerance for the singular values 
%                        based observability tests 
%                        (Default: tol = prod(size(C)) * eps(norm(C',1)))
%       OPTIONS.mindeg - minimum degree option of denominator SYSM
%                        true  - determine minimum degree denominator
%                        false - both factors have the same order (default) 
%       OPTIONS.mininf - specifies option for removal of simple infinite
%                        eigenvalues of SYSN and SYSM
%                        true  - remove simple infinite eigenvalues
%                        false - keep simple infinite eigenvalues (default) 
%
%       See also GRCFID.

%  Author:    A. Varga 21.01.2016.
%  Revisions: A. Varga 05.12.2016, 07.06.2017. 
%
%  Method:  An extension of the recursive factorization approach of [1]  
%           is used (see [2] for details). 
%
%  References:
%  [1] A. Varga. 
%      Computation of coprime factorizations of rational matrices.
%      Linear Algebra and Its Applications, vol. 271, pp.88-115, 1998.
%  [2] Varga A.
%      On recursive computation of coprime factorizations of rational 
%      matrices. arXiv:1703.07307, https://arxiv.org/abs/1703.07307, 2017.

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

if nargout == 1
   sysn = grcfid(sys.',options);
else
   [sysn,sysm] = grcfid(sys.',options);
end
   
sysn = xperm(sysn,order(sysn):-1:1).';
if nargout == 2
   sysm = xperm(sysm,order(sysm):-1:1).';
end

% end GLCFID
end
