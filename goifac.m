function [sysi,syso,info] = goifac(sys,options)
%GOIFAC Co-outer-inner/RQ-like factorization.
% [SYSI,SYSO,INFO] = GOIFAC(SYS,OPTIONS) calculates, for a descriptor 
% system SYS = (A-lambda*E,B,C,D) with the transfer-function matrix 
% G(lambda), the square inner factor SYSI = (Ai-lambda*Ei,Bi,Ci,Di) with 
% the transfer function matrix Gi(lambda) and the minimum-phase 
% quasi-co-outer or the full column rank factor SYSO = (A-lambda*E,B,Co,Do) 
% with the transfer function matrix Go(lambda) such that 
%
%                 G(lambda) = Go(lambda)*Gi(1:r,:)(lambda) ,    (*)
%
% where r is the normal rank of G(lambda). The resulting proper and stable 
% inner factor satisfies Gi(lambda)*Gi'(lambda) = I. The resulting factor 
% Go(lambda) has full column rank r. Depending on the selected options,
% Go(lambda) is either minimum phase, excepting possibly zeros on the 
% boundary of the appropriate stability domain, or contains all zeros 
% of G(lambda), in which case (*) is the extended RQ-like factorization of
% G(lambda).
% Note: SYSO may generally contain free inner factors which can be 
%       eliminated by removing the finite unobservable eigenvalues. 
%
% The OPTIONS structure allows to specify various user options, 
% as follows:
% OPTIONS.tol      - tolerance for rank determinations and observability 
%                    tests (Default: internally computed)
% OPTIONS.minphase - option to compute a minimum phase quasi-co-outer factor:
%                    true  - compute a minimum phase quasi-co-outer factor, 
%                            with all zeros stable, excepting possibly 
%                            zeros on the boundary of the appropriate 
%                            stability domain (Default);
%                    false - compute a full column rank factor, which 
%                            includes all zeros of G(lambda).  
% OPTIONS.offset   - stability boundary offset OFFSET, to specify the 
%                    marginally stable finite zeros as follows: 
%                    in the continuous-time case are the finite zeros 
%                    having real parts in the interval [-OFFSET, OFFSET], 
%                    while in the discrete-time case having moduli in the 
%                    interval [1-OFFSET,1+OFFSET] 
%                    (Default: OFFSET = sqrt(eps)).
%                    
% OPTIONS.balance  - specifies the balancing option for the Riccati 
%                    equation solvers: 
%                    true  - use balancing (Default);
%                    false - disable balancing.  
%   
% INFO is a structure containing additional information, as follows: 
% INFO.rankG  - normal rank r of the transfer function matrix G(lambda);
% INFO.nfuz   - number of finite unstable zeros of SYSO; these are the 
%               finite zeros of SYS lying on the boundary of the stability 
%               region;
% INFO.niuz   - number of infinite zeros of SYSO; these are the infinite  
%               zeros of SYS in the continuous-time case and 0 in the 
%               discrete-time case
% INFO.ricrez - diagnosis flag, as provided provided by the generalized 
%               Riccati equation solvers care and dare; if non-negative, 
%               this value represents the Frobenius norm of relative 
%               residual of the Riccati equation, while a negative value 
%               indicates failure of solving the Riccati equation. 
%
%
% See also GIOFAC.

% Copyright 2015-2018 A. Varga
% Author:      A. Varga, 22.12.2015.
% Revision(s): A. Varga, 26.07.2017, 21.08.2018, 09.09.2018.
%
% References:
% [1] C. Oara and A. Varga.
%     Computation of the general inner-outer and spectral factorizations.
%     IEEE Trans. Autom. Control, vol. 45, pp. 2307--2325, 2000.
% [2] C. Oara.
%     Constructive solutions to spectral and inner–outer factorizations 
%     respect to the disk. Automatica,  41, pp. 1855–1866, 2005. 

narginchk(1,2)
nargoutchk(0,3)

if ~isa(sys,'ss')
   error('The input system SYS must be a state space system')
end

if nargin < 2
   options = struct('tol',0,'minphase',true);
else
   if isa(options,'double')
      % support for old syntax (undocumented)
      validateattributes(options, {'double'},{'real','scalar','>=',0},'','TOL') 
      options = struct('tol',options,'minphase',true);
   else    
      validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
   end
end


[sysi,syso,info] = giofac(sys.',options);
sysi = sysi.';
syso = syso.';

% end of GOIFAC
end