function [sysn,sysm] = gnrcf(sys,options)
%GNRCF  Generalized normalized right coprime factorization.
%    [SYSN,SYSM] = GNRCF(SYS,OPTIONS) calculates for the transfer 
%    function matrix
%                                       -1
%               G(lambda) =  C(lambda*E-A) B + D
%
%    of a given system SYS = (A-lambda*E,B,C,D) a normalized
%    right coprime factorization
%                                              -1
%               G(lambda) = N(lambda)*M(lambda)
% 
%    where N(lambda) and M(lambda) are stable transfer functions 
%    such that [N(lambda);M(lambda)] is inner.
%
%    The OPTIONS structure allows to specify various user options, 
%    as follows:
%    OPTIONS.tol      - tolerance for rank determinations
%                       (Default: internally computed)
%    OPTIONS.ss       - option to compute standard state-space 
%                       realizations of the factors:
%                       true  - standard state-space realizations, 
%                       false - descriptor system realizations (default).  
%    OPTIONS.balance  - specifies the balancing option for the Riccati 
%                       equation solvers: 
%                       true  - use balancing (Default);
%                       false - disable balancing.  
%
%    See also GNLCF.

%  Copyright 2018 A. Varga 
%  Author: A. Varga, 22.08.2018.
%  Revision(s): 
%
%  References:
%  [1]  A. Varga. 
%       A note on computing range space bases of rational matrices. 
%       2017. https://arxiv.org/abs/1707.00489.
%  [2] Oara, C. 
%      Constructive solutions to spectral and inner–outer factorizations 
%      with respect to the disk. Automatica, 41:1855–1866, 2005.

narginchk(1,2)
nargoutchk(0,2)

if ~isa(sys,'ss')
   error('The input system SYS must be an SS object')
end

if nargin < 2
   options = struct('tol',0);
else
   if isa(options,'double')
      % support for old syntax (undocummented)
      validateattributes(options, {'double'},{'real','scalar','>=',0},'','TOL') 
      options = struct('tol',options);
   else    
      validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
   end
end

% decode options

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

% standard state-space realization option
if isfield(options,'ss')
   standard = options.ss;
   validateattributes(standard, {'logical'},{'binary'},'','OPTIONS.ss')
else
   standard = false;
end

% balancing option for Riccati equation solvers
if isfield(options,'balance')
   balopt = options.balance;
   validateattributes(balopt, {'logical'},{'binary'},'','OPTIONS.balance')
else
   balance = true;
end

% compute the inner range of [G(lambda);I] 
[p,m]=size(sys);
sysr = grange([sys;eye(m)],struct('tol',tol,'inner',true,'balance',balance));

if standard
    % compute a standard state-space realization
    sysr = gss2ss(sysr,tol);
end

% extract the factors
sysn = sysr(1:p,1:m); sysm = sysr(p+1:end,1:m);

% end GNRCF
