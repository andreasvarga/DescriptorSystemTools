function [poles,info] = gpole(sys,tol,offset)
%GPOLE  Poles of a LTI descriptor system.
%   [POLES,INFO] = GPOLE(SYS,TOL,OFFSET) computes for the system 
%   SYS = (A-lambda*E,B,C,D), the complex vector POLES, representing 
%   the zeros of the linear pencil A-lambda*E. These are the poles of the 
%   transfer-function matrix of SYS, if A-lambda*E is regular and the 
%   descriptor system realization (A-lambda*E,B,C,D) is irreducible. 
%   If the pencil A-lambda*E is singular, POLES also contains NaN elements,
%   whose number is the rank deficiency of the pencil  A-lambda*E.
%   TOL is a tolerance used for rank determinations 
%   (Default: internally computed if not specified or 0).
%   OFFSET is the stability boundary offset, to be used to assess the  
%   finite eigenvalues which belong to the boundary of the stability domain 
%   as follows: in the continuous-time case, these are the finite 
%   eigenvalues having real parts in the interval [-OFFSET, OFFSET], while 
%   in the discrete-time case, these are the finite eigenvalues having moduli 
%   in the interval [1-OFFSET,1+OFFSET]. (Default: OFFSET = sqrt(eps)).
%
%   INFO is a structure containing additional information, as follows: 
%   INFO.nfev is the number of finite eigenvalues of the pencil A-lambda*E,
%            and also the number of finite poles of SYS.
%   INFO.niev is the number of infinite eigenvalues of the pencil A-lambda*E. 
%   INFO.nisev is the number of infinite simple eigenvalues of the 
%           pencil A-lambda*E. 
%   INFO.nip is the number of infinite poles of the system SYS. 
%   INFO.nfsev is the number of finite stable eigenvalues.
%   INFO.nfsbev is the number of finite eigenvalues on the boundary of the 
%           stability domain.
%   INFO.nfuev is the number of finite unstable eigenvalues.
%   INFO.nhev is the number of hidden eigenvalues 
%           (can be nonzero only if the pencil A-lambda*E is singular).  
%   INFO.nrank is the normal rank of the pencil A-lambda*E. 
%   INFO.miev is an integer vector, which contains the multiplicities 
%           of the infinite eigenvalues of the pencil A-lambda*E 
%           (also the dimensions of the elementary infinite blocks in the
%           Kronecker form of A-lambda*E).
%   INFO.mip is an integer vector, which contains the information on the  
%           multiplicities of the infinite zeros of A-lambda*E as follows: 
%           A-lambda*E has INFO.mip(i) infinite zeros of multiplicity i. 
%           INFO.mip is empty if A-lambda*E has no infinite zeros.
%   INFO.kr is an integer vector, which contains the right Kronecker 
%           indices of the pencil A-lambda*E. For a regular pencil, this
%           vector should be empty. 
%   INFO.kl is an integer vector, which contains the left Kronecker 
%           indices of the pencil A-lambda*E. For a regular pencil, this
%           vector should be empty. 
%   INFO.regular is set to:
%           true,  if the pencil A-lambda*E is regular, or  
%           false, if the pencil A-lambda*E is singular.
%   INFO.proper is 
%           true, if the pencil A-lambda*E is regular and all its infinite 
%                 eigenvalues are simple (has only non-dynamic modes), or 
%           false, if the pencil A-lambda*E is singular or it is regular, 
%                  but has non-simple infinite eigenvalues.
%   INFO.stable is 
%           true, if the pencil A-lambda*E is regular and has only stable  
%                 finite eigenvalues and all its infinite eigenvalues are
%                 simple (has only non-dynamic modes), or 
%           false, if the pencil A-lambda*E is singular or it is regular, 
%                  but has unstable finite eigenvalues or non-simple 
%                  infinite eigenvalues.
%
%   Remark: The multiplicities of non-simple infinite eigenvalues of the 
%           pencil A-lambda E are in excess with one to the multiplicities  
%           of the infinite zeros. 
%   
%    See also  GZERO, GKLF.


%  Author:    A. Varga, 22.12.2015.
%  Revisions: A. Varga, 15.12.2016, 23.04.2017, 12.09.2018.
%
%  References:
%
%  [1] Misra P., Van Dooren, P., Varga, A.:
%      Computation of structural invariants of generalized state space systems. 
%      Automatica, vol. 30, pp. 1921-1936, 1994.
%
%  [2] Svaricek, F.
%      Computation of the structural invariants of linear
%      multivariable systems with an extended version of
%      the program ZEROS.
%      System & Control Letters, vol. 6, pp. 261-266, 1985.

narginchk(1,3)
nargoutchk(0,2)

if ~isa(sys,'ss')
    error('The input system SYS must be an SS object')
end

if nargin == 1
  tol = 0;
  offset = sqrt(eps);
else
  validateattributes(tol, {'double'}, {'real','scalar', '>=', 0,'<',1},'','TOL',2) 
  if nargin == 2
     offset = sqrt(eps);
  else
     validateattributes(offset, {'double'},{'real','scalar','>',0,'<',1},'','OFFSET',3) 
  end
end

[a,~,~,~,e,Ts] = dssdata(sys);
discr = (Ts ~= 0);
n = size(a,1);

if isequal(e,eye(size(a,1)))
    % standard system
    poles = eig(a);
    if nargout > 1
       if discr
          temp = abs(poles);
          nfsev = sum(temp < 1-offset);
          nfsbev = sum(abs(temp-1) < offset);
          stable = all(temp < 1-offset);
       else
          temp = real(poles);
          nfsev = sum(temp < -offset);
          nfsbev = sum(abs(temp) < offset);
          stable = all(temp < -offset);
       end
       nfuev = n-nfsev-nfsbev;
       info = struct('nfev',n,'niev',0,'nisev',0,'nip',0,'nfsev',nfsev,...
                     'nfsbev',nfsbev,'nfuev',nfuev,'nhev',0,'nrank',n,...
                     'miev',[],'mip',[],'kr',[],'kl',[],...
                     'regular',true,'proper',true,'stable',stable); 
    end
else
    % descriptor system
    [pf,ni,mip,kr,miev,kl] = sl_gzero(a,e,[],[],[],tol); 
    nip = ni(1); % number of inifinite poles
    nr = ni(2);  % normal rank of the pole pencil
    poles = [pf;Inf*ones(nip,1); zeros(n-nr,1)*NaN];
    regular = (n == nr);
    if nargout > 1
       nfev = length(pf);        % number of finite eigenvalues/poles
       niev = sum(miev);         % number of infinite eigenvalues
       nhev = sum(kr)+sum(kl);  % number of hidden eigenvalues
       % nr = nfev + niev + nfspev must be fulfilled 
       nisev = sum(miev == 1);   % number of simple infinite eigenvalues
       proper = regular && (niev == nisev);
       if discr
          temp = abs(pf);
          nfsev = sum(temp < 1-offset);
          nfsbev = sum(abs(temp-1) < offset);
          stable = proper && all(temp < 1-offset);
       else
          temp = real(pf);
          nfsev = sum(temp < -offset);
          nfsbev = sum(abs(temp) < offset);
          stable = proper && all(temp < -offset);
       end
       nfuev = nfev-nfsev-nfsbev;
       info = struct('nfev',nfev,'niev',niev,'nisev',nisev,'nip',nip,...
                     'nfsev', nfsev,'nfsbev',nfsbev,'nfuev',nfuev,...
                     'nhev',nhev,'nrank',nr,'miev',miev','mip',mip',...
                     'kr',kr','kl',kl','regular',regular,'proper',...
                     proper,'stable',stable);  
    end
    if nargout < 2 && ~regular
       disp('Warning: The pole pencil of the system is singular')
    end    
end

% end GPOLE
end
