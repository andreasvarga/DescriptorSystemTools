function [z,info] = gzero(sys,tol,offset)
%GZERO  Zeros of a LTI descriptor system. 
%   [Z,INFO] = GZERO(SYS,TOL,OFFSET) computes for the system 
%   SYS = (A-lambda*E,B,C,D), the complex vector Z, representing 
%   the (invariant) zeros of the of the system pencil
%        S(lambda) = ( A - lambda E  B )
%                    (     C         D ). 
%   These are the transmission zeros of the LTI system SYS, if the 
%   descriptor realization (A-lambda E,B,C,D) is irreducible. 
%   TOL is a tolerance used for rank determinations 
%   (Default: internally computed if not specified or 0).
%   OFFSET is the stability boundary offset, to be used to assess the  
%   finite zeros which belong to the boundary of the stability domain 
%   as follows: in the continuous-time case, these are the finite 
%   zeros having real parts in the interval [-OFFSET, OFFSET], while 
%   in the discrete-time case, these are the finite zeros having moduli 
%   in the interval [1-OFFSET,1+OFFSET]. (Default: OFFSET = sqrt(eps)).
%
%   INFO is a structure containing additional information, as follows: 
%   INFO.nfz is the number of finite eigenvalues of the pencil S(lambda),
%            and also the number of finite zeros of SYS.
%   INFO.niev is the number of infinite eigenvalues of the pencil S(lambda). 
%   INFO.nisev is the number of infinite simple eigenvalues of the 
%           pencil S(lambda). 
%   INFO.niz is the number of infinite zeros of the system SYS.
%   INFO.nfsz is the number of finite stable zeros 
%   INFO.nfsbz is the number of finite zeros on the boundary of the 
%           stability domain
%   INFO.nfuz is the number of finite unstable zeros 
%   INFO.nrank is the normal rank of the pencil S(lambda). 
%   INFO.miev is an integer vector, which contains the multiplicities 
%           of the infinite eigenvalues of the pencil S(lambda) 
%           (also the dimensions of the elementary infinite blocks in the
%           Kronecker form of S(lambda))
%   INFO.miz is an integer vector, which contains the information on the  
%           multiplicities of the infinite zeros of S(lambda) as follows: 
%           S(lambda) has INFO.miz(i) infinite zeros of multiplicity i. 
%           INFO.miz is empty if A-lambda*E has no infinite zeros.
%   INFO.kr is an integer vector, which contains the right Kronecker 
%           indices of the pencil S(lambda). For a regular pencil, this
%           vector should be empty. 
%   INFO.kl is an integer vector, which contains the left Kronecker 
%           indices of the pencil S(lambda). For a regular pencil, this
%           vector should be empty. 
%   INFO.minphase is 
%           true, if the pencil S(lambda) has only stable finite zeros, or 
%           false, if the pencil S(lambda) has unstable finite zeros or  
%                  infinite zeros.
%
%   Note: The finite zeros and the finite eigenvalues of the pencil
%   S(lambda) are the same, but the multiplicities of infinite eigenvalues 
%   are in excess with one to the multiplicities of infinite zeros. 
%
%   See also  GPOLE, GKLF.

%  Author:    A. Varga, 11.01.2016.
%  Revisions: A. Varga, 02.10.2016, 31.10.2016, 13.09.2018.
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
%
%  [3] Emami-Naeini, A. and Van Dooren, P.
%      Computation of zeros of linear multivariable systems.
%      Automatica, vol. 18, pp. 415-430, 1982.

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

[a,b,c,d,e] = dssdata(sys);
[zf,ni,miz,kr,miev,kl] = sl_gzero(a,e,b,c,d,tol); 
niz = ni(1);  % number of inifinite zeros
nr  = ni(2);  % normal rank of the system matrix pencil
z = [zf;Inf*ones(ni(1),1)];

if nargout > 1
   nfz = length(zf);         % number of finite eigenvalues/zeros
   niev = sum(miev);         % number of infinite eigenvalues
   nisev = sum(miev == 1);   % number of simple infinite eigenvalues
   if sys.Ts ~= 0
      temp = abs(zf);
      nfsz = sum(temp < 1-offset);       % number of finite stable zeros
      nfsbz = sum(abs(temp-1) < offset); % number of finite zeros on the unit circle
      stable = all(temp < 1-offset);     
   else
      temp = real(zf);
      nfsz = sum(temp < -offset);        % number of finite stable zeros
      nfsbz = sum(abs(temp) < offset);   % number of finite zeros on the imaginary axis
      stable = (niev == nisev) && all(temp < -offset);
   end
   nfuz = nfz-nfsz-nfsbz;
   info = struct('nfz',nfz,'niev',niev,'nisev',nisev,'niz',niz,...
                 'nfsz',nfsz,'nfsbz',nfsbz,'nfuz',nfuz,'nrank',nr,...
                 'miev',miev','miz',miz','kr',kr','kl',kl',...
                 'minphase',stable);  
end


% end GZERO
