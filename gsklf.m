function [At,Et,dimsc,Q,Z] = gsklf(sys,tol,varargin)
%GSKLF Special Kronecker-like form of a system matrix pencil. 
% [At,Et,DIMSC,Q,Z] = GSKLF(SYS,TOL) computes, for a given LTI system
% with an infinite controllable descriptor system representation 
% (A-lambda*E,B,C,D), a special Kronecker-like form (SKLF)
% of the system matrix pencil, such that
%                         ( A-lambda*E B )
%       At-lambda*Et = Q'*(--------------)*Z 
%                         (     C      D )
%
% where Q and Z are orthogonal transformation matrices and the reduced 
% pencil At-lambda*E has the special structure
%
%
%                 ( Arg-lambda*Erg     *            *   *   )
%                 (   0            Abl-lambda*Ebl  Bbl  *   )
%  At-lambda*Et = (   0                0            0   Bn  ) .    (1)
%                 (-----------------------------------------)
%                 (   0               Cbl          Dbl  *   )
%
% TOL is an optional tolerance used for rank tests 
% (Default: internally computed, if TOL = 0 or not provided). 
%
% The reduced pencil At-lambda*Et contains the following blocks:
%   [ Arg-lambda*Erg ], is a full row rank pencil, which contains the  
%      right Kronecker structure and all or a part of the system zeros; 
%      DIMSC(1) contains the column dimension of this pencil; 
%   Abl-lambda*Ebl is a square pencil of order DIMSC(2), with Ebl
%      invertible;
%   Bbl is a DIMSC(2) x DIMS(3) matrix, where DIMSC(3) is the normal rank 
%      of the transfer function matrix of SYS; 
%   Bn is an invertible matrix of order DIMSC(4).
%   If non-negative, DIMSC(5) contains the number of marginally stable 
%   zeros of SYS on the boundary of the stability region (i.e., imaginary
%   axis in continuous-time case, unit circle in the origine in 
%   discrete-time case). 
%   If non-negative, DIMSC(6) contains, in the continuous-time case,
%   the number of infinite zeros of SYS and is set zero in the 
%   discrete-time case.
%
%   The subpencil
%                           ( Abl-lambda*Ebl  Bbl )
%                           (      Cbl        Ddl )
%
%   is full column rank and the pair (Abl-lambda*Ebl,Bbl) is stabilizable
%   if the pair (A-lambda*E,B) is finite stabilizable. 
%
% [At,Et,DIMSC,Q,Z] = GSKLF(SYS,TOL,'none') is the same as above.
%
% [At,Et,DIMSC,Q,Z] = GSKLF(SYS,TOL,'unstable') computes the SKLF (1) 
% such that Abl-lambda*Ebl contains all unstable zeros of SYS, while 
% Arg-lambda*Erg includes the rest of (stable) zeros. 
% For a continuous-time system, the selected unstable zeros either are 
% infinite or have real parts greater than sqrt(eps), while for a  
% discrete-time system the selected unstable zeros have moduli greater 
% than 1+sqrt(eps) (i.e., including also the infinite zeros). 
% DIMSC(5) contains the number of marginally stable zeros of SYS.
% DIMSC(6) contains the number of infinite zeros of SYS in the 
% continuous-time case and zero in the discrete-time case. 
%
% [At,Et,DIMSC,Q,Z] = GSKLF(SYS,TOL,'s-unstable') computes the SKLF (1) 
% such that Abl-lambda*Ebl contains all strictly unstable zeros of SYS,  
% while Arg-lambda*Erg includes the rest of strictly stable zeros. 
% For a continuous-time system, the selected strictly unstable zeros  
% are infinite or have real parts greater than -sqrt(eps), while for a  
% discrete-time system the selected strictly unstable zeros have moduli  
% greater than 1-sqrt(eps) (i.e., including also the infinite zeros). 
% DIMSC(5) contains the number of marginally stable zeros of SYS.
% DIMSC(6) contains the number of infinite zeros of SYS in the 
% continuous-time case and zero in the discrete-time case. 
% 
% [At,Et,DIMSC,Q,Z] = GSKLF(SYS,TOL,'stable') computes the SKLF (1) such
% that Abl-lambda*Ebl contains all (strictly) stable zeros of SYS, while 
% Arg-lambda*Erg includes the rest of (unstable) zeros. 
% For a continuous-time system, the selected stable zeros have real 
% parts less than -sqrt(eps), while for a discrete-time system the
% selected stable zeros have moduli less than 1-sqrt(eps). 
% DIMSC(5) contains the number of marginally stable zeros of SYS.
% DIMSC(6) contains the number of infinite zeros of SYS in the 
% continuous-time case and zero in the discrete-time case. 
% 
% [At,Et,DIMSC,Q,Z] = GSKLF(SYS,TOL,'all') computes the SKLF (1) such
% that Abl-lambda*Ebl contains all zeros of SYS. 
% 
% [At,Et,DIMSC,Q,Z] = GSKLF(SYS,TOL,'finite') computes the SKLF (1) such
% that Abl-lambda*Ebl contains all finite zeros of SYS, while 
% Arg-lambda*Erg includes the rest of (infinite) zeros. 
% 
% [At,Et,DIMSC,Q,Z] = GSKLF(SYS,TOL,'infinite') computes the SKLF (1) such
% that Abl-lambda*Ebl contains all infinite zeros of SYS, while 
% Arg-lambda*Erg includes the rest of (finite) zeros. 
%
% [At,Et,DIMSC,Q,Z] = GSKLF(SYS,TOL,...,OFFSET) uses OFFSET, the stability 
% boundary offset, to be used to assess the finite zeros which belong to 
% the boundary of the stability domain as follows: 
% in the continuous-time case, these are the finite zeros having real parts
% in the interval [-OFFSET, OFFSET], while in the discrete-time case, 
% these are the finite zeros having moduli in the interval 
% [1-OFFSET,1+OFFSET]. (Default: OFFSET = sqrt(eps)).
% 
% [At,Et,INFO,Q,Z] = GSKLF(SYS,TOL,...,'noQ') computes the SKLF (1)
% without internally accumulating the transformation matrix Q. 
% In this case, Q = [] (the empty matrix). 
%
% See also GKLF, SL_KLF.

% Copyright 2017-2020 A. Varga
% Author:      A. Varga, 01.07.2017.
% Revision(s): A. Varga, 09.11.2017, 21.08.2018, 15.09.2018, 2.11.2020. 
%
% Method: The reduction algorithm of [1] has been adapted to deal with
% several zero selection options. The computation of the involved
% Kronecker-like form is based on the algorithm of [2]. 
%
% References
% [1] Oara, C. 
%     Constructive solutions to spectral and inner–outer factorizations 
%     with respect to the disk. Automatica, 41:1855–1866, 2005.
% [2] Beelen, Th. and Van Dooren, P.
%     An improved algorithm for the computation of Kronecker's
%     canonical form of a singular pencil.
%     Linear Algebra and Applications, vol. 105, pp. 9-65, 1988.


narginchk(1,5)
nargoutchk(0,5)

if ~isa(sys,'ss')
   error('The input system SYS must be a state-space system object')
end
discr = (sys.Ts ~= 0);

if nargin < 2
   tol = 0; 
else
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','TOL') 
end

opts = {'none','unstable','all','infinite','finite','stable','s-unstable','noQ'};

withQ = nargout > 3; 
withZ = (nargout > 4); 
if nargin < 3
   job = 0; 
else
   % nargin >= 3 
   validateattributes(varargin{1}, {'char'},{'nonempty'},'','argument 3') 
   job = find(strcmp(varargin{1},opts),1)-1;
   if isempty(job) || (job > 6 && nargin > 3)
      error('No such zero selection option')
   else
      if job > 6 
         job = 0;
         withQ = false;
      end
      if nargin == 3
         offset = sqrt(eps);
      end
   end
   if nargin == 4
      % nargin = 4
      if ischar(varargin{2})
         if strcmp(varargin{2},'noQ')
            withQ = false;
            offset = sqrt(eps);
         else
            error('No such option')
         end
      else
         validateattributes(varargin{2}, {'double'},{'real','scalar','>',0,'<',1},'','OFFSET',4)
         offset = varargin{2};
      end
   elseif nargin > 4
      % nargin = 5
      validateattributes(varargin{2}, {'double'},{'real','scalar','>',0,'<',1},'','OFFSET',4)
      offset = varargin{2};
      validateattributes(varargin{3}, {'char'},{'nonempty'},'','argument 5') 
      if strcmp(varargin{3},'noQ')
         withQ = false;
      else
         error('No such option for input argument 5')
      end
   end
end

[a,b,c,d,e] = dssdata(sys); 
[n,m]=size(b); 
sstype = isequal(e,eye(n));

%  Step 0: Compute orthogonal Q and U0 such that 
%                                       
%                 ( A-s*E | B )        ( A11-s*E11  A12-s*E12  A13-s*E13 )
%    diag(Q',I) * (-------|---) * U0 = (    0          0         B2      )
%                 (   C   | D )        (    C1         C2        D1      )
%      
%           with E11(n1,n1) and B2(n-n1,n-n1) invertible.

nm = n+m;
if sstype
   AB=[a b]; CD = [c d]; EF = eye(n,nm); 
   if withZ, Z = eye(nm); end
   if withQ, Q = eye(n); end
   n1 = n; nsinf = 0; i1 = 1:n;  
else
   [q,se,z] = svd(e);
   tolr = tol; se = diag(se);
   if tol == 0
      tolr = n * max(se) * eps;
   end
   n1 = sum(se > tolr); nsinf = n-n1; i1 = 1:n1; 
   if nsinf
      i2 = n1+1:n; 
      if withQ, Q = q; end
      if withZ, Z = blkdiag(z,eye(m)); end
      AB = [q'*a*z q'*b]; CD = [c*z d];
      sb = svd(AB(i2,n+1:nm));
      if m < nsinf || sb(nsinf) <= max(n,m)*eps(norm(b,1))
         error(' The system is not infinite controllable')
      end
      EF = [ diag(se(i1)) zeros(n1,nsinf+m) ; zeros(nsinf,nm) ]; 
      [q2,ab2,z2] = svd(AB(i2,:));
      jp = [nsinf+1:nm 1:nsinf];
      z2 =z2(:,jp); 
      AB(i2,:) = ab2(:,jp); AB(i1,:) = AB(i1,:)*z2;
      EF(i1,:) = EF(i1,:)*z2;  CD = CD*z2;
      if withZ, Z = Z*z2; end
      if withQ, Q(:,i2) = Q(:,i2)*q2; end
   else
      AB=[a b]; CD = [c d]; EF = [e zeros(n,m)]; 
      if withZ, Z = eye(nm); end
      if withQ, Q = eye(n); end
   end
end

%  Step 1: Compress [C1 C2] to [C1 C2]*U1 = [0 D2] with D2(p,rcd) monic.
%
%          Compute  [ A11-s*E11  A12-s*E12 ]*U1 = [ A1-s*E1   B1-s*F1 ]

n1m = n1+m; j1 = 1:n1m;
[~,se,z]=svd(CD(:,j1)); 
if min(size(se)) > 1
   se = diag(se);
else
   se = se(:);
end 
if tol == 0
   tolr = n1m * max(se) * eps;
else
   tolr = max(se) * tol; 
end
rcd = sum(se > tolr);
jp = [ rcd+1:n1m 1:rcd ]; z = z(:,jp);
CD(:,j1) = CD(:,j1)*z; 
EF(i1,j1) = EF(i1,j1)*z; AB(i1,j1) = AB(i1,j1)*z;
if withZ, Z(:,j1) = Z(:,j1)*z; end

% Step 2: Compute the  Kronecker-like staircase form of A1-s*E1 such that
%
%                              (A11-s*E11    X         X      |    X   )
% Q2*(A1-s*E1 | B1-s*F1 )*Z2 = (    0     A22-s*E22    X      |    X   )
%                              (    0        0      A33-s*E33 | B3-s*F3)

n3 = n1m-rcd;
j1 = 1:n3; j2 = n3+1:n1m; j3 = n1m+1:nm; 

if n3
   if job == 0 || job == 3
       % compute the KLF with (Arf,Ainf,Al) diagonal blocks 
       [AB(i1,j1),EF(i1,j1),info1,q,z] = gklf(AB(i1,j1),EF(i1,j1),tol,'left'); 
   elseif job == 1
      if discr
         % compute the KLF with (Ar,Af,Ainf,Al) diagonal blocks
         [AB(i1,j1),EF(i1,j1),info1,q,z] = gklf(AB(i1,j1),EF(i1,j1),tol,'reverse'); 
      else
         % compute the KLF with (Ar,Ainf,Af,Al) diagonal blocks
         [AB(i1,j1),EF(i1,j1),info1,q,z] = gklf(AB(i1,j1),EF(i1,j1),tol); 
      end
   elseif job == 2 || job == 4 || job == 5
      % compute the KLF with (Ar,Ainf,Af,Al) diagonal blocks
      [AB(i1,j1),EF(i1,j1),info1,q,z] = gklf(AB(i1,j1),EF(i1,j1),tol); 
   elseif job == 6
       % compute the KLF with (Ar,Af,Ainf,Al) diagonal blocks
       [AB(i1,j1),EF(i1,j1),info1,q,z] = gklf(AB(i1,j1),EF(i1,j1),tol,'reverse'); 
   end
   AB(i1,[j2 j3]) = q'*AB(i1,[j2 j3]);
   EF(i1,[j2 j3]) = q'*EF(i1,[j2 j3]);
   if withQ, Q(:,i1) = Q(:,i1)*q; end
   if withZ, Z(:,j1) = Z(:,j1)*z; end
   nc = abs(sum(info1.mr)); mc = abs(sum(info1.nr));
   ni = sum(info1.minf); 
   nl = sum(info1.ml); nf = info1.mf;
else
   nc=0; mc = 0; ni = 0; nf = 0; nl = n1;  
end

nmszer = -1; nizer = -1;
if job == 0
   % no zeros included
   nus = 0; ns = nf+ni;
elseif job == 1 
   % include all unstable zeros excluding those on the boundary of the
   % stability region
   if nf == 0
      if discr
         nus = ni; ns = 0; nmszer = 0; nizer = 0;
      else
         nus = 0; ns = ni; nmszer = 0; nizer = ni;
      end
   else    
      if discr
         % the finite zeros are in the leading block
         izf = nc+1:nc+nf; jzf = mc+1:mc+nf;
      else
         % the finite zeros are in the trailing block
         izf = nc+ni+1:nc+ni+nf; jzf = mc+ni+1:mc+ni+nf;
      end
      [AB(izf,jzf),EF(izf,jzf),q,z] = qz(AB(izf,jzf),EF(izf,jzf),'real');
      ev = ordeig(AB(izf,jzf),EF(izf,jzf)); 
      if discr
         select = abs(ev) <= 1+offset;
         nizer = 0; 
         nmszer = length(find(abs(ev(select)) >= 1-offset));
      else
         select = (real(ev) <= offset);
         nizer = ni;
         nmszer = length(find(real(ev(select)) >= -offset));
      end
      ns = length(select(select==true)); 
      [AB(izf,jzf),EF(izf,jzf),q,z] = ordqz(AB(izf,jzf),EF(izf,jzf),q,z,select); 
   
      if discr
         % apply q to the trailing columns including infinite zeros
         j4 = mc+nf+1:nm;  
         izf1 = 1:nc;
         nus = nf-ns+ni;   % include infinite zeros among unstable zeros
      else
         % apply q to the trailing columns without infinite zeros
         j4 = mc+ni+nf+1:nm;  
         izf1 = 1:(nc+ni);
         nus = nf-ns;      % include only the finite unstable zeros
         ns = ns+ni;       % include infinite zeros among stable zeros
      end
      AB(izf,j4) = q*AB(izf,j4); EF(izf,j4) = q*EF(izf,j4);
      AB(izf1,jzf) = AB(izf1,jzf)*z; EF(izf1,jzf) = EF(izf1,jzf)*z;
      if withQ, Q(:,izf) = Q(:,izf)*q'; end
      if withZ, Z(:,jzf) = Z(:,jzf)*z; end
   end
elseif job == 2   
   % include all zeros 
   nus = nf+ni; ns = 0;
elseif job == 3   
   % include all infinite zeros 
   nus = ni; ns = nf;
elseif job == 4   
   % include all finite zeros 
   nus = nf; ns = ni;
elseif job == 5 
   % include all stable zeros 
   if nf == 0
      nus = 0; ns = ni;
      if discr
         nmszer = 0; nizer = 0;
      else
         nmszer = 0; nizer = ni;
      end
   else    
      % separate stable/unstable zeros
      % the finite zeros are in the trailing block
      izf = nc+ni+1:nc+ni+nf; jzf = mc+ni+1:mc+ni+nf;
      [AB(izf,jzf),EF(izf,jzf),q,z] = qz(AB(izf,jzf),EF(izf,jzf),'real');
      ev = ordeig(AB(izf,jzf),EF(izf,jzf)); 
      if discr
         %select = abs(ev) <= 1-offset;
         select = abs(ev) >= 1-offset;
         nizer = 0; 
         % nmszer = length(find(abs(ev(~select)) <= 1+offset));
         nmszer = length(find(abs(ev(select)) <= 1+offset));
      else
         %select = (real(ev) <= -offset);
         select = (real(ev) >= -offset);
         nizer = ni;
         % nmszer = length(find(real(ev(~select)) <= offset));
         nmszer = length(find(real(ev(select)) <= offset));
      end
      %nus = length(select(select==true)); 
      nus = length(select(select==false)); 
      [AB(izf,jzf),EF(izf,jzf),q,z] = ordqz(AB(izf,jzf),EF(izf,jzf),q,z,select); 
      ns = nf-nus+ni; 
      izf1 = 1:(nc+ni);
      j4 = mc+ni+nf+1:nm; 
      AB(izf,j4) = q*AB(izf,j4); EF(izf,j4) = q*EF(izf,j4);
      AB(izf1,jzf) = AB(izf1,jzf)*z; EF(izf1,jzf) = EF(izf1,jzf)*z;
      if withQ, Q(:,izf) = Q(:,izf)*q'; end
      if withZ, Z(:,jzf) = Z(:,jzf)*z; end
   end
elseif job == 6 
   % include all strictly unstable zeros and infinite zeros
   if nf == 0
      nus = ni; ns = 0;
      if discr
         nmszer = 0; nizer = 0;
      else
         nmszer = 0; nizer = ni;
      end
   else    
      % the finite zeros are in the leading block
      izf = nc+1:nc+nf; jzf = mc+1:mc+nf;
      [AB(izf,jzf),EF(izf,jzf),q,z] = qz(AB(izf,jzf),EF(izf,jzf),'real');
      ev = ordeig(AB(izf,jzf),EF(izf,jzf)); 
      if discr
         select = abs(ev) <= 1-offset;
         nizer = 0;
         nmszer = length(find(abs(ev(~select)) <= 1+offset));
      else
         select = (real(ev) <= -offset);
         nizer = ni;
         nmszer = length(find(real(ev(~select)) <= offset));
      end
      ns = length(select(select==true)); 
      [AB(izf,jzf),EF(izf,jzf),q,z] = ordqz(AB(izf,jzf),EF(izf,jzf),q,z,select); 
   
      % apply q to the trailing columns including infinite zeros
      j4 = mc+nf+1:nm;  
      izf1 = 1:nc;
      nus = nf-ns+ni;   % include infinite zeros among unstable zeros
      AB(izf,j4) = q*AB(izf,j4); EF(izf,j4) = q*EF(izf,j4);
      AB(izf1,jzf) = AB(izf1,jzf)*z; EF(izf1,jzf) = EF(izf1,jzf)*z;
      if withQ, Q(:,izf) = Q(:,izf)*q'; end
      if withZ, Z(:,jzf) = Z(:,jzf)*z; end
   end
end

% Step 3: Compress [E33 F3] to [E33 F3]*U2 = [E2 0] with E2(nl+nus,nl+nus) invertible.
%         Compute  [A33 B3]*U2 = [A2 B2].
%         Compute  [0 D1]*U2 = [C2 D2]  with D2 monic.

ir = nc+ns+1:n1; jr = mc+ns+1:n1m;
nbl = nl+nus;
if nbl
   ir1 = 1:(nc+ns); 
   [q,EF(ir,jr),z] = svd(EF(ir,jr));
   j4 = n1m+1:nm;  
   AB(ir,[jr j4]) = q'*AB(ir,[jr j4]); EF(ir,j4) = q'*EF(ir,j4);
   AB([ir1 ir],jr) = AB([ir1 ir],jr)*z; EF(ir1,jr) = EF(ir1,jr)*z;
   CD(:,jr) = CD(:,jr)*z; 
   if withQ, Q(:,ir) = Q(:,ir)*q; end
   if withZ, Z(:,jr) = Z(:,jr)*z; end
end 

% if withQ && withZ
%    norm(Q'*[a b]*Z-AB), norm(Q'*[e zeros(n,m)]*Z-EF), norm([c d]*Z-CD), 
% end

%
%  Build the matrices of the transformed system matrix pencil
%
At = [AB;CD]; Et = [EF; zeros(size(c,1),nm)];
mrg = mc+ns; mbl = n1m-mrg-nbl; % normal rank
dimsc = [mrg nbl mbl nsinf nmszer nizer];    % column dimensions of blocks

if nargout > 3 && ~withQ
    Q = [];
end

% if withQ && withZ
%    norm([Q'*[a b];c d]*Z-At), norm([Q'*[e zeros(n,m)]; zeros(size(CD))]*Z-Et),
% end

% end GSKLF
end

