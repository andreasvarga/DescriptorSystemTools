function [sysn,sysm] = grcfid(sys,options)
%GRCFID Right coprime factorization with inner denominator.
%       [SYSN,SYSM] = GRCFID(SYS,OPTIONS) calculates for the LTI descriptor
%       system SYS = (A-lambda E,B,C,D) with the transfer function matrix 
%                                        -1
%               G(lambda) =  C(lambda*E-A) B + D
%
%       a right coprime factorization
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
%                        based controllability tests 
%                        (Default: tolmin = prod(size(B)) * eps(norm(B,1)))
%       OPTIONS.mindeg - minimum degree option of denominator SYSM
%                        true  - determine minimum degree denominator
%                        false - both factors have the same order (default) 
%       OPTIONS.mininf - specifies option for removal of simple infinite
%                        eigenvalues of SYSN and SYSM
%                        true  - remove simple infinite eigenvalues
%                        false - keep simple infinite eigenvalues (default) 
%
%       See also GLCFID.

%  Author:    A. Varga, 20.01.2016.
%  Revisions: A. Varga, 04.12.2016, 12.06.2017.
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

[At,Bt,Ct,Dt,Et,Ts] = dssdata(sys);

[n,m] = size(Bt); 
discr = Ts ~= 0;
standsys = isequal(Et,eye(n));  

% decode options and set default values

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end
if tol == 0 && ~standsys
   % set absolute tolerance for determination of rank of E 
   tol = n*n*eps(norm(Et,1));
end

% tolerance for controllability tests
if isfield(options,'tolmin')
   tolc = options.tolmin;
   validateattributes(tolc, {'double'},{'real','scalar','>=',0},'','OPTIONS.tolmin') 
else
   tolc = 0;
end
if tolc == 0 
   % set absolute tolerance for negligible elements in A and B 
   tolb = n*m*eps(max(1,norm(Bt,1)));
else
   tolb = tolc;
end

% set stability margin and tolerance for eigenvalues
toleig = sqrt(eps); 
if discr
   smarg = 1-toleig;
else
   smarg = -toleig;
end

% minimum degree denominator option
if isfield(options,'mindeg')
   mindeg = options.mindeg;
   validateattributes(mindeg, {'logical'},{'binary'},'','OPTIONS.mindeg') 
else
   mindeg = false;
end

% simple infinite eigenvalues removal option
if isfield(options,'mininf')
   mininf = options.mininf;
   validateattributes(mininf, {'logical'},{'binary'},'','OPTIONS.mininf') 
else
   mininf = false;
end

% reduce the pair (A,E) to the specially ordered generalized real Schur
% form (At,Et) = (Q'*A*Z,Q'*E*Z), with
%
%             [ Ai  *   *   *  ]          [ 0  *    *   *  ]
%        At = [  0  Ag  *   *  ]  ,  Et = [ 0  Eg   *   *  ] ,
%             [  0  0  Abf  *  ]          [ 0  0   Ebf  *  ]
%             [  0  0   0  Abi ]          [ 0  0    0  Ebi ]
% 
% where
%    (Ai,0)    contains the firt order inifinite eigenvalues
%    (Ag,Eg)   contains the "good" finite eigenvalues
%    (Abf,Ebf) contains the "bad" finite eigenvalues
%    (Abi,Ebi) contains the "bad" (higher order) infinite eigenvalues

options_gsorsf = struct('tol',tol,'disc',discr,'smarg',smarg,...
    'sepinf',true,'fast',true);
[At,Et,Q,Z,dims] = gsorsf(At,Et,options_gsorsf);
Bt = Q.'*Bt; Ct = Ct*Z; 
nsinf = dims(1); ng = dims(2); nb = dims(3); nbi = dims(4);

% finish if SYS is stable and Et is invertible or
% SYS is stable and simple infinite eigenvalues have not to be removed 
if nb+nsinf == 0 || (nb == 0 && ~mininf)
   sysn = sys;
   if nargout == 2
      sysm = ss(eye(size(Dt,2))); sysm.Ts = Ts;
   end
   return
end

i3 = n-nb+1:n;
evb = ordeig(At(i3,i3),Et(i3,i3));

% initialization of the recursive factorization 
CN = Ct; CM = zeros(m,n); DN =Dt; DM = eye(m);
% 
nrmBt = max(1,norm(Bt,1));
fnrmtol = 100*max(1,norm(At,1))/nrmBt;
fwarn = 0;

while nb
   select = false(nb,1); 
   if imag(evb(end)) ~= 0
      k = 2; 
   else
      k = 1;
   end
   kk = n-k+1:n;
   if norm(Bt(kk,:)) < tolb
      nb = nb-k; n = n-k; 
      if nbi, nbi = nbi-1; end
      i1 = 1:n; At = At(i1,i1); Bt = Bt(i1,:); 
      if ~standsys , Et = Et(i1,i1); end
      CN = CN(:,i1); CM = CM(:,i1); 
      evb = evb(1:nb);
   else
      a2 = At(kk,kk); b2 = Bt(kk,:); 
      if standsys, e2 = eye(k); else e2 = Et(kk,kk);  end
      if ~standsys && k == 1 && nbi
         if ~discr
             error('Improper continuous-time system')
         end
         % move infinite eigenvalue to zero
         [q2,r2]=qr(b2.',0); 
         f2 = -(q2/r2)*a2; W =  eye(m)-q2*q2.';
         % update matrices
         i1 = 1:n-1;
         At(i1,kk) = At(i1,kk)+Bt(i1,:)*f2; At(n,n) = 0; Et(n,n)=-a2;
         Bt(i1,:) = Bt(i1,:)*W; 
         % enforce positive Et(n,n) (needed for MATLAB 2020b)
         if a2 > 0, Et(n,n) = -Et(n,n); Bt(n,:) = - Bt(n,:); end
         CN(:,kk) = CN(:,kk)+DN*f2; CM(:,kk) = CM(:,kk)+DM*f2; 
         DN = DN*W; DM = DM*W;
         select(end) = true; i2 = n-nb+1:n; i1 = 1:n-nb;
         [At(i2,i2),Et(i2,i2),Q3,Z3] = ...
                    ordqz(At(i2,i2),Et(i2,i2),eye(nb),eye(nb),select);
         At(i1,i2) = At(i1,i2)*Z3; Et(i1,i2) = Et(i1,i2)*Z3;
         Bt(i2,:) = Q3*Bt(i2,:); 
         CN(:,i2) = CN(:,i2)*Z3; CM(:,i2) = CM(:,i2)*Z3;
         nb = nb-1; nbi = nbi-1;
         i3 = n-nb+1:n;
         evb = ordeig(At(i3,i3),Et(i3,i3));
      else
         % check existence conditions
         if discr
            if abs(abs(evb(end))-1) < toleig
               error('Eigenvalue(s) on the unit circle present')
            end
         else
            if abs(real(evb(end))) < toleig
               error('Eigenvalue(s) on the imaginary axis present')
            end  
         end
         if k == 1 
            % reflect single eigenvalue 
            if discr
               % solve e2*y*e2'-a2*y*a2'+b2*b2' for y=S*S'
               S = sl_glme(4,e2,a2,b2,[1 1],1);                
               f2 = -(S'\(S\(a2\b2))).'; % f2 = -b2.'/(y*a2.');
               At(kk,kk) = e2*(e2/a2); 
               x = (e2*S)\b2; w = inv(cholupdate(eye(m),x(:)));
            else
               % solve -a2*y*e2'-e2*y*a2'+b2*b2' for y=S*S'
               S = sl_glme(4,-a2,e2,b2,[0 1],1); 
               f2 = -(S'\(S\(e2\b2))).'; % f2 = -b2.'/(y*e2.');
               At(kk,kk) = -At(kk,kk);
               w = eye(m);
            end
            iupd = 1:n-1;
         else
            % reflect a pair of eigenvalues 
            if discr
                % solve e2*y*e2'-a2*y*a2'+b2*b2' for y=S*S'
                S = sl_glme(4,e2,a2,b2,[1 0],1); 
                f2 = -(S'\(S\(a2\b2))).'; % f2 = -b2.'/(y*a2.');
                x = (e2*S)\b2; 
                w = inv(cholupdate(cholupdate(eye(m),x(1,:)'),x(2,:)'));
             else
               % solve -a2*y*e2'-e2*y*a2'+b2*b2' for y=S*S'
                S = sl_glme(4,-a2,e2,b2,[0 0],1); 
                f2 = -(S'\(S\(e2\b2))).'; % f2 = -b2.'/(y*e2.');
                w = eye(m);
            end
            iupd = 1:n;
         end
         if norm(f2,inf) > fnrmtol
            fwarn = fwarn+1;
         end
         % update matrices
         At(iupd,kk) = At(iupd,kk)+Bt(iupd,:)*f2; 
         Bt = Bt*w; 
         CN(:,kk) = CN(:,kk)+DN*f2; 
         CM(:,kk) = CM(:,kk)+DM*f2; 
         DN = DN*w; DM = DM*w; 
         select(end-k+1:end) = true; i2 = n-nb+1:n; i1 = 1:n-nb;
         if standsys
            if k == 2
               % standardize Schur form (needed for MATLAB 2020b)
               [Q2,At(kk,kk)] = schur(At(kk,kk));
               At(1:n-2,kk) = At(1:n-2,kk)*Q2; Bt(kk,:) = Q2.'*Bt(kk,:); 
               CN(:,kk) = CN(:,kk)*Q2; CM(:,kk) = CM(:,kk)*Q2;
            end
            [Q3,At(i2,i2)] = ordschur(eye(nb),At(i2,i2),select);
            At(i1,i2) = At(i1,i2)*Q3; Bt(i2,:) = Q3.'*Bt(i2,:); 
            CN(:,i2) = CN(:,i2)*Q3; CM(:,i2) = CM(:,i2)*Q3;
         else 
            [At(i2,i2),Et(i2,i2),Q3,Z3] = ...
                    ordqz(At(i2,i2),Et(i2,i2),eye(nb),eye(nb),select);
            At(i1,i2) = At(i1,i2)*Z3; Et(i1,i2) = Et(i1,i2)*Z3;
            Bt(i2,:) = Q3*Bt(i2,:); 
            CN(:,i2) = CN(:,i2)*Z3; CM(:,i2) = CM(:,i2)*Z3;
         end
         nb = nb-k;
         i3 = n-nb+1:n;
         if standsys
            evb = ordeig(At(i3,i3));
         else
            evb = ordeig(At(i3,i3),Et(i3,i3));
         end
      end
   end
end

if mininf && nsinf 
    % remove the simple infinite eigenvalues
    i1 = 1:nsinf; i2 = nsinf+1:n;
    ca = CN(:,i1)/At(i1,i1); ee = Et(i1,i2)/Et(i2,i2);
    DN = DN-ca*Bt(i1,:)+ca*ee*Bt(i2,:); 
    CN = CN(:,i2)-ca*At(i1,i2)+ca*ee*At(i2,i2); CM = CM(:,i2);
    At = At(i2,i2); Et = Et(i2,i2); Bt = Bt(i2,:);
    n = n-nsinf; 
end

if standsys
  sysn = ss(At,Bt,CN,DN,Ts);
  if nargout == 2
    if mindeg
       i2 = ng+1:n;
       sysm = ss(At(i2,i2),Bt(i2,:),CM(:,i2),DM,Ts);
    else
       sysm = ss(At,Bt,CM,DM,Ts);
    end
  end
else
  sysn = dss(At,Bt,CN,DN,Et,Ts);
  if nargout == 2
    if mindeg
       i2 = ng+1:n;
       sysm = dss(At(i2,i2),Bt(i2,:),CM(:,i2),DM,Et(i2,i2),Ts);
    else
       sysm = dss(At,Bt,CM,DM,Et,Ts);
    end
  end
end

if fwarn 
   warning('Possible loss of numerical accuracy due to high feedback gain')
end

% end GRCFID
end



