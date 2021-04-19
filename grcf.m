function [sysn,sysm] = grcf(sys,options)
%GRCF   Right coprime factorization with proper and stable factors.
%       [SYSN,SYSM] = GRCF(SYS,OPTIONS) calculates for the LTI descriptor
%       system SYS = (A-lambda E,B,C,D) with the transfer function matrix 
%                                        -1
%               G(lambda) =  C(lambda*E-A) B + D
%
%       a right coprime factorization
%                                              -1
%               G(lambda) = N(lambda)*M(lambda)   ,
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
%                        based controllability tests 
%                        (Default: tolmin = prod(size(B)) * eps(norm(B,1)))
%       OPTIONS.smarg  - sets the stability margin
%                        (Default: -sqrt(eps) for a continuous-time system 
%                         and 1-sqrt(eps) for a discrete-time system) 
%       OPTIONS.sdeg   - specifies a prescribed stability degree for the 
%                        eigenvalues of the factors
%                        (Default:  -0.05 for a continuous-time system and
%                                    0.95  for a discrete-time system) 
%       OPTIONS.poles  - specifies a complex conjugated set of desired 
%                        poles to be assigned for the factors (Default: []) 
%       OPTIONS.mindeg - minimum degree option of denominator SYSM
%                        true  - determine minimum degree denominator
%                        false - both factors have the same order (default) 
%       OPTIONS.mininf - specifies option for removal of simple infinite
%                        eigenvalues of SYSN and SYSM
%                        true  - remove simple infinite eigenvalues
%                        false - keep simple infinite eigenvalues (default) 
%
%  See also GLCF.

%  Author:    A. Varga 22.11.2015.
%  Revisions: A. Varga 21.01.2016, 28.09.2016, 08.12.2016, 12.06.2017.
%
%  Method:  The Procedure GRCF from [1] is implemented, which represents
%  an extension of the recursive factorization approach of [2] to cope with  
%  infinite poles. All infinite poles are assigned to finite real values. 
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

 
[At,Bt,Ct,Dt,Et,Ts] = dssdata(sys);

[n,m] = size(Bt); 
discr = (Ts ~= 0);
standsys = isequal(Et,eye(n));  
if discr
    smax = 1-sqrt(eps); smin = 0;
else
    smax = -sqrt(eps);  smin = -inf;
end

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
   tola = n*n*eps(max(1,norm(At,1)));
   tolb = n*m*eps(max(1,norm(Bt,1)));
else
   tola = tolc;
   tolb = tolc;
end

% stability margin
if isfield(options,'smarg')
   smarg = options.smarg;
   if ~isempty(smarg)
      validateattributes(smarg, {'double'},{'real','scalar','<=',smax,'>=',smin},'','OPTIONS.smarg') 
   end
else
   smarg = [];
end
if isempty(smarg)
   % set stability margin
   if discr
      smarg = 1-sqrt(eps);
   else
      smarg = -sqrt(eps);
   end
end    

% desired stability degree
if isfield(options,'sdeg')
   sdeg = options.sdeg;
   if ~isempty(sdeg)
      validateattributes(sdeg, {'double'},{'real','scalar','<=',smax,'>=',smin},'','OPTIONS.sdeg') 
   end
else
   sdeg = [];
end
if isempty(sdeg)
   % set desired stability degree
   if discr
      sdeg = 0.95;
   else
      sdeg = -0.05;
   end
end    

% desired poles
if isfield(options,'poles')
   poles = options.poles;
   if ~isempty(poles)
      validateattributes(poles, {'double'},{'vector'},'','OPTIONS.poles')
   end
else
   poles = [];
end

% standardize desired poles
if ~isempty(poles)
    tempr = poles(imag(poles)==0);
    tempi = sort(poles(imag(poles)>0));
    if ~isequal(tempi,sort(conj(poles(imag(poles)<0)))) 
        error('OPTIONS.poles must be a self-conjugated complex vector')
    end 
    ti = [tempi(:).'; conj(tempi(:).')];
    poles = [tempr(:); ti(:)];
    if (discr && ~isempty(find(abs(poles) > 1-sqrt(eps), 1))) || ...
       (~discr && ~isempty(find(real(poles) > -sqrt(eps), 1))) 
          error('The elements of OPTIONS.pole must lie in the stability region')
    end
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

% initialization of the recursive factorization 
CN = Ct; CM = zeros(m,n); DN =Dt; DM = eye(m);
% 
nrmBt = max(tolb,norm(Bt,1));
fnrmtol = 10000*max(1,norm(At,1))/nrmBt;  % set threshold for high feedback
fwarn = 0;

while nb
   select = false(nb,1); 
   if nb == 1 || abs(At(n,n-1)) == 0
      k = 1; 
   else
      k = 2;
   end
   kk = n-k+1:n;
   if standsys
      evb = ordeig(At(kk,kk));
   else
      evb = ordeig(At(kk,kk),Et(kk,kk));
   end
   if norm(Bt(kk,:)) <= tolb
      nb = nb-k; n = n-k; 
      if nbi, nbi = nbi-1; end
      i1 = 1:n; At = At(i1,i1); Bt = Bt(i1,:); 
      if ~standsys , Et = Et(i1,i1); end
      CN = CN(:,i1); CM = CM(:,i1); 
   else
      a2 = At(kk,kk); b2 = Bt(kk,:); 
      if standsys, e2 = eye(k); else e2 = Et(kk,kk);  end
      if ~standsys && k == 1 && nbi
         % move infinite eigenvalue to a stable (finite) position
         [gamma,poles] = eigselect1(poles,sdeg,evb(end));
         if isempty(gamma), gamma = sdeg; end
         [q2,r2]=qr(b2.',0); 
         f2 = -(q2/r2)*a2; W =  eye(m)-q2*q2.';
         % update matrices
         i1 = 1:n-1;
         At(i1,kk) = At(i1,kk)+Bt(i1,:)*f2; At(n,n) = gamma; Et(n,n)=1;
         CN(:,kk) = CN(:,kk)+DN*f2; CM(:,kk) = CM(:,kk)+DM*f2; 
         Bt(i1,:) = Bt(i1,:)*W; DN = DN*W; DM = DM*W;
         select(end) = true; i2 = n-nb+1:n; i1 = 1:n-nb;
         [At(i2,i2),Et(i2,i2),Q3,Z3] = ...
                    ordqz(At(i2,i2),Et(i2,i2),eye(nb),eye(nb),select);
         At(i1,i2) = At(i1,i2)*Z3; Et(i1,i2) = Et(i1,i2)*Z3;
         Bt(i2,:) = Q3*Bt(i2,:); 
         CN(:,i2) = CN(:,i2)*Z3; CM(:,i2) = CM(:,i2)*Z3;
         nb = nb-1; nbi = nbi-1;
      else
         if k == 1 
            % assign a single eigenvalue 
            [gamma,poles] = eigselect1(poles,sdeg,evb(end));
            if isempty(gamma)
               % no real pole available, adjoin a new 1x1 block if possible
               if nb == 1
                 % incompatible poles with the eigenvalue structure
                 % assign the last real pole to SDEG (if possible)
                 warning(['No real eigenvalue available for assignment: assigning instead SDEG = ', num2str(sdeg)])
                 f2 = -b2\(a2-e2*sdeg);
               else
                 % adjoin a real block or interchange the last two blocks
                 k = 2; kk = n-k+1:n; 
                 if n == 2 || abs(At(n-1,n-2)) == 0
                    % adjoin blocks and update evb 
                    if standsys
                       evb = ordeig(At(kk,kk));
                    else
                       evb = ordeig(At(kk,kk),Et(kk,kk));
                    end
                 else
                    % interchange last two blocks
                    i1 = 1:n-3; i2 = n-2:n; 
                    if standsys
                      [Q3,At(i2,i2)] = ordschur(eye(3),At(i2,i2),[false,false,true]);
                      At(i1,i2) = At(i1,i2)*Q3; Bt(i2,:) = Q3.'*Bt(i2,:); 
                      Z(:,i2) = Z(:,i2)*Q3; 
                      CN(:,i2) = CN(:,i2)*Q3; CM(:,i2) = CM(:,i2)*Q3;
                    else
                      [At(i2,i2),Et(i2,i2),Q3,Z3] = ordqz(At(i2,i2),Et(i2,i2),eye(3),eye(3),[false,false,true]);
                      At(i1,i2) = At(i1,i2)*Z3; Et(i1,i2) = Et(i1,i2)*Z3;
                      Bt(i2,:) = Q3*Bt(i2,:); 
                      Q(i2,:) = Q3*Q(i2,:); Z(:,i2) = Z(:,i2)*Z3; 
                      CN(:,i2) = CN(:,i2)*Z3; CM(:,i2) = CM(:,i2)*Z3;
                    end
                    if standsys
                       evb = ordeig(At(kk,kk));
                    else
                       evb = ordeig(At(kk,kk),Et(kk,kk));
                    end
                 end
                 a2 = At(kk,kk); b2 = Bt(kk,:); 
                 if standsys, e2 = eye(k); else e2 = Et(kk,kk);  end
               end
            else
               f2 = -b2\(a2-e2*gamma);
            end
         end
         if k == 2
            % assign a pair of eigenvalues 
            [px,poles] = eigselect2(poles,sdeg,evb(end),discr);
            [f2,u,v] = galoc2(a2,e2,b2,px,tola,tolb);
            if isempty(f2)
               At(kk,kk) = u'*At(kk,kk); At(:,kk) = At(:,kk)*v;
               if ~standsys, 
                  Et(kk,kk) = u'*Et(kk,kk); Et(:,kk) = Et(:,kk)*v; 
               end
               Bt(kk,:) = u'*Bt(kk,:);
               Q(kk,:) = u'*Q(kk,:);
               CN(:,kk) = CN(:,kk)*v; CM(:,kk) = CM(:,kk)*v;
            end
         end
         if isempty(f2)
             % delete last uncontrollable eigenvalue
             i1 = 1:n-1; At = At(i1,i1); Bt = Bt(i1,:); 
             if ~standsys, Et = Et(i1,i1); end
             CN = CN(:,i1); CM = CM(:,i1); 
             nb = nb-1; n = n-1;
         else
             if norm(f2,inf) > fnrmtol
                fwarn = fwarn+1;
             end
             % update matrices
             At(:,kk) = At(:,kk)+Bt*f2; 
             CN(:,kk) = CN(:,kk)+DN*f2; 
             CM(:,kk) = CM(:,kk)+DM*f2; 
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
         end
      end
   end
end

if mininf && nsinf > 0
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
   warning('ACC:TEST','Possible loss of numerical accuracy due to high feedback gain')
end

% end GRCF
end