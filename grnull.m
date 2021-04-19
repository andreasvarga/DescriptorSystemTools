function [sysrnull,info] = grnull(sys,options)
%GRNULL  Right nullspace basis of a transfer function matrix.
%   [SYSRNULL,INFO] = GRNULL(SYS,OPTIONS) computes for a partitoned 
%   system SYS = [ SYS1; SYS2 ], the right nullspace basis NR1 of SYS1. 
%   In terms of the corresponding transfer-function matrices 
%   we have SYS1(lambda)*NR1(lambda) = 0. 
%   The returned SYSRNULL contains the compound system [ NR1; SYS2*NR1 ]. 
%   The system SYS = [ SYS1; SYS2 ] has the partitioned descriptor  
%   realization (A-lambda*E,B,[C1;C2],[D1;D2]) and the resulting 
%   SYSRNULL = [ NR1; SYS2*NR1 ] has the partitoned controllable descriptor  
%   realization (Ar-lambda*Er,Br,[Cr1;Cr2],[Dr1;Dr2]). For a minimal proper 
%   right nullspace basis, the pencil [Br Ar-lambda*Er] is obtained 
%   in an controllability staircase form, with Er upper triangular and 
%   invertible. For a simple proper minimal basis, Ar, Er and Br are block 
%   diagonal and Er is upper triangular. The i-th basis vector NR1(:,i) can 
%   be explicitly constructed as (Ari-lambda*Eri,Bri,Cr1,Dr1(:,i)), where 
%   Ari, Eri and Bri are the i-th diagonal blocks of Ar, Er, and Br, 
%   respectively, and Dr1(:,i) is the i-th column of Dr1. The corresponding
%   SYS2*NR1(:,i) can be constructed as (Ari-lambda*Eri,Bri,Cr2,Dr2(:,i)).
%   The OPTIONS structure allows to specify various user options, 
%   as follows:
%   OPTIONS.tol     - tolerance for rank determinations 
%                     (Default: internally computed)
%   OPTIONS.p2      - P2, the number of outputs of SYS2 (Default: P2 = 0) 
%   OPTIONS.simple  - specifies the option to compute a simple proper basis
%                     true  - compute a simple basis; the orders of the  
%                             basis vectors are provided in INFO.deg
%                     false - no simple basis computed (default) 
%   OPTIONS.inner   - specifies the option to compute an inner proper basis
%                     true  - compute an inner basis; if the option for
%                             simple basis has been selected then each 
%                             basis vector results inner and the orders of   
%                             the basis vectors are provided in INFO.deg 
%                             (the resulting basis may not be inner) 
%                     false - no inner basis computed (default) 
%   OPTIONS.offset  - stability boundary offset OFFSET, to specify the 
%                     marginally stable finite zeros as follows: 
%                     in the continuous-time case are the finite zeros 
%                     having real parts in the interval [-OFFSET, OFFSET], 
%                     while in the discrete-time case having moduli in the 
%                     interval [1-OFFSET,1+OFFSET] 
%                     (Default: OFFSET = sqrt(eps)).
%   OPTIONS.tcond   - maximum allowed value for the condition numbers of  
%                     the employed non-orthogonal transformation matrices 
%                     (Default: 10^4) (only used if OPTIONS.simple = true)
%   OPTIONS.sdeg    - specifies a prescribed stability degree for the 
%                     eigenvalues of SYSRNULL (Default: [])
%   OPTIONS.poles   - specifies a complex conjugated set of desired poles
%                     to be assigned as eigenvalues of SYSRNULL 
%                     (Default: []) 
%   OPTIONS.balance - specifies the balancing option for the Riccati 
%                     equation solvers: 
%                     true  - apply balancing (Default)
%                     false - disable balancing  
%   
%   INFO is a structure containing additional information, as follows: 
%   INFO.nrank - normal rank of SYS1;
%   INFO.stdim - dimensions of the diagonal blocks of Ar-lambda*Er: 
%                if OPTIONS.simple = false, these are the column dimensions
%                of the full row rank diagonal blocks of the pencil 
%                [Br Ar-lambda*Er] in controllability staircase form; 
%                if OPTIONS.simple = true, these are the orders of the 
%                state-space realizations of the proper rational vectors of 
%                the computed simple proper rational right nullspace basis 
%                of SYS1;
%   INFO.degs  - increasingly ordered degrees of the vectors of a  
%                polynomial right nullspace basis of SYS1, representing  
%                the right Kronecker indices of SYS1;
%                also the orders of the realizations of the proper rational 
%                vectors of a simple proper rational right nullspace 
%                basis of SYS1. If OPTIONS.simple = true, INFO.deg(i) is
%                the dimension of the i-th diagonal blocks of Al and El.
%   INFO.tcond - maximum of the condition numbers of the employed 
%                non-orthogonal transformation matrices; a warning is 
%                issued if INFO.tcond >= OPTIONS.tcond.
%   INFO.fnorm - the norm of the employed state-feedback used for 
%                stabilization (if any of OPTION.sdeg or OPTIONS.pole is 
%                not empty) or if the inner basis option has been selected.  
%   Note: Large values of INFO.tcond and/or INFO.fnorm indicate possible
%         loss of numerical stability. 
%
%   Notes: 
%   1. The resulting realization of SYSRNULL is minimal provided the
%      realization of SYS is minimal. However, NR1 is a minimal proper 
%      basis only if the realization (A-lambda E,B,C1,D1) of SYS1 is 
%      minimal. In this case, INFO.degs are the degrees of the vectors of a 
%      minimal polynomial basis or, if OPTIONS.simple = true, of the  
%      resulting minimal simple proper basis. 
%   2. If SYS2 has poles on the boundary of the stability domain which are  
%      not poles of SYS1, then there exists no inner NR1 which stabilizes 
%      SYS2*NR1. 
%
%   See also GLNULL. 

%  Copyright 2016-2018 A. Varga
%  Author:      A. Varga, 24.10.2016.
%  Revision(s): A. Varga, 25.11.2016, 08.07.2017, 17.08.2018, 17.09.2018,
%                         07.07.2019. 
%
%  Method:
%  (a) Computation of a minimal proper right nullspace basis
%  Let (A-lambda E,B,[C1;C2],[D1;D2]) be a descriptor realization of the   
%  partitioned system SYS = [SYS1;SYS2].  A rational basis NR1(lambda) of
%  the left nullspace of SYS1(lambda) is computed as 
%                NLR(lambda) = [ 0 I] *NR1e(lambda),  
%  where NR1e(lambda) is a rational basis of the right nullspace of the  
%  system matrix        
%           S(lambda) = [ A-lambda*E  B1 ] 
%                       [      C      D1 ]
%  The left nullspace NL1e(lambda) is computed using the Kronecker-like   
%  staircase form of S(lambda) computed as             
%        Q'*S(lambda)*Z = [ Br  Ar-lambda*Er     *        ]
%                         [ 0      0         Al-lambda*El ] 
%  where Q and Z are orthogonal transformation matrices, 
%  the subpencil Al-lambda*El contains the left Kronecker structure
%  and regular part of S(lambda), and the subpencil [Br  Ar-lambda*Er]
%  contains the right Kronecker structure of S(lambda). 
%  [Br  Ar-lambda*Er] is obtained in an controllability staircase form, 
%  with Er upper triangular and nonsingular. 
% 
%  The rational basis NR1(lambda) is determined with the controllable  
%  state-space realization (Ar+Br*F-lambda*Er,Br,Cr1+Dr1*F,Dr1), where    
%  [ * Cr1 Dr1 ] = [0 I]*Z and F is a stabilizing state feedback. 
%  The resulting basis is column proper, that is, Dr1 has full column rank.
%
%  SYS2*NR1 is computed as (Ar+Br*F-lambda*Er,Br,Cr2+Dr2*F,Dr2), where 
%   [ * Cr2 Dr2 ] = [ C2 D2]*Z.
%
%  (b) Computation of a minimal simple proper right nullspace basis
%      The method of [3] is emloyed to compute a simple basis from a
%      minimal proper basis as computed at (a). 
%
%  (c) Computation of an inner proper right nullspace basis
%      The outer factor Nro of the inner-outer factorization NR1 = Nri*Nro 
%      is determined as Nro = (Ar+Br*F-lambda*Er,Br,-H*F,H), where F is 
%      the stabilizing state feedback computed by solving an appropriate 
%      control Riccati equation and H is a suitable feedthrough matrix 
%      (see [4]). The updating of NR1 <- NR1*inv(Nro) and 
%      SYS2*NR1 <- SYS2*NR1*inv(Nro) is explicitly obtained as 
%      NR1 = (Ar+Br*F-lambda*Er,Br/H,Cr1+Dr1*F,Dr1/H) and
%      SYS2*NR1 = (Ar+Br*F-lambda*Er,Br/H,Cr2+Dr2*F,Dr2/H).
%      For a simple basis, the above approach is applied separately to each
%      basis vector. 
%
%  References:
%  [1] Beelen, T.G.J.:
%      New algorithms for computing the Kronecker structure of a pencil 
%      with applications to systems and control theory. 
%      Ph. D. Thesis, Technical University Eindhoven, 1987.
%  [2] Varga, A.:
%      On computing least order fault detectors using rational nullspace bases. 
%      IFAC SAFEPROCESS'03 Symposium, Washington DC, USA, 2003.
%  [3] Varga, A.:
%      On computing nullspace bases – a fault detection perspective. 
%      Proc. IFAC 2008 World Congress, Seoul, Korea, pages 6295–6300, 2008.
%  [4] K. Zhou, J. C. Doyle, and K. Glover. 
%      Robust and Optimal Control. Prentice Hall, 1996.

narginchk(1,2)
nargoutchk(0,2)

% check for state-space system
if ~isa(sys,'ss')
   error('The input system SYS must be an SS object')
end

if nargin == 1
   options = struct('tol',0,'p2',0,'simple',false);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

discr = (sys.Ts ~= 0);
if discr
    smax = 1-sqrt(eps); smin = 0;
else
    smax = -sqrt(eps);  smin = -inf;
end

% decode options

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

% number of outputs of SYS2
[p,m] = size(sys);
if isfield(options,'p2')
   p2 = options.p2;
   validateattributes(p2, {'double'},{'integer','scalar','>=',0,'<=',p},'','OPTIONS.p2') 
else
   p2 = 0;
end

% option for simple basis
if isfield(options,'simple')
   simple = options.simple; 
   validateattributes(simple, {'logical'},{'binary'},'','OPTIONS.simple') 
else
   simple = false;
end

% option for inner basis
if isfield(options,'inner')
   inner = options.inner; 
   validateattributes(inner, {'logical'},{'binary'},'','OPTIONS.inner') 
else
   inner = false;
end

% offset for the stability region boundary
if isfield(options,'offset')
   offset = options.offset;
   validateattributes(offset, {'double'},{'real','scalar','>',0,'<',1},'','OPTIONS.offset') 
else
   offset = sqrt(eps);
end

% maximum allowed condition number 
if isfield(options,'tcond')
   tcond = options.tcond;
   validateattributes(tcond, {'double'},{'real','scalar','>=',1},'','OPTIONS.tcond') 
else
   tcond = 1.e4;
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

% desired poles
if isfield(options,'poles')
   poles = options.poles;
   if ~isempty(poles)
      validateattributes(poles, {'double'},{'vector'},'','OPTIONS.poles')
   end
else
   poles = [];
end
if ~isempty(poles)
    t = poles;
    if (min(size(t)) ~= 1 || ~isequal(sort(t(imag(t)>0)),sort(conj(t(imag(t)<0)))) )
        error('OPTIONS.poles must be a self-conjugated complex vector')
    end 
end  

% balancing option for Riccati equation solvers
if isfield(options,'balance')
   balopt = options.balance;
   validateattributes(balopt, {'logical'},{'binary'},'','OPTIONS.balance')
   if balopt
      balance = 'balance';
   else
      balance = 'nobalance';
   end    
else
   balance = 'balance';
end

% set stabilization option
if isempty(sdeg) && isempty(poles)
   stabilize = false;
else
   stabilize = true;
end

[a,b,c,d,e,Ts] = dssdata(sys);
n = size(a,1); 
p1 = p-p2; ip1 = 1:p1; ip2 = p1+1:p;

% set absolute tolerance for Kronecker-like form computation
tolc = tol; 
if tol == 0
   tolc = n*n*eps*max([norm(a,1) norm(b,1),norm(c,inf) norm(d,1) norm(e,1)]); 
   if tolc == 0, tolc = eps; end  
   % nonzero value of tola enforced to fix an error in sl_klf
end

% separate the left Kronecker structure using the Kronecker-like form of 
% the system pencil S(z) = [ A-zE B; C1 D1]:
%
%        Q'*S(z)*Z = [ Br Ar - z Er      *     ]
%                    [ 0     0       Al - z El ] 
% 
A = [a b ; c(ip1,:) d(ip1,:)]; E = [e zeros(n,m);zeros(p1,n+m)]; 
mode = 2; % perform separation of infinite part
qz = 2;   % accumulate only Z
[A,E,~,Z,nur,mur] = sl_klf(A,E,tolc,mode,qz); 
     
% determine main dimensional parameters
nr = sum(mur);    % order of left structure (also order of Ar)
mr = sum(nur)-nr; % number of basis elements
nrank = m - mr;   % normal rank of the transfer matrix SYS1

maxtcond = 1;     % set condition number of used transformations
maxfnorm = 0; 

if nrank == m 
   % full row rank: the null space is empty
   sysrnull = ss(zeros(m+p2,0)); degs = [];
else   
   % in the case nrank = 0, a null space basis is the identity matrix 
   % and therefore sysrnull = [eye(p); sys(ip2,:)] and degs = zeros(1,m);
   % however, to handle properly this case, the right Kronecker indices 
   % must be determined in degs; the basis is generally not minimal
   
   % compute the right Kronecker indices, which are also the degrees of 
   % the vectors of a polynomial right nullspace basis
   nvi = nur-mur;    % there are nvi(i) vectors of degree i-1
   nb = length(nur); % number of blocks in [ Br Ar ]
   degs = zeros(1,mr);
   k = 1;
   for i = 1:nb
    for j=1:nvi(i)
        degs(k) = i-1;
        k = k+1;
    end
   end
   
   % form  CR1 = [0 I]*Z = [ * Cr1 Dr1 ] and 
   %       CR2 = [ C2 D2]*Z = [ *  Cr2 Dr2 ]    
   CR1 =  Z(n+1:end,:); 
   CR2 = [ c(ip2,:) d(ip2,:)]*Z; 

   iar = 1:nr; jar = mr+1:mr+nr; jbr = 1:mr; 

   ar = A(iar,jar); cr = [ CR1(:,jar); CR2(:,jar) ];  er = E(iar,jar); 
   br = A(iar,jbr); dr = [ CR1(:,jbr); CR2(:,jbr) ]; 
   f = [];
   if inner && p2
      % check stabilizability of SYS2*NR1
      [~,infoz] = gzero(dss(ar,br,cr(1:m,:),dr(1:m,:),er,Ts),tol,offset);
      if infoz.nfsbz 
         mes = sprintf(['No inner nullspace NR1 exists such that SYS2*NR1 is stable -',...
                        '\n         standard nullspace computed instead']);
         warning(mes)
         inner = false;  % perform standard nullspace computation
      end
   end
       
   if inner && ~simple
      % compute the outer factor Nro of the inner-outer factorization 
      % NR1 = Nri*Nro and form explicitly [ NR1; SYS2*NR1]*inv(Nro)
      C = cr(1:m,:); D = dr(1:m,:); 
      if discr
         [X,~,f,ricrez] = dare(ar,br,C.'*C,D.'*D,C.'*D,er,balance); f = -f;
         if ricrez < 0 || isnan(ricrez)
            if ricrez == -1 
               error('Solution of the DARE failed: Symplectic matrix has eigenvalues on the unit circle')
            else
               error('Solution of the DARE failed: no finite stabilizing solution exists')
            end
         end
         H = chol(D.'*D+br.'*X*br); 
      else
         [~,~,f,ricrez] = care(ar,br,C.'*C,D.'*D,C.'*D,er,balance); f = -f;
         if ricrez < 0 || isnan(ricrez)
            if ricrez == -1 
               error('Solution of the CARE failed: Hamiltonian has jw-axis eigenvalues')
            else
               error('Solution of the CARE failed: no finite stabilizing solution exists')
            end
         end                                                                                                 
         [~,H] = qr(D); H = H(1:mr,:);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
      end
      maxfnorm = norm(f,'fro');
      sysrnull = dss(ar+br*f,br/H,cr+dr*f,dr/H,er,Ts); 
  else
      if stabilize && ~inner
         opt_poles = struct('tol',tol);
         if ~simple || (simple && mr == 1)                                            
            % make poles stable
            if nr
               f = gsfstab(ar,er,br,poles,sdeg,opt_poles);
               ar = ar+br*f; cr = cr+dr*f;
            end
         end
      end
      maxfnorm = norm(f,'fro');
      sysrnull = dss(ar,br,cr,dr,er,Ts);
   end
   
   if mr > 1 && simple && max(degs)
      % compute a simple basis  
      sysn = ss(zeros(size(dr)),sysrnull); 
      maxtcond = 1; 
      for i = 1:mr
          im = [i, 1:i-1, i+1:mr];  % row permutation indices
          [sysni,info1] = grmcover1(sysrnull(:,im),1,tol);
          if order(sysni) < degs(i)
%             warning('GRNULL: Resulting order from GRMCOVER1 less than expected')
             degs(i) = sum(info1.stdim);
          elseif order(sysni) > degs(i)
             warning('GRNULL: Resulting order from GRMCOVER1 larger than expected')
          end
          maxtcond = max([maxtcond; info1.fnorm; info1.tcond]);
          [ar,br,cr,dr,er] = dssdata(sysni);
          if inner
             C = cr(1:m,:); D = dr(1:m,:); 
             if discr
                [X,~,f,ricrez] = dare(ar,br,C.'*C,D.'*D,C.'*D,er,balance); f = -f;
                if ricrez < 0
                   if ricrez == -1 
                      error('Solution of the DARE failed: Symplectic matrix has eigenvalues on the unit circle')
                else
                   error('Solution of the DARE failed: no finite stabilizing solution exists')
                   end
                end
                H = chol(D'*D+br'*X*br); 
             else
                [~,~,f,ricrez] = care(ar,br,C.'*C,D.'*D,C.'*D,er,balance); f = -f;
                if ricrez < 0
                   if ricrez == -1 
                      error('Solution of the CARE failed: Hamiltonian has jw-axis eigenvalues')
                   else
                      error('Solution of the CARE failed: no finite stabilizing solution exists')
                   end
                end                                                                                                 
                [~,H] = qr(D); H = H(1,1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
             end
             maxfnorm = max([maxfnorm; norm(f,'fro')]);
             sysni = dss(ar+br*f,br/H,cr+dr*f,dr/H,er,Ts); 
          else
             if stabilize  
                f = gsfstab(ar,er,br,poles,sdeg,opt_poles);
                sysni.a = ar + br*f; sysni.c = sysni.c + sysni.d*f;
                maxfnorm = max([maxfnorm; norm(f,'fro')]);
             end
          end
          sysn(:,i) = sysni;
      end
      if maxtcond > tcond
          disp(['Possible loss of numerical stability',...
               ' due to ill-conditioned transformations'])
      end
      sysrnull = sysn;
   end
end


if nargout == 2
   if simple 
      info = struct('nrank',nrank,'degs',degs,'stdim',degs,...
                   'tcond',maxtcond,'fnorm',maxfnorm);
   else
      info = struct('nrank',nrank,'degs',degs,'stdim',mur(mur>0)',...
                   'tcond',maxtcond,'fnorm',maxfnorm);
   end
end

% end GRNULL
end

