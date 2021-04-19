function [syslnull,info] = glnull(sys,options)
%GLNULL  Left nullspace basis of a transfer function matrix.
%  [SYSLNULL,INFO] = GLNULL(SYS,OPTIONS) computes for a partitoned
%  system SYS = [ SYS1 SYS2 ], the left nullspace basis NL1 of SYS1.
%  In terms of the corresponding transfer-function matrices
%  we have NL1(lambda)*SYS1(lambda) = 0.
%  The returned SYSLNULL contains the compound system [ NL1  NL1*SYS2 ].
%  The system SYS = [ SYS1 SYS2 ] has the partitioned descriptor
%  realization (A-lambda E,[B1,B2],C,[D1,D2]) and the resulting
%  SYSLNULL = [ NL1  NL1*SYS2 ] has the partitoned observable descriptor
%  realization (Al-lambda*El,[Bl1,Bl2],Cl,[Dl1,Dl2]). For a minimal proper
%  left nullspace basis, the pencil [Al-lambda*El; Cl] is obtained
%  in an observability staircase form, with El upper triangular and
%  invertible. For a simple proper minimal basis, Al, El and Cl are block
%  diagonal and El is upper triangular. The i-th basis vector NL1(i,:) can
%  be explicitly constructed as (Ali-lambda*Eli,Bl1,Cli,Dl1(i,:)), where
%  Ali, Eli and Cli are the i-th diagonal blocks of Al, El, and Cl,
%  respectively, and Dl1(i,:) is the i-th row of Dl1. The corresponding
%  NL1(i,:)*SYS2 is constructed as (Ali-lambda*Eli,Bl2,Cli,Dl2(i,:)).
%  The OPTIONS structure allows to specify various user options,
%  as follows:
%  OPTIONS.tol     - tolerance for rank determinations
%                    (Default: internally computed)
%  OPTIONS.m2      - M2, the number of inputs of SYS2 (Default: M2 = 0)
%  OPTIONS.simple  - specifies the option to compute a simple proper basis
%                    true  - compute a simple basis; the orders of the
%                            basis vectors are provided in INFO.deg
%                    false - no simple basis computed (default)
%  OPTIONS.coinner - specifies the option to compute a coinner proper basis
%                    true  - compute an coinner basis; if the option for
%                            simple basis has been selected then each
%                            basis vector results coinner and the orders of
%                            the basis vectors are provided in INFO.deg
%                            (the resulting basis may not be coinner)
%                    false - no coinner basis computed (default)
%  OPTIONS.offset  - stability boundary offset OFFSET, to specify the 
%                    marginally stable finite zeros as follows: 
%                    in the continuous-time case are the finite zeros 
%                    having real parts in the interval [-OFFSET, OFFSET], 
%                    while in the discrete-time case having moduli in the 
%                    interval [1-OFFSET,1+OFFSET] 
%                    (Default: OFFSET = sqrt(eps)).
%  OPTIONS.tcond   - maximum allowed value for the condition numbers of the
%                    employed non-orthogonal transformation matrices
%                    (Default: 10^4) (only used if OPTIONS.simple = true)
%  OPTIONS.sdeg    - specifies a prescribed stability degree for the
%                    eigenvalues of SYSLNULL (Default: [])
%  OPTIONS.poles   - specifies a complex conjugated set of desired poles
%                    to be assigned as eigenvalues of SYSLNULL
%                    (Default: [])
%  OPTIONS.balance - specifies the balancing option for the Riccati 
%                    equation solvers: 
%                    true  - apply balancing (Default)
%                    false - disable balancing  
%
%  INFO is a structure containing additional information, as follows:
%  INFO.nrank - normal rank of SYS1;
%  INFO.stdim - dimensions of the diagonal blocks of Al-lambda*El:
%               if OPTIONS.simple = false, these are the column dimensions
%               of the full column rank subdiagonal blocks of the pencil
%               [Al-lambda*El;Cl] in observability staircase form;
%               if OPTIONS.simple = true, these are the orders of the
%               state-space realizations of the proper rational vectors of
%               the computed simple proper rational left nullspace
%               basis of SYS1;
%  INFO.degs  - increasingly ordered degrees of the vectors of a
%               polynomial left nullspace basis of SYS1 representing
%               the left Kronecker indices of SYS1;
%               also the orders of the realizations of the proper rational
%               vectors of a simple proper rational left nullspace
%               basis of SYS1. If OPTIONS.simple = true, INFO.deg(i) is
%               the dimension of the i-th diagonal blocks of Al and El.
%  INFO.tcond - maximum of the condition numbers of the employed
%               non-orthogonal transformation matrices; a warning is
%               issued if INFO.tcond >= OPTIONS.tcond.
%  INFO.fnorm - the norm of the employed output injection used for
%               stabilization (if any of OPTION.sdeg or OPTIONS.pole is
%               not empty) or if the coinner basis option has been selected. 
%  Note: Large values of INFO.tcond and/or INFO.fnorm indicate possible
%        loss of numerical stability. 
%
%  Note: The resulting realization of SYSLNULL is minimal provided the
%  realization of SYS is minimal. However, NL1 is a minimal proper basis
%  only if the realization (A-lambda E,B1,C,D1) of SYS1 is minimal. In
%  this case, INFO.degs are the degrees of the vectors of a minimal
%  polynomial basis or, if OPTIONS.simple = true, of the resulting
%  minimal simple proper basis.
%  Notes: 
%  1. The resulting realization of SYSLNULL is minimal provided the
%     realization of SYS is minimal. However, NL1 is a minimal proper 
%     basis only if the realization (A-lambda E,B1,C,D1) of SYS1 is 
%     minimal. In this case, INFO.degs are the degrees of the vectors of a 
%     minimal polynomial basis or, if OPTIONS.simple = true, of the  
%     resulting minimal simple proper basis. 
%  2. If SYS2 has poles on the boundary domain of the stability which are  
%     not poles of SYS1, then there exists no coinner NL1 which stabilizes 
%     NL1*SYS2. 

%
%  See also GRNULL.

%  Copyright 2015-2018 A. Varga
%  Author:      A. Varga, 21.11.2015.
%  Revision(s): A. Varga, 25.11.2016, 08.07.2017,16.08.2018, 18.09.2018,
%                         07.07.2019. 
%
%  Method:
%  (a) Computation of a minimal proper left nullspace basis
%  Let (A-lambda E,[B1,B2],C,[D1,D2]) be a descriptor realization of the   
%  partitioned system SYS = [SYS1,SYS2].  A rational basis NL1(lambda) of
%  the left nullspace of SYS1(lambda) is computed as 
%                NL1(lambda) = NL1e(lambda)*[0; I],  
%  where NL1e(lambda) is a rational basis of the left nullspace of the  
%  system matrix        
%           S(lambda) = [ A-lambda*E  B1 ] 
%                       [      C      D1 ]
%  The left nullspace NL1e(lambda) is computed using the Kronecker-like   
%  staircase form of S(lambda) computed as             
%                         [ Ar-lambda*Er     *        ]
%        Q'*S(lambda)*Z = [     0        Al-lambda*El ] ,
%                         [     0           Cl        ]
%  where Q and Z are orthogonal transformation matrices, 
%  the subpencil Ar-lambda*Er contains the right Kronecker structure
%  and the regular part of S(lambda), and the subpencil [ Al-lambda*El; Cl]
%  contains the left Kronecker structure of S(lambda). 
%  [Al-lambda El; Cl] is obtained in an observability staircase form, with
%  El upper triangular and nonsingular. 
% 
%  The rational basis NL1(lambda) is determined with the observable  
%  state space realization (Al+K*Cl-lambda El,Bl1+K*Dl1,Cl,Dl1), where    
%  [ *; Bl1; Dl1 ] = Q'*[0; I] and K is a stabilizing output injection 
%  matrix. The resulting basis is row proper, that is, Dl1 has full row rank.
%
%  NL1*SYS2 is computed as (Al+K*Cl-lambda El,Bl2+K*Dl2,Cl,Dl2), where 
%   [ *; Bl2; Dl2 ] = Q'*[ B2; D2].
%
%  (b) Computation of a minimal simple proper left nullspace basis
%      The method of [3] is employed to compute a simple basis from a
%      minimal proper basis as computed at (a). 
%
%  (c) Computation of an coinner proper left nullspace basis
%      The co-outer factor Nlo of the co-outer-coinner factorization  
%      NL1 = Nlo*Nli is determined as Nlo = (Al+K*Cl-lambda*El,-K*H,Cl,H), 
%      where where K is the stabilizing output injection computed by 
%      solving an appropriate filter Riccati equation and H is a suitable 
%      feedthrough matrix (see [4]). The updating of 
%      NL1 <- inv(Nlo)*NL1 and NL1*SYS2 <- inv(Nlo)*NL1*SYS2 is explicitly
%      obtained as NL1 = (Al+K*Cl-lambda*El,Bl1+K*Dl1,H\Cl,H\Dl1) and
%      SYS2*NR1 = (Al+K*Cl-lambda*El,Bl2+K*Dl2,H\Cl,H\Dl2).
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
   options = struct('tol',0,'m2',0,'simple',false);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

discr = (sys.Ts ~= 0);
if discr
    smax = 1-sqrt(eps); smin = 0;
else
    smax = -sqrt(eps);  smin = -inf;
end

% decode and check options

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

[p,m] = size(sys);
% number of inputs of SYS2
if isfield(options,'m2')
   m2 = options.m2;
   validateattributes(m2, {'double'},{'integer','scalar','>=',0,'<=',m},'','OPTIONS.m2') 
else
   m2 = 0;
end

% option for simple basis
if isfield(options,'simple')
   simple = options.simple; 
   validateattributes(simple, {'logical'},{'binary'},'','OPTIONS.simple') 
else
   simple = false;
end

% option for inner basis
if isfield(options,'coinner')
   coinner = options.coinner; 
   validateattributes(coinner, {'logical'},{'binary'},'','OPTIONS.coinner') 
else
   coinner = false;
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

if isempty(sdeg) && isempty(poles)
   stabilize = false;
else
   stabilize = true;
end

[a,b,c,d,e,Ts] = dssdata(sys);
n = size(a,1); 
m1 = m-m2; im1 = 1:m1; im2 = m1+1:m; 

% set absolute tolerance for Kronecker-like form computation
tolc = tol; 
if tol == 0
   tolc = n*n*eps*max([norm(a,1) norm(b,1),norm(c,inf) norm(d,1) norm(e,1)]); 
   if tolc == 0, tolc = eps; end  
   % nonzero value of tola enforced to fix an error in sl_klf
end

% separate the left Kronecker structure using the Kronecker-like form of 
% the system pencil S(z) = [ A-zE B1; C D1]:
%
%        Q'*S(z)*Z = [ Ar - z Er      *     ]
%                    [     0      Al - z El ] 
%                    [     0          Cl    ]
% 
A = [a b(:,im1); c d(:,im1)]; E = [e zeros(n,m1);zeros(p,n+m1)]; 
mode = 2; % perform separation of infinite part
qz = 2;   % accumulate only Q
[A,E,~,Q,nul,mul] = sl_klf(A',E',tolc,mode,qz); 
A = flip(flip(A,1),2)'; E = flip(flip(E,1),2)'; Q = flip(Q,2); 
mul = flip(mul); nul = flip(nul); 
     
% determine main dimensional parameters
nl = sum(mul);    % order of left structure (also order of Al)
pl = sum(nul)-nl; % number of basis elements
nrank = p - pl;   % normal rank of the transfer matrix SYS1
% 

% if nrank == 0 && nl == 0
%    % correct output from sl_klf by concatenating 1x0 blocks
%    mul = 0; nul = sum(nul);
% end

maxtcond = 1;     % set condition number of used transformations
maxfnorm = 0; 

if nrank == p 
   % full row rank: the null space is empty
   syslnull = ss(zeros(0,m)); degs = [];
else   
   % in the case nrank = 0, a null space basis is the identity matrix 
   % and therefore syslnull = [eye(p) sys(:,im2)] and degs = zeros(1,p);
   % however, to handle properly this case, the left Kronecker indices 
   % must be determined in degs; the basis is generally not minimal
   
   % compute the leftt Kronecker indices, which are also the degrees of 
   % the vectors of a polynomial left nullspace basis
   nvi = flip(nul-mul);  % there are nvi(i) vectors of degree i-1
   nb = length(nul);     % number of blocks in [ Al; Cl ]
   degs = zeros(1,pl);
   k = 1;
   for i = 1:nb
    for j=1:nvi(i)
        degs(k) = i-1;
        k = k+1;
    end
   end
   
   % form  BL1 = Q'*[0; I] = [ *; Bl; Dl ] and 
   %       BL2 = Q'*[ B2; D2] = [ *; Bl2; Dl2 ]    
   BL1 =  Q(n+1:end,:).'; 
   BL2 = Q.'*[b(:,im2); d(:,im2)]; 

   ial = n+p-nl-pl+1:n+p-pl; jal = n+m1-nl+1:n+m1; icl = n+p-pl+1:n+p; 

   al = A(ial,jal); bl = [ BL1(ial,:) BL2(ial,:) ];  el = E(ial,jal); 
   cl = A(icl,jal); dl = [ BL1(icl,:) BL2(icl,:) ]; 
   
   f = [];
   if coinner && m2
      % check stabilizability of SYS2*NR1
      [~,infoz] = gzero(dss(al,bl(:,1:p),cl,dl(:,1:p),el,Ts),tol,offset);
      if infoz.nfsbz 
         mes = sprintf(['No inner nullspace NL1 exists such that NL1*SYS2 is stable -',...
                        '\n         standard nullspace computed instead']);
         warning(mes)
         coinner = false;  % perform standard nullspace computation
      end
   end
       
   if coinner && ~simple
      % compute the co-outer factor Nlo of the co-outer-co-inner 
      % factorization NL1 = Nlo*Nli and form explicitly 
      % inv(Nro)*[NL1; NL1*SYS2]
      B = bl(:,1:p); D = dl(:,1:p); 
      if discr
         [X,~,f,ricrez] = dare(al.',cl.',B*B.',D*D.',B*D.',el.',balance); f = -f.';
         if ricrez < 0 || isnan(ricrez)
            if ricrez == -1 
               error('Solution of the DARE failed: Symplectic matrix has eigenvalues on the unit circle')
            else
               error('Solution of the DARE failed: no finite stabilizing solution exists')
            end
         end
         H = chol(D*D.'+cl*X*cl.','lower'); 
      else
         [~,~,f,ricrez] = care(al.',cl.',B*B.',D*D.',B*D.',el.',balance); f = -f.';
         if ricrez < 0 || isnan(ricrez)
            if ricrez == -1 
               error('Solution of the CARE failed: Hamiltonian has jw-axis eigenvalues')
            else
               error('Solution of the CARE failed: no finite stabilizing solution exists')
            end
         end                                                                                                 
         [~,H] = qr(D.'); H = H(1:pl,:).';                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
      end
      maxfnorm = norm(f,'fro');
      syslnull = dss(al+f*cl,bl+f*dl, H\cl, H\dl,el,Ts); 
  else
      if stabilize && ~coinner
         opt_poles = struct('tol',tol);
         if ~simple || (simple && pl == 1)  
            % make poles stable
            if nl
               f = gsfstab(al.',el.',cl.',poles,sdeg,opt_poles).';
               al = al+f*cl; bl = bl+f*dl;
            end
         end
      end
      maxfnorm = norm(f,'fro');
      syslnull = dss(al,bl,cl,dl,el,Ts);
   end
   
   if pl > 1 && simple && max(degs)
      % compute a simple basis  
      sysn = ss(zeros(size(dl)),syslnull); 
      maxtcond = 1;
      i = pl;
      for k = 1:pl
          ip = [i, 1:i-1, i+1:pl];  % row permutation indices
          [sysnk,info1] = glmcover1(syslnull(ip,:),1,tol);
          if order(sysnk) < degs(k)
%             warning('GLNULL: Resulting order from GLMCOVER1 less than expected')
             degs(k) = sum(info1.stdim);
          elseif order(sysnk) > degs(k)
             warning('GLNULL: Resulting order from GLMCOVER1 larger than expected')
          end
          maxtcond = max([maxtcond; info1.fnorm; info1.tcond]);
          [al,bl,cl,dl,el] = dssdata(sysnk);
          if coinner
             B = bl(:,1:p); D = dl(:,1:p); 
             if discr
                [X,~,f,ricrez] = dare(al.',cl.',B*B.',D*D.',B*D.',el.',balance); f = -f.';
                if ricrez < 0
                   if ricrez == -1 
                      error('Solution of the DARE failed: Symplectic matrix has eigenvalues on the unit circle')
                else
                   error('Solution of the DARE failed: no finite stabilizing solution exists')
                   end
                end
                H = chol(D*D.'+cl*X*cl.','lower'); 
             else
                [~,~,f,ricrez] = care(al.',cl.',B*B.',D*D.',B*D.',el.',balance); f = -f.';
                if ricrez < 0
                   if ricrez == -1 
                      error('Solution of the CARE failed: Hamiltonian has jw-axis eigenvalues')
                   else
                      error('Solution of the CARE failed: no finite stabilizing solution exists')
                   end
                end                                                                                                 
                [~,H] = qr(D.'); H = H(1,1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
             end
             maxfnorm = max([maxfnorm; norm(f,'fro')]);
             sysnk = dss(al+f*cl,bl+f*dl, H\cl, H\dl,el,Ts);  
          else
             if stabilize  
                f = gsfstab(al.',el.',cl.',poles,sdeg,opt_poles).';
                sysnk.a = al + f*cl; sysnk.b = sysnk.b + f*sysnk.d;
                maxfnorm = max([maxfnorm; norm(f,'fro')]);
             end
          end
          sysn(k,:) = sysnk;
          i = i-1;
      end
      if maxtcond > tcond
          disp(['Possible loss of numerical stability',...
               ' due to ill-conditioned transformations'])
      end
      syslnull = sysn;
   end
end


if nargout == 2
   if simple 
      info = struct('nrank',nrank,'degs',degs,'stdim',degs,...
                   'tcond',maxtcond,'fnorm',maxfnorm);
   else
      info = struct('nrank',nrank,'degs',degs,'stdim',mul(mul>0)',...
                   'tcond',maxtcond,'fnorm',maxfnorm);
   end
end

% end GLNULL
end

