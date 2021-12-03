function [f,info] = gsfstab(a,e,b,poles,sdeg,options)
%GSFSTAB  Generalized state-feedback stabilization.
%   F = GSFSTAB(A,E,B,POLES,SDEG,OPTIONS) computes a state-feedback 
%   matrix F such that the controllable finite generalized eigenvalues 
%   of the closed-loop pair (A+B*F,E) lie in the stability domain Cs 
%   specified by SDEG. For a continuous-time setting, with SDEG < 0, 
%   Cs is the set of complex numbers with real parts at most SDEG,
%   while for a discrete-time setting, with 0 <= SDEG < 1, Cs is the set 
%   of complex numbers with moduli at most SDEG. 
%   POLES is a self-conjugated complex vector, which contains a set of 
%   desired closed-loop generalized eigenvalues for the pair (A+B*F,E), 
%   to replace the finite generalized eigenvalues of (A,E) lying outside
%   of the specified stability domain Cs. If the number of specified  
%   eigenvalues in POLES is less than the number of controllable 
%   generalized eigenvalues of (A,E) outside of Cs, then the rest of 
%   generalized eigenvalues of (A+BF,E) are assigned to the nearest values 
%   on the boundary of Cs (defined by  SDEG) or are kept unmodified if 
%   SDEG = []. (Default values: POLES = [])
%   The default value for SDEG, used if SDEG = [] and POLES = [], is -0.2. 
%
%   The OPTIONS structure allows to specify user options, as follows:
%   OPTIONS.tol    - tolerance for rank determinations  
%                    (Default: internally computed) 
%   OPTIONS.sepinf - specifies the option for a preliminary separation 
%                    of the infinite eigenvalues:
%                    false - no separation of infinite generalized 
%                            eigenvalues (default)
%                    true  - perform preliminary separation of the infinite
%                            generalized eigenvalues from the finite ones.
%
%   [F,INFO] = GSFSTAB(A,E,B,POLES,SDEG,OPTIONS) provides additionally
%   in the INFO structure information on the closed-loop matrices 
%   Acl := Q*(A+B*F)*Z, Ecl := Q*E*Z, and Bcl := Q*B, where Q and Z are the  
%   orthogonal matrices used to obtain the pair (Acl,Ecl) in a generalized 
%   real Schur form (GRSF). The fields of the INFO structure are:
%   INFO.Acl      - contains Acl in a quasi-upper triangular form;
%   INFO.Ecl      - contains Ecl in an upper triangular form;
%   INFO.Bcl      - contains Bcl;
%   INFO.Q        - contains the orthogonal transformation matrix Q;
%   INFO.Z        - contains the orthogonal transformation matrix Z;
%   INFO.ninf     - number of infinite generalized eigenvalues of (A,E)
%   INFO.nfg      - number of finite generalized eigenvalues of (A,E) 
%                   lying in the stability domain Cs 
%   INFO.naf      - number of assigned finite generalized eigenvalues in Cs
%   INFO.nuf      - number of uncontrollable finite generalized eigenvalues 
%                   lying outside Cs 
%
%   The resulting matrices Acl, Ecl and Bcl have the forms
%
%          ( Ai   *   *   *  )         ( Ei   *   *   *  )          ( * )
%    Acl = ( 0   Afg  *   *  ) , Ecl = ( 0   Efg  *   *  ) ,  Bcl = ( * ) , 
%          ( 0    0  Afa  *  )         ( 0    0  Efa  *  )          ( * )
%          ( 0    0   0  Afu )         ( 0    0   0  Efu )          ( 0 )
%
%   where: (Ai,Ei) with Ai upper triangular and invertible and Ei upper 
%   triangular and nilpotent, contains the INFO.ninf infinite generalized 
%   eigenvalues of (A,E); (Afg,Efg), in GRSF, contains the INFO.nfg finite 
%   generalized eigenvalues of (A,E) in Cs; (Afa,Efa), in GRSF, contains 
%   the INFO.naf assigned finite generalized eigenvalues in Cs; and, 
%   (Afu,Efu), in GRSF, contains the uncontrollable finite generalized 
%   eigenvalues of (A,E) lying outside Cs. 

%  Author:      A. Varga, 01.02.2016
%  Revision(s): A. Varga, 04.11.2016, 13.06.2017, 27.11.2021. 
%
%  Method:  For a standard system pair (A,I) (for E = []), the Schur method
%  of [1] is used, while for a generalized system pair (A,E) the 
%  generalized Schur method of [2] is used.
%
%  References:
%  [1] A. Varga. 
%      A Schur method for pole assignment.
%      IEEE Trans. on Automatic Control, vol. 26, pp. 517-519, 1981.
%  [2] A. Varga. 
%      On stabilization methods of descriptor systems.
%      Systems & Control Letters, vol. 24, pp.133-138, 1995.
%

narginchk(3,6)
nargoutchk(0,2)

validateattributes(a, {'double'}, {'real','finite', 'square'});

n = size(a,1); 

if ~isempty(e)
    validateattributes(e, {'double'}, {'real','finite', 'square','size',[n,n]});
end

validateattributes(b, {'double'}, {'real','finite'});
m = size(b,2); 
validateattributes(b, {'double'}, {'size',[n,m]});
% quick exit if possible 
if n == 0
   f = zeros(m,n); 
   if nargout > 1
      info = struct('Acl',[],'Ecl',[],'Bcl',[],'Q',[],'Z',[],...
                 'ninf',0,'nfg',0,'naf',0,'nuf',0);
   end
   return
end

if m == 0
   f = zeros(m,n); 
   if nargout > 1
      info = struct('Acl',a,'Ecl',a,'Bcl',b,'Q',eye(n),'Z',eye(n),...
                 'ninf',0,'nfg',0,'naf',0,'nuf',n);
   end
   return
end

% desired poles
if nargin < 4
   poles = [];
else
   if ~isempty(poles)
      validateattributes(poles, {'double'},{'vector'},'','POLES')
   end
end

if nargin < 5 
    sdeg = []; 
end
if ~isempty(sdeg)
   validateattributes(sdeg, {'double'},{'real','scalar'},'','SDEG') 
end

if ~isempty(sdeg)
   discr = (sdeg >= 0);
   if discr && sdeg >= 1
      error('Invalid value for SDEG')
   end
else
   discr = false;
end

% standardize desired poles
if ~isempty(poles)
    tempr = poles(imag(poles)==0);
    tempi = sort(poles(imag(poles)>0));
    if ~isequal(tempi,sort(conj(poles(imag(poles)<0)))) 
        error('POLES must be a self-conjugated complex vector')
    end 
    ti = [tempi(:).'; conj(tempi(:).')];
    poles = [tempr(:); ti(:)];
    if ~isempty(sdeg)
       if (discr && ~isempty(find(abs(poles) > sdeg, 1))) || (~discr && ~isempty(find(real(poles) > sdeg, 1))) 
          error('The elements of POLES must lie in the stability region')
       end
    end
end    


if nargin < 6
    options = struct('tol',0);
else
    validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

% set default values of SDEG if POLES = []
if isempty(sdeg)
    if isempty(poles) 
       sdeg = -0.2; 
       smarg = sdeg;         
    else
       smarg = -inf;
    end
else
    smarg = sdeg;  
end

% decode options and set default values

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end
if tol == 0 
   % set absolute tolerances for negligible elements in A and B 
   tola = n*n*eps(max(1,norm(a,1)));
   tolb = n*m*eps(max(1,norm(b,1)));
else
   tola = tol;
   tolb = tol;
end

% infinite eigenvalue separation
if isfield(options,'sepinf')
   sepinf = options.sepinf;
   validateattributes(sepinf, {'logical'},{'binary'},'','OPTIONS.sepinf') 
else
   sepinf = false;
end


standsys = isempty(e) || isequal(e,eye(n));  

% separate stable and unstable parts with respect to sdeg
if  standsys
    % compute orthogonal Q  such that
    %
    %      Z^T*A*Z = [ Ag   * ]
    %                [  0  Ab ]
    %
    % where Ag has eigenvalues within the stability degree region
    % and Ab has eigenvalues outside the stability degree region.
    
    [Z,At] = schur(a,'real');
    ev = ordeig(At); 
    if discr
        select = (abs(ev) <= smarg);
    else
        select = (real(ev) <= smarg);
    end
    [Z,At] = ordschur(Z,At,select);
    Bt = Z.'*b;  
    nb = sum(double(select == false));
    Et = eye(n);
    ninf = 0;
else
    if sepinf 
       mode = 2; qzopt = 3;
       [At,Et,Q,Z,nr,mr,minf] = sl_klf(a,e,tol,mode,qzopt); Q = Q';
       if ~isempty(nr) || ~isempty(mr)
           error('The pencil A-lambda*E is singular')
       end
       ninf = sum(minf);
       i1 =1:ninf; i2 = ninf+1:n;
       [At(i2,i2),Et(i2,i2),Q2,Z2] = qz(At(i2,i2),Et(i2,i2),'real');
       At(i1,i2) = At(i1,i2)*Z2; Et(i1,i2) = Et(i1,i2)*Z2;
       Q(i2,:) = Q2*Q(i2,:); Z(:,i2) = Z(:,i2)*Z2; 
    else
       ninf = 0;
       [At,Et,Q,Z] = qz(a,e,'real');
       i2 = ninf+1:n;
    end
    % compute orthogonal Q and Z such that
    %
    %      Q1*(Af-lambda Ef)*Z1 = [   Ag-lambda Eg      *       ]
    %                             [       0        Ab-lambda Eb ]
    %
    % where Ag-lambda Eg has eigenvalues within the stability degree region
    % and Ab-lambda Eb has eigenvalues outside the stability degree region.
    
    % separate stable and unstable parts
    ev = ordeig(At,Et); 
    if discr
        select = [ false(ninf,1); (abs(ev(i2)) <= smarg)];
    else
        select = [ false(ninf,1); ((real(ev(i2)) <= smarg) | (abs(ev(i2))*tola > 1))];
    end
    [At,Et,Q,Z] = ordqz(At,Et,Q,Z,select);
    Bt = Q*b; 
    nb = sum(double(select == false))-ninf;
end       

nfg = n-nb-ninf;
nfb = nb;
nrmBt = max(1,norm(Bt,1));
fnrmtol = 100*max(1,norm(At,1))/nrmBt;

nuf = 0;
nc = n;
f = zeros(m,n); 
while nb
   if nb == 1 || abs(At(nc,nc-1)) == 0
      k = 1;
   else
      k = 2;
   end
   kk = nc-k+1:nc;
   if standsys
      evb = ordeig(At(kk,kk));
   else
      evb = ordeig(At(kk,kk),Et(kk,kk));
   end
   select = false(nb,1); 
   if norm(Bt(kk,:)) < tolb
      % deflate uncontrollable stable block
      nb = nb-k; nc = nc-k; 
      nuf = nuf+1;
   else
      a2 = At(kk,kk); b2 = Bt(kk,:); 
      if standsys, e2 = eye(k); else, e2 = Et(kk,kk);  end
      if k == 1 
         % assign a single eigenvalue 
         [gamma,poles] = eigselect1(poles,sdeg,evb(end));
         if isempty(gamma)
             % no real pole available, adjoin a new 1x1 block if possible
             if nb == 1
                 % incompatible poles with the eigenvalue structure
                 % assign the last real pole to SDEG (if possible)
                 if isempty(sdeg)
                    warning('No real eigenvalue available for assignment')
                    gamma = evb(nb);
                 else
                    warning(['No real eigenvalue available for assignment: assigning instead SDEG = ', num2str(sdeg)])
                    gamma = sdeg;
                 end
                 f2 = -b2\(a2-e2*gamma);
             else
                 % adjoin a real block or interchange the last two blocks
                 k = 2; kk = nc-k+1:nc; 
                 if nc == 2 || abs(At(nc-1,nc-2)) == 0
                    % adjoin blocks and update evb 
                    if standsys
                       evb = max(ordeig(At(kk,kk)));
                    else
                       evb = max(ordeig(At(kk,kk),Et(kk,kk)));
                    end
                 else
                    % interchange last two blocks
                    i1 = 1:nc-3; i2 = nc-2:nc; lcol = nc+1:n;
                    if standsys
                      [Q3,At(i2,i2)] = ordschur(eye(3),At(i2,i2),[false,false,true]);
                      At(i1,i2) = At(i1,i2)*Q3; Bt(i2,:) = Q3.'*Bt(i2,:); 
                      At(i2,lcol) = Q3.'*At(i2,lcol); 
                      Z(:,i2) = Z(:,i2)*Q3; 
                    else
                      [At(i2,i2),Et(i2,i2),Q3,Z3] = ordqz(At(i2,i2),Et(i2,i2),eye(3),eye(3),[false,false,true]);
                      At(i1,i2) = At(i1,i2)*Z3; Et(i1,i2) = Et(i1,i2)*Z3;
                      Bt(i2,:) = Q3*Bt(i2,:); 
                      At(i2,lcol) = Q3*At(i2,lcol); Et(i2,lcol) = Q3*Et(i2,lcol); 
                      Q(i2,:) = Q3*Q(i2,:); Z(:,i2) = Z(:,i2)*Z3; 
                    end
                    if standsys
                       evb = ordeig(At(kk,kk));
                    else
                       evb = ordeig(At(kk,kk),Et(kk,kk));
                    end
                end
                 a2 = At(kk,kk); b2 = Bt(kk,:); 
                 if standsys, e2 = eye(k); else, e2 = Et(kk,kk);  end
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
            if ~standsys 
               Et(kk,kk) = u'*Et(kk,kk); Et(:,kk) = Et(:,kk)*v; 
%               evb(kk) = ordeig(At(kk,kk),Et(kk,kk));
%             else
%                evb(kk) = ordeig(At(kk,kk));
            end
            Bt(kk,:) = u'*Bt(kk,:);
            Q(kk,:) = u'*Q(kk,:);
         end
      end
      if isempty(f2)
         nb = nb-1; nc = nc-1; nuf = nuf+1; 
      else
         if norm(f2,inf) > fnrmtol
            %warning('Possible loss of numerical reliability due to high feedback gain')
         end
         % update matrices Acl and F
         At(:,kk) = At(:,kk)+Bt*f2; 
         f = f + f2*Z(:,kk)'; 
         % reorder eigenvalues 
         select(end-k+1:end) = true; i2 = nc-nb+1:nc; i1 = 1:nc-nb; 
         lcol = nc+1:n;
         if standsys
            if k == 2
               % standardize Schur form (needed for MATLAB 2020b)
               [Q2,At(kk,kk)] = schur(At(kk,kk));
                At(1:nc-2,kk) = At(1:nc-2,kk)*Q2; Bt(kk,:) = Q2.'*Bt(kk,:); 
                At(kk,lcol) = Q2.'*At(kk,lcol); 
                Z(:,kk) = Z(:,kk)*Q2; 
            end
            [Q3,At(i2,i2)] = ordschur(eye(nb),At(i2,i2),select);
            At(i1,i2) = At(i1,i2)*Q3; Bt(i2,:) = Q3.'*Bt(i2,:); 
            At(i2,lcol) = Q3.'*At(i2,lcol); 
            Z(:,i2) = Z(:,i2)*Q3; 
         else 
            if k == 2
               % standardize generalized Schur form (needed for MATLAB 2020b)
               [At(kk,kk),Et(kk,kk),Q2,Z2] = qz(At(kk,kk),Et(kk,kk),'real');
                At(1:nc-2,kk) = At(1:nc-2,kk)*Z2; Et(1:nc-2,kk) = Et(1:nc-2,kk)*Z2; 
                Bt(kk,:) = Q2*Bt(kk,:); 
                At(kk,lcol) = Q2*At(kk,lcol); Et(kk,lcol) = Q2*Et(kk,lcol);  
                Q(kk,:) = Q2*Q(kk,:); Z(:,kk) = Z(:,kk)*Z2; 
            elseif Et(kk,kk) < 0
                At(1:nc,kk) = -At(1:nc,kk);
                Et(1:nc,kk) = -Et(1:nc,kk);
                Z(:,kk) = -Z(:,kk);
            end
            [At(i2,i2),Et(i2,i2),Q3,Z3] = ...
                   ordqz(At(i2,i2),Et(i2,i2),eye(nb),eye(nb),select);
            At(i1,i2) = At(i1,i2)*Z3; Et(i1,i2) = Et(i1,i2)*Z3;
            Bt(i2,:) = Q3*Bt(i2,:); 
            At(i2,lcol) = Q3*At(i2,lcol); Et(i2,lcol) = Q3*Et(i2,lcol); 
            Q(i2,:) = Q3*Q(i2,:); Z(:,i2) = Z(:,i2)*Z3; 
         end
         nb = nb-k;
      end
   end
end
if nargout > 1
   info = struct('Acl',At,'Ecl',[],'Bcl',Bt,'Q',[],'Z',Z,...
                 'ninf',ninf,'nfg',nfg,'naf',nfb-nuf,'nuf',nuf);
   if standsys
      info.Q = Z.';
      info.Ecl = eye(n);
   else
      info.Q = Q; 
      info.Ecl = Et;
   end
end

% end GSFSTAB
end



