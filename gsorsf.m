 function [at,et,q,z,dims,ni] = gsorsf(a,e,options)
%GSORSF Specially ordered generalized real Schur form.
%   [At,Et,Q,Z,DIMS,NI] = GSORSF(A,E,OPTIONS) calculates, for a given
%   matrix pair (A,E), the transformed pair (At,Et) and the orthogonal
%   transformation matrices Q and Z such that
%           At = Q'*A*Z,   Et = Q'*E*Z,
%   and At and Et are in a specially ordered generalized real Schur
%   form with the following block-upper triangular structure
%
%         [ Ai  *    *    *  ]          [ 0  *    *    *  ]
%    At = [  0  Ag   *    *  ]  ,  Et = [ 0  Eg   *    *  ]  .       (*)
%         [  0  0   Abf   *  ]          [ 0  0   Ebf   *  ]
%         [  0  0    0   Abi ]          [ 0  0    0   Ebi ]
%
%    For a regular pencil A-lambda*E, Ai is nonsingular and upper
%    triangular, the pair (Ag,Eg) has only finite generalized eigenvalues
%    in a "good" region Cg of the complex plane, the pair (Abf,Ebf)
%    has only finite generalized eigenvalues outside of Cg, and the pair
%    (Abi,Ebi) has only infinite generalized eigenvalues. By default,
%    the "good" region Cg is the open left half complex plane.
%    The OPTIONS structure allows to specify various user options,
%    as follows:
%    OPTIONS.tol   - specifies the tolerance for the singular values
%                    based rank determination of E
%                    (Default: tol = prod(size(E)) * eps(norm(E,1)))
%    OPTIONS.disc  - specifies the disc option for the "good" region Cg
%                    true  - Cg is a disc centered in the origin
%                    false - Cg is a left half complex plane (default)
%    OPTIONS.smarg - specifies the desired stability margin for the
%                    "good" generalized eigenvalues of the pair (A,E)
%                    (Default: -sqrt(eps)  if OPTIONS.disc = false, and
%                              1-sqrt(eps) if OPTIONS.disc = true)
%    OPTIONS.reverse - specifies the option for reverse ordering
%                    false - the diagonal blocks are as in (*)(default)
%                    true  - the alternative ordering is computed
%
%                         [ Abi  *   *    *  ]         [ Ebi  *   *    *  ]
%                    At = [  0  Abf  *    *  ],   Et = [  0  Ebf  *    *  ]
%                         [  0   0   Ag   *  ]         [  0   0   Eg   *  ]
%                         [  0   0   0    Ai ]         [  0   0   0    0  ]
%    OPTIONS.sepinf - specifies the option for separation of higher order
%                    infinite eigenvalues 
%                    true  - separate higher order generalized infinite 
%                            eigenvalues in the trailing positions if 
%                            OPTIONS.reverse = false, or in the leading
%                            positions if OPTIONS.reverse = true (default);
%                    false - no separation of higher order infinite
%                            eigenvalues. 
%    OPTIONS.fast  - specifies the option for fast separation of higher 
%                    order infinite eigenvalues (to be used in conjuntion 
%                    with the OPTIONS.sepinf = true), as follows:  
%                    true  - fast separation of higher order infinite 
%                            eigenvalues using orthogonal pencil
%                            manipulation techniques with QR decompositon
%                            based rank determinations (default)
%                    false - separation of higher order infinite
%                            eigenvalues, by using singular value 
%                            decomposition based rank determinations 
%                            (potentially more reliable, but slower).
%
%    DIMS(1:4) contain the orders of the diagonal blocks (Ai,Ag,Abf,Abi),
%    if OPTIONS.reverse = false, or of the diagonal blocks (Abi,Abf,Ag,Ai),
%    if OPTIONS.reverse = true.
%    NI(i) contains, in the case when OPTIONS.sepinf = true, the dimension 
%    of the i-th nilpotent diagonal block of the matrix Eb in the 
%    leading positions, if OPTIONS.reverse = true, or trailing positions,
%    if OPTIONS.reverse = false. Otherwise, NI is an empty matrix.
%
%    WARNING: For a non-regular pencil A-lambda*E, the computed results 
%             may be erroneous!
%

%    Author:  A. Varga 26.09.2016.
%    Revised: A. Varga 31.10.2016, 13.06.2017.
%
%  Method:  The Procedure GSORSF from [1] is implemented.  
%
%  References:
%  [1] Varga A.
%      On recursive computation of coprime factorizations of rational 
%      matrices. arXiv:1703.07307, https://arxiv.org/abs/1703.07307, 2017.

narginchk(1,3)
nargoutchk(0,6)

validateattributes(a, {'double'}, {'real','finite', 'square'});

n = size(a,1); 

if nargin == 1
    e = [];
else
    validateattributes(e, {'double'}, {'real','finite', 'square','size',[n,n]});
end

if nargin <= 2
   options = struct('tol',0);
end

standev = isempty(e) || isequal(e,eye(n));  

% decode options and set default values

% tolerance for rank determination
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end
if tol == 0 && ~standev
   % set absolute tolerance for determination of rank of E 
   tol = n*n*eps(norm(e,1));
end


% ordering option with respect to a disc
if isfield(options,'disc')
   discr = options.disc;
   validateattributes(discr, {'logical'},{'binary'},'','OPTIONS.disc') 
else
   discr = false;
end
if discr
    smax = 1-sqrt(eps); smin = 0;
else
    smax = -sqrt(eps);  smin = -inf;
end


% stability margin
if isfield(options,'smarg')
   smarg = options.smarg;
   if ~isempty(smarg)
%      validateattributes(smarg, {'double'},{'real','scalar','<=',smax,'>=',smin},'','OPTIONS.smarg') 
      validateattributes(smarg, {'double'},{'real','scalar','>=',smin},'','OPTIONS.smarg') 
   end
else
   smarg = [];
end
if isempty(smarg)
   % set stability margin
   smarg = smax;
end

% reverse ordering option
if isfield(options,'reverse')
   reverse = options.reverse;
   validateattributes(reverse, {'logical'},{'binary'},'','OPTIONS.reverse') 
else
   reverse = false;
end

% infinite eigenvalue separation
if isfield(options,'sepinf')
   sepinf = options.sepinf;
   validateattributes(sepinf, {'logical'},{'binary'},'','OPTIONS.sepinf') 
else
   sepinf = true;
end

% fast separation using the reduction to a Kronecker-like form 
if isfield(options,'fast')
   fast = options.fast;
   validateattributes(fast, {'logical'},{'binary'},'','OPTIONS.fast') 
else
   fast = true;
end

% separate stable and unstable parts
if  standev
    % compute orthogonal Q  such that
    %
    %      Q^T*A*Q = [ Ag   * ]
    %                [  0  Ab ]
    %
    % where Ag has eigenvalues within the stability margin region
    % and Ab has eigenvalues outside the stability margin region.
    
    [q,at] = schur(a,'real');
    ev = ordeig(at); 
    if discr
        select = abs(ev)<=smarg;
    else
        select = real(ev)<=smarg;
    end
    nb = sum(~select);  
    ng = n-nb;
    nsinf = 0;
    nbi = 0;
    if reverse
        select = ~select;
    end
    [q,at] = ordschur(q,at,select);
    et = e; z = q; ni = [];
else
    % compute orthogonal Q and Z such that
    %
    %                           [ A11      *             *       ]
    %      Q^T*(A-lambda E)*Z = [  0   Ag-lambda Eg      *       ]
    %                           [  0       0        Ab-lambda Eb ]
    %
    % where Ag-lambda Eg has eigenvalues within the stability margin region
    % and Ab-lambda Eb has eigenvalues outside the stability margin region.
    
    % separate the first order infinite eigenvalues
    %        At = [ Ai   *  ] ,  Et =  [ 0  *   ]  if REVERSE = false;
    %             [  0  Agb ]          [ 0  Egb ] 
    %
    %        At = [ Agb  *  ] ,  Et =  [ Egb  * ]  if REVERSE = true.
    %             [  0   Ai ]          [ 0    0 ]  
    [at,et,q,z,nsinf] = gsinf0(a,e,tol,reverse);
    if nsinf 
       nbi = 0;  ni = [];
       if sepinf 
          % separate higher order infinite eigenvalues of (Agb,Egb) in
          % trailing positions, if REVERSE = false, and in the leading
          % positions, if REVERSE = true.
          nsinf1 = nsinf;  
          if reverse
             n1 = n-nsinf; i3 = n1+1:n;  
             if fast
                i2 = 1:n1;  
                [a1,e1,info1,q1,z1] = gklf(at(i2,i2),et(i2,i2),tol);                       
                ni = info1.minf;  
                nbi = sum(ni);
                if nbi
                   at(i2,i2) = a1; et(i2,i2) = e1; 
                   at(i2,i3) = q1.'*at(i2,i3); et(i2,i3) = q1.'*et(i2,i3); 
                   q(:,i2) = q(:,i2)*q1; z(:,i2) = z(:,i2)*z1;
                end
             else
                while nsinf1
                  i2 = nbi+1:n1; i1 = 1:nbi; 
                  [at(i2,i2),et(i2,i2),q1,z1,nsinf1] = ...
                                gsinf0(at(i2,i2),et(i2,i2),tol,false);  
                  if nsinf1, ni = [ni,nsinf1]; end
                  nbi = nbi+nsinf1; 
                  at(i2,i3) = q1.'*at(i2,i3); et(i2,i3) = q1.'*et(i2,i3); 
                  at(i1,i2) = at(i1,i2)*z1; et(i1,i2) = et(i1,i2)*z1; 
                  q(:,i2) = q(:,i2)*q1; z(:,i2) = z(:,i2)*z1;
                end
             end
          else
             if fast
                i2 = nsinf+1:n; i1 = 1:nsinf; 
                [a1,e1,info1,q1,z1] = ...
                                gklf(at(i2,i2),et(i2,i2),tol,'reverse');                       
                ni = info1.minf;  
                nbi = sum(ni);
                if nbi
                   at(i2,i2) = a1; et(i2,i2) = e1; 
                   at(i1,i2) = at(i1,i2)*z1; et(i1,i2) = et(i1,i2)*z1; 
                   q(:,i2) = q(:,i2)*q1; z(:,i2) = z(:,i2)*z1;
                end
             else
                i1 = 1:nsinf; n1 = n; 
                while nsinf1
                  i2 = nsinf+1:n1; i3 = n1+1:n; 
                  [at(i2,i2),et(i2,i2),q1,z1,nsinf1] = ...
                                  gsinf0(at(i2,i2),et(i2,i2),tol,true);  
                  if nsinf1, ni = [nsinf1,ni]; end
                  nbi = nbi+nsinf1; n1 = n1-nsinf1;
                  at(i2,i3) = q1.'*at(i2,i3); et(i2,i3) = q1.'*et(i2,i3); 
                  at(i1,i2) = at(i1,i2)*z1; et(i1,i2) = et(i1,i2)*z1;        
                  q(:,i2) = q(:,i2)*q1; z(:,i2) = z(:,i2)*z1; 
                end
             end
          end
       end
       if reverse
          re = n-nsinf; 
          i1 = 1:nbi; i2 = nbi+1:re; i3 = re+1:n;
       else
          i1 = 1:nsinf; i2 = nsinf+1:n-nbi; i3 = n-nbi+1:n;
       end
       % separate stable and unstable parts of (Agb,Egb)
       [at(i2,i2),et(i2,i2),qt,zt] = qz(at(i2,i2),et(i2,i2),'real');
       if ( rcond(at(i1,i1)-rand*et(i1,i1)) < sqrt(eps) || ...
            rcond(at(i2,i2)-rand*et(i2,i2)) < sqrt(eps) ) && ...
            rcond(a-rand*e) < sqrt(eps)
            warning('A-lambda*E likely non-regular: results may be erroneous.')
       end
       ev = ordeig(at(i2,i2),et(i2,i2)); 
       if discr
          select = abs(ev) <= smarg;
       else
          select = (real(ev) <= smarg) & (abs(ev)*tol < 1);
       end
       nb = nbi + sum(~select);  
       ng = n-nb-nsinf;
       % it should be
       % nb = sum(~select);  
       % ng = n-nb-nsinf-nbi;
       
       if reverse
          select = ~select;
       end
       [at(i2,i2),et(i2,i2),qt,zt] = ordqz(at(i2,i2),et(i2,i2),qt,zt,select);
       if reverse
          at(i2,i3) = qt*at(i2,i3); et(i2,i3) = qt*et(i2,i3);
          if ~isempty(i1)
              at(i1,i2) = at(i1,i2)*zt; et(i1,i2) = et(i1,i2)*zt; 
          end
       else
          at(i2,i3) = qt*at(i2,i3); et(i2,i3) = qt*et(i2,i3);
          if ~isempty(i1)
              at(i1,i2) = at(i1,i2)*zt; et(i1,i2) = et(i1,i2)*zt;
          end
       end
       q(:,i2) = q(:,i2)*qt.'; z(:,i2) = z(:,i2)*zt;
   else
       % separate stable and unstable parts of (Agb,Egb) 
       [at,et,q,z] = qz(a,e,'real');
       ev = ordeig(at,et); 
       if discr
          select = abs(ev) <= smarg;
       else
          select = (real(ev) <= smarg) & (abs(ev)*tol < 1);
       end
       nb = sum(double(select == false));
       ng = n-nb-nsinf;
       ni = [];
       if reverse
          select = ~select;
       end
       [at,et,q,z] = ordqz(at,et,q,z,select); q = q.';
    end    
end        

if reverse
   dims = [sum(ni) nb ng nsinf];
   % dims = [nbi nb ng nsinf];  % this is a correct alternative
else
   dims = [nsinf ng nb sum(ni)];
   % dims = [nsinf ng nb nbi];  % this is a correct alternative
end

% end GSORCF
end

function [at,et,q,z,nsinf] = gsinf0(a,e,tol,reverse)
% GSINF0 Separation of the first order infinite eigenvalues
%        [At,Et,Q,Z,NSINF] = GSINF0(A,E,TOL,REVERSE) determines At and Et, 
%        and the transformation matrices Q and Z such that
%               At = Q'*A*Z,   Et = Q'*E*Z,
%        At and Et have the following block-upper triangular structure:
%        
%        At = [ Ai   *  ] ,  Et =  [ 0  *   ]  if REVERSE = false;
%             [  0  Agb ]          [ 0  Egb ] 
%
%        At = [ Agb  *  ] ,  Et =  [ Egb  * ]  if REVERSE = true.
%             [  0   Ai ]          [ 0    0 ]  
%
%        The NSINF x NSINF matrix Ai is non-singular and upper triangular.
%        TOL is a tolerance for the determination of the rank of E.

[q,et,z] = svd(e);  n = size(a,1);
re = sum(diag(et) > tol); nsinf = n-re;
if re < n
   if reverse
      at = q.'*a*z;
      i2 = 1:re; i1 = re+1:n;
      [q2,r] = qr(at(i1,:).'); 
      z2 = flip(q2,2);
      at = [at(i2,:)*z2; flip(flip(r.',2),1)]; et(i2,:) = et(i2,:)*z2;
      z = z*z2; q(:,i1) = flip(q(:,i1),2); et(i1,i1) = 0;
   else
      perm = [re+1:n 1:re ]; q = q(:,perm); z = z(:,perm);
      et = et(perm,perm); at = q.'*a*z; 
      i1 = 1:nsinf; i2 = nsinf+1:n;
      [q2,r] = qr(at(:,i1)); 
      at = [ r q2.'*at(:,i2)]; et(:,i2) = q2.'*et(:,i2);
      q = q*q2; et(i1,i1) = 0;
   end
else
   at = a; et = e; q = eye(n); z = eye(n); 
end
end
