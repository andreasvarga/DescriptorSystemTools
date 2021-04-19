function [At,Et,info,Q,Z] = gklf(A,E,tol,varargin)
%GKLF Kronecker-like staircase forms of a linear matrix pencil.
% [At,Et,INFO,Q,Z] = GKLF(A,E,TOL) computes for a given matrix pencil 
% A-s*E the Kronecker-like staircase form 
%
%                           ( Br Ar-s*Er     *           *         *     )
%                           ( 0    0     Ainf-s*Einf     *         *     )
%  At-s*Et = Q'*(A-s*E)*Z = ( 0    0         0       Afin-s*Efin   *     ),
%                           ( 0    0         0           0       Al-s*El )
%                           ( 0    0         0           0          Cl   )
%
% where Q and Z are orthogonal transformation matrices. TOL is an optional 
% tolerance to be used for rank tests. The INFO structure provides
% information on the structure of the pencil At-s*Et as follows:
%   [ Br Ar-s*Er ], with Er invertible and upper triangular, contains the  
%      right Kronecker structure and is in a controllability staircase form 
%      with INFO.mr(i) x INFO.nr(i) full row rank diagonal blocks;
%      INFO.mr and INFO.nr are empty if no right structure exists; 
%   Ainf-s*Einf, with Ainf invertible and Einf nilpotent, contains the 
%     infinite structure and is in a block upper triangular form with 
%     INFO.minf(i) x INFO.minf(i) upper triangular diagonal blocks;
%     INFO.minf is empty if no infinite structure exists; 
%   Afin-s*Efin contains the finite structure and is INFO.mf x INFO.mf;
%   [ Al-s*El; Cl ], with El invertible and upper triangular, contains the 
%     left Kronecker structure and is in an observability staircase form 
%     with INFO.ml(i) x INFO.nl(i) full column rank subdiagonal blocks;
%     INFO.ml and INFO.nl are empty if no left structure exists. 
%
% [At,Et,INFO,Q,Z] = GKLF(A,E,TOL,'reverse') computes the Kronecker-like 
% staircase form with reverse ordering of the blocks of the regular part
%
%                           ( Br Ar-s*Er    *            *           *    )
%                           ( 0    0    Afin-s*Efin      *           *    ) 
%  At-s*Et = Q'*(A-s*E)*Z = ( 0    0        0        Ainf-s*Einf     *    )
%                           ( 0    0        0            0        Al-s*El )
%                           ( 0    0        0            0           Cl   )
%
%
% [At,Et,INFO,Q,Z] = GKLF(A,E,TOL,'right') computes the Kronecker-like 
% staircase form exhibiting the right and infinite Kronecker structures
%
%                           ( Br Ar-s*Er     *             *      )
%  At-s*Et = Q'*(A-s*E)*Z = ( 0    0     Ainf-s*Einf       *      ) ,
%                           ( 0    0        0         Af,l-s*Ef,l )
%
% where Af,l-s*Ef,l contains the finite and left Kronecker structures. 
% In this case, INFO.ml and INFO.nl are set to the negative values of the
% row and column dimensions of the pencil Af,l-s*Ef,l, respectively, and
% INFO.mf = 0. 
%
% [At,Et,INFO,Q,Z] = GKLF(A,E,TOL,'left') computes the Kronecker-like 
% staircase form exhibiting the left and infinite Kronecker structures
%
%                           ( Ar,f-s*Er,f     *          *    )
%  At-s*Et = Q'*(A-s*E)*Z = (     0       Ainf-s*Einf    *    ) ,
%                           (     0           0       Al-s*El )
%                           (     0           0          Cl   )
%
% where Ar,f-sEr,f contains the right and finite Kronecker structures. 
% In this case, INFO.mr and INFO.nr are set to the negative values of the
% row and column dimensions of the pencil Ar,f-s*Er,f, respectively, and
% INFO.mf = 0.  
%
% [At,Et,INFO,Q,Z] = GKLF(A,E,TOL,...,'noQ') computes the respective
% Kronecker-like staircase forms without internally accumulating the  
% transformation matrix Q. In this case, Q = [] (the empty matrix). 
%
% See also SL_KLF.

% Author:      A. Varga, 24.12.2015.
% Revision(s): A. Varga, 26.10.2016, 11.06.2017. 
%
% References
% [1] Beelen, Th. and Van Dooren, P.
%     An improved algorithm for the computation of Kronecker's
%     canonical form of a singular pencil.
%     Linear Algebra and Applications, vol. 105, pp. 9-65, 1988.

narginchk(2,5)
nargoutchk(0,5)

if nargin < 3
   tol = 0; 
else
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','TOL') 
end

withQ = isempty(find(strcmp(varargin,'noQ'),1)) && (nargout > 3);
withZ = (nargout > 4); 

if withQ 
   if withZ, qzopt = 3; else qzopt = 1; end
else
   if withZ, qzopt = 2; else qzopt = 0; end
end

if nargin < 4
   job = 0; %opt = 'normal';
else
   if ~isempty(find(strcmp(varargin,'right'),1)); 
      job = 2; % opt = 'right'
   elseif ~isempty(find(strcmp(varargin,'left'), 1)); 
      job = 3; % opt = 'left'
   elseif ~isempty(find(strcmp(varargin,'reverse'), 1)); 
      job = 1; % opt = 'reverse';
   else
      job = 0; % opt = 'normal'
   end
end

[m,n] = size(A); 
[me,ne] = size(E); 
if m ~= me || n ~= ne
    error('A and E must have the same dimensions')
end

switch job
  case 0
    % reduce A-sE to the Kronecker-like form
    %                ( Br Ar-sEr      *           *      )
    %  Q'*(A-sE)*Z = ( 0    0    Ainf-sEinf       *      ) ,
    %                ( 0    0        0        Af,l-sEf,l )
    % where [Br Ar-sEr] and Ainf-sEinf) contain the right and infinite 
    % Kronecker structures and Af,l-sEf,l contains the finite and left 
    % Kronecker structure. 
    %
    mode = 2; 
    [At,Et,Q,Z,nr,mr,minf] = sl_klf(A,E,tol,mode,qzopt); 

    % compute the dimensions of A(inf) and of the right and infinite parts
    di = sum(minf); mrinf = sum(mr)+di; nrinf = sum(nr)+di; 
    % compute the dimensions of the finite and left part
    mfl = m-mrinf; nfl = n-nrinf; 

    if nfl == 0
       % case of no finite and left parts 
       mf = 0; 
       if mfl
          ml = ones(1,mfl); nl = zeros(1,mfl);
       else
          ml = []; nl = [];
       end
    elseif nfl == mfl
       % case of no left part
       mf = nfl; ml = []; nl = [];
    else
       % separate finite and left parts by using the pertransposed 
       % pencil (Af,l-sEf,l)^P
       ifl = mrinf+1:m; jfl = nrinf+1:n;
       Afl = flip(flip(At(ifl,jfl)),2); Efl = flip(flip(Et(ifl,jfl)),2);
       mode = 1; 
       if withQ, qzopt0 = 3; else qzopt0 = 1; end
       [A2,B2,Z2,Q2,ml,nl] = sl_klf(Afl',Efl',tol,mode,qzopt0);
       At(ifl,jfl) = flip(flip(A2,1),2)'; Et(ifl,jfl) = flip(flip(B2,1),2)';
       ml = flip(ml); nl = flip(nl); 
       % compute dimension of finite part
       mf = mfl-sum(ml);
       % update the rest of submatrices and transformations 
       Z2 = flip(flip(Z2,1),2);
       irinf = 1:mrinf; 
       At(irinf,jfl) = At(irinf,jfl)*Z2; Et(irinf,jfl) = Et(irinf,jfl)*Z2;
       if withQ
           Q2 = flip(flip(Q2,1),2);
           Q(:,ifl)=Q(:,ifl)*Q2; 
       end
       if withZ
          Z(:,jfl)=Z(:,jfl)*Z2; 
       end
    end
  case 1
    % apply GKLF to the transposed pair (A',E') and pertranspose the
    % results   
    if ~withZ, optZ = 'noQ'; else optZ = ''; end
    %[At,Et,Z,Q,nl,ml,minf,mf,nr,mr] = gklf(A',E',tol,optQ,optZ);
    if withQ
       [At,Et,info,Z,Q] = gklf(A',E',tol,optZ);
    else
       [At,Et,info,Z] = gklf(A',E',tol,optZ);
    end
    nl = info.mr; ml = info.nr; 
    nr = info.ml; mr = info.nl;
    At = flip(flip(At,1),2)'; Et = flip(flip(Et,1),2)'; 
    if withQ, Q = flip(Q,2); end
    if withZ, Z = flip(Z,2); end
    info.ml = flip(ml); info.nl = flip(nl); 
    info.mr = flip(mr); info.nr = flip(nr); 
    info.minf = flip(info.minf); 
    return
  case 2  
    % reduce A-sE to the Kronecker-like form
    %                ( Br Ar-sEr      *           *      )
    %  Q'*(A-sE)*Z = ( 0    0    Ainf-sEinf       *      ) ,
    %                ( 0    0        0        Af,l-sEf,l )
    % where [Br Ar-sEr] and Ainf-sEinf) contain the right and infinite 
    % Kronecker structures and Af,l-sEf,l contains the finite and left 
    % Kronecker structure. 
    mode = 2; 
    [At,Et,Q,Z,nr,mr,minf] = sl_klf(A,E,tol,mode,qzopt); 
    % compute the negatives of the dimensions of the left and regular parts
    di = sum(minf); ml = di+sum(mr)-m; nl = di+sum(nr)-n; mf = 0;
  case 3
    % reduce A-sE to the Kronecker-like form
    %                ( Ar,f-sEr,f      *           *    )
    %  Q'*(A-sE)*Z = (     0       Ainf-sEinf      *    ) ,
    %                (     0           0         Al-sEl ) 
    %                (     0           0           Cl   )
    %
    % where Ar,f-sEr,f contains the right and finite Kronecker structures.
    %
    % apply SL_KLF to the transposed pair (A',E') and pertranspose the
    % results
    mode = 2; qzopt3 = qzopt;
    if qzopt == 2, qzopt3 = 1; elseif qzopt == 1,  qzopt3 = 2; end
    [At,Et,Z,Q,ml,nl,minf] = sl_klf(A',E',tol,mode,qzopt3); 
    At = flip(flip(At,1),2)'; Et = flip(flip(Et,1),2)'; 
    Q = flip(Q,2); Z = flip(Z,2);
    ml = flip(ml); nl = flip(nl); minf = flip(minf); 
    % compute the negatives of the dimensions of the right and regular parts
    di = sum(minf); mr = di+sum(ml)-m; nr = di+sum(nl)-n;  mf = 0;      
end
info.mr = mr'; info.nr = nr'; info.minf = minf'; info.mf = mf;
info.ml = ml'; info.nl = nl';

% end GKLF
end