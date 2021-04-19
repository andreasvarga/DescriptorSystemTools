function [sysx,info,sysy] = grmcover2(sys1,sys2,tol)
% GRMCOVER2  Right minimum dynamic cover of Type 2 based order reduction. 
% [SYSX,INFO,SYSY] = GRMCOVER2(SYS1,SYS2,TOL) determines for the proper  
%   system SYS = [ SYS1 SYS2 ] with the descriptor system realization 
%   (A-lambda*E,[B1,B2],C,[D1,D2]) with E invertible, a state-feedback F 
%   and a feedforward  gain G such that the descriptor system realizations 
%   SYSX = (A+B2*F-lambda*E,B1+B2*G,C+D2*F,D1+D2*G) and 
%   SYSY = (A+B2*F-lambda*E,B1,F,G) have a maximum number of 
%   uncontrollable eigenvalues and the corresponding transfer-function 
%   matrices satisfy 
%       SYSX(lambda) = SYS1(lambda) + SYS2(lambda)*SYSY(lambda). 
%   A minimum dynamic cover of Type 2 is determined to obtain controllable
%   realizations of SYSX and SYSY of the form 
%   SYSX = (Ar-lambda*Er,Br,Cr1,Dr1) and SYSY = (Ar-lambda*Er,Br,Cr2,Dr2),
%   with the pair (Ar-lambda*Er,Br) in a controllability staircase form. 
%   TOL is a tolerance for rank determinations 
%   (default: TOL = 0, i.e., internally computed).
%
%   The resulting INFO structure contains additional information:
%   INFO.stdim is a vector which contains the dimensions of the diagonal 
%        blocks of Ar-lambda*Er, which are the row dimensions of the full 
%        row rank diagonal blocks of the pencil [Br Ar-lambda*Er] in 
%        controllability staircase form.
%   INFO.tcond is the maximum of the Frobenius-norm condition numbers of 
%        the employed non-orthogonal transformation matrices.
%   INFO.fnorm is the Frobenius-norm of the employed state-feedback F. 
%   INFO.gnorm is the Frobenius-norm of the employed feedforward gain G. 
%   Note: Large values of INFO.tcond, INFO.fnorm or INFO.gnorm indicate 
%         possible loss of numerical accuracy. 
%
% [SYSX,INFO,SYSY] = GRMCOVER2(SYS,M1,TOL) uses the compound realization 
%     of SYS = [SYS1 SYS2], where SYS1 has M1 inputs.
%
% Note: GRMCOVER2 also works if the original descriptor system realization
%       has E singular, but SYS2 is proper. In this case, the order 
%       reduction is performed working only with the proper part of SYS1. 
%       The improper part of SYS1 is included without modification in the  
%       resulting realization of SYSX. 
%
% See also GRMCOVER1. 

%  Author:      A. Varga, 25.01.2016.
%  Revision(s): A. Varga, 31.01.2016, 30.10.2016, 09.05.2017, 27.07.2019.
%
%  Method: The state feedback gain F and the feedforward gain G are 
%  determined using the method  of [1] to compute Type 2 minimum dynamic 
%  covers for standard systems (E = I) and the method of [2] for proper  
%  descriptor systems. 
%
%  Note: The resulting McMillan degree of SYSX is the least achievable one
%  provided the realization (A-lambda*E,B2,C,D2) is maximally observable 
%  (i.e., the pair (A+B2*F-lambda*E,C+D2*F) is observable for any F). 
%
%  References:
%  [1] A. Varga, 
%      Reliable algorithms for computing minimal dynamic covers,
%      Proc. CDC'03, Maui, Hawaii, 2003.
%  [2] A. Varga. 
%      Reliable algorithms for computing minimal dynamic covers for 
%      descriptor systems.
%      Proc. MTNS Symposium, Leuven, Belgium, 2004. 

narginchk(2,3)
nargoutchk(0,3)

if nargin < 3
    tol = 0;
else
   validateattributes(tol, {'double'}, {'real','scalar', '>=', 0},'','TOL') 
end

if ~isa(sys1,'ss')
   error('The input system SYS1 must be a state-space system object')
end

[p,m] = size(sys1);
if ~isa(sys2,'ss')
   if ~isa(sys2,'double') 
      error('SYS2 must be an SS object or a positive integer')
   else
      validateattributes(sys2,{'double'},{'integer','scalar','>=',0,'<=',m},'','M1')
   end
end  

Ts = sys1.Ts;
discr = (Ts ~= 0); 

if isa(sys2,'double') 
    m1 = sys2; 
    m2 = m-m1;
    [A,B,C,D,E,Ts] = dssdata(sys1);
else
   if ~isa(sys2,'ss')
      error('The input system SYS2 must be a state-space system object')
   end
   if p ~= size(sys2,1)
      error('The systems SYS1 and SYS2 must have the same number of outputs')
   end
   if (discr && (sys2.Ts <= 0)) || (~discr && (sys2.Ts > 0))
      error('The systems SYS1 and SYS2 must have the same type')
   end
   if discr && (Ts ~= sys2.Ts) 
      error('The systems SYS1 and SYS2 must have the same sampling period')
   end
   m1 = m; 
   m2 = size(sys2,2);
   [A,B,C,D,E] = dssdata(gir([sys1,sys2],tol));
   m = m1+m2; 
end    


% exchange the roles of B1 and B2
nsys1 = size(A,1); 
isys1 = m1+1:m; isys2 = 1:m1; 
sstype = isequal(E,eye(nsys1));

% handle symple cases
if nsys1 == 0
    sysx = sys1(:,isys2);
    info.stdim = [];
    info.tcond = 1;
    info.fnorm = 0;
    info.gnorm = 0;
    if nargout == 3
        sysy = ss(zeros(m2,m1)); 
        sysy.Ts = Ts;
    end
    return
end

if  m1 == 0
    sysx = ss(zeros(p,m1)); 
    info.stdim = [];
    info.tcond = 1;
    info.fnorm = 0;
    info.gnorm = 0;
    if nargout == 3
        sysy = ss(zeros(m2,m1)); 
        sysy.Ts = Ts;
    end
    return
end

% reduce to the special controllability form of [1]
[At,Et,Bt,Ct,ncont,tau] = sl_gstra(10,A,E,[B(:,isys1), B(:,isys2)],C,m2,tol); 

if ncont == 0
    sysx = ss(D(:,isys2));
    sysx.Ts = Ts;
    info.stdim = [];
    info.tcond = 1;
    info.fnorm = 0;
    if nargout == 3
       sysy = ss(zeros(m2,m1)); 
       sysy.Ts = Ts;
    end
    return
end

if m2 == 0
   if sstype
      sysx = ss(At,Bt,Ct,D,Ts);
   else
      sysx = dss(At,Bt,Ct,D,Et,Ts);
   end
   info.stdim = tau(tau > 0);
   info.tcond = 1;
   info.fnorm = 0;
   info.gnorm = 0;
   if nargout == 3
      sysy = ss(zeros(m2,m1)); 
      sysy.Ts = Ts;
   end
   return
end

% handle improper case (more efficient implementation is possible)
if rcond(Et) < eps
   [sysf,sysi] = gsdec(dss(A,B,C,D,E,Ts),struct('tol',tol,'job','finite'));
   % check properness of SYS2
   sys2i = gss2ss(gir(sysi(:,m1+1:end),tol,'infinite'),tol);
   if size(sys2i,'order')
       error('The system SYS2 must be proper')
   else
       % add constant term in the finite part (without increasing order)
       sysf.d = [sysf(:,1:m1).d sysf(:,m1+1:end).d+sys2i.d];
   end
   [sysxf,info,sysy] = grmcover2(sysf,m1,tol); 
   sysx = sysxf + sysi(:,1:m1); 
   return
end

% perform permutation for type II cover
cind = length(tau)/2; 
ind = 1:ncont;
i1 = 1:2:length(tau);
i2 = i1+1;
ind1 = tau(i1);
ind2 = tau(i2);
n1 = sum(ind1);
n2 = sum(ind2);
i2 = 1; i1 = 1; ioff1 = 0; ioff2 = 0; 
for i = 1:cind
   ioff1 =  ioff1 + ind1(i);
   ind(i2:i2+ind2(i)-1) =  (i2:i2+ind2(i)-1)+ioff1; 
   i2 = i2+ind2(i);
   ind((i1:i1+ind1(i)-1)+n2) =  (i1:i1+ind1(i)-1)+ioff2; 
   ioff2 =  ioff2 + ind2(i);
   i1 = i1+ind1(i);
end
At = At(ind,ind); Bt = Bt(ind,:); Ct = Ct(:,ind); 

if ~sstype
   Et = Et(ind,ind); 
   % check invertibility of the leading part of Et 
   i1 = 1:n2;
   if ~isempty(i1) && (rcond(Et(i1,i1)) < eps || ...
      norm(Et(i1,i1),1) < numel(Et)*eps(norm(Et,1)))
      error('The system [SYS1 SYS2] is possibly non-proper')
   end
end

% anihilate lower left blocks of At and Et
tcond = ncont; % the Frobenius condition number of identity matrix
nlowc1 = n2-ind2(cind);  nlowc2 = n1+n2-ind1(cind);
nlowr1 = n1+n2;          nlowr2 = nlowr1-ind1(cind);
for i=cind:-1:1
   if i == 1
      ice1 = 1:n2; ir1 = n2+1:n2+ind1(i);
   else
      nlowc1 = nlowc1-ind2(i-1); nlowc2 = nlowc2-ind1(i-1); 
      nlowr1 = nlowr1-ind1(i); nlowr2 = nlowr2-ind1(i-1);
      ic1 = nlowc1+1:n2; ic2 = nlowc2+1:nlowc2+ind1(i-1);
      ir1 =  nlowr1+1:nlowr1+ind1(i); 
      ice1 = nlowc1+ind2(i-1)+1:n2; 
   end
   if ~sstype
      y = -Et(ir1,ice1)/Et(ice1,ice1); 
      tcond = max(tcond,ncont+norm(y,'fro')^2);
      Et(ir1,:) = Et(ir1,:)+y*Et(ice1,:);
      At(ir1,:) = At(ir1,:)+y*At(ice1,:);
      Bt(ir1,:) = Bt(ir1,:)+y*Bt(ice1,:);
   end
   if i == 1, break, end
   x = -At(ir1,ic2)\At(ir1,ic1); 
   tcond = max(tcond,ncont+norm(x,'fro')^2);
   ir = 1:nlowr1;
   At(ir,ic1) = At(ir,ic2)*x+At(ir,ic1); 
   At(ir1,ic1) = zeros(ind1(i),n2-nlowc1);
   Ct(:,ic1) = Ct(:,ic2)*x+Ct(:,ic1); 
   if sstype
      x = -x; At(ic2,:) = At(ic2,:)+x*At(ic1,:);
      Bt(ic2,:) = Bt(ic2,:)+x*Bt(ic1,:);
   else
      Et(ir,ic1) = Et(ir,ic2)*x+Et(ir,ic1); 
   end
end


% form the reduced system
ic = 1:n2; i3 = n2+1:n2+ind1(1); 
f2 = -Bt(i3,1:m2)\At(i3,ic);
g = -Bt(i3,1:m2)\Bt(i3,m2+1:m); 
info.fnorm = norm(f2,'fro');
info.gnorm = norm(g,'fro');
info.tcond = tcond; 
info.stdim = ind2(ind2 > 0)';
   
if sstype
   sysx = ss(At(ic,ic),Bt(ic,m2+1:m),Ct(:,ic)+D(:,isys1)*f2,D(:,isys1)*g+D(:,isys2),Ts);
else
   sysx = dss(At(ic,ic),Bt(ic,m2+1:m),Ct(:,ic)+D(:,isys1)*f2,D(:,isys1)*g+D(:,isys2),Et(ic,ic),Ts);
end
if nargout == 3
   if sstype
      sysy = ss(At(ic,ic),Bt(ic,m2+1:m),f2,g,Ts);
   else
      sysy = dss(At(ic,ic),Bt(ic,m2+1:m),f2,g,Et(ic,ic),Ts);
   end
end

% end GRMCOVER2
end

