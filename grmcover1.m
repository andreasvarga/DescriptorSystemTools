function [sysx,info,sysy] = grmcover1(sys1,sys2,tol)
% GRMCOVER1  Right minimum dynamic cover of Type 1 based order reduction. 
% [SYSX,INFO,SYSY] = GRMCOVER1(SYS1,SYS2,TOL) determines for the proper  
%   system SYS = [ SYS1 SYS2 ] with the descriptor system realization 
%   (A-lambda*E,[B1,B2],C,[D1,D2]) with E invertible, a state-feedback F 
%   such that the descriptor system realizations 
%   SYSX = (A+B2*F-lambda*E,B1,C+D2*F,D1) and 
%   SYSY = (A+B2*F-lambda*E,B1,F,0) have a maximum number of uncontrollable 
%   eigenvalues and the corresponding transfer-function matrices satisfy 
%       SYSX(lambda) = SYS1(lambda) + SYS2(lambda)*SYSY(lambda). 
%   A minimum dynamic cover of Type 1 is determined to obtain controllable
%   realizations of SYSX and SYSY of the form 
%   SYSX = (Ar-lambda*Er,Br,Cr1,D1) and SYSY = (Ar-lambda*Er,Br,Cr2,0),
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
%   Note: Large values of INFO.tcond or INFO.fnorm indicate possible
%         loss of numerical accuracy. 
%
% [SYSX,INFO,SYSY] = GRMCOVER1(SYS,M1,TOL) uses the compound realization 
%   of SYS = [SYS1 SYS2], where SYS1 has M1 inputs.
%
% Note: GRMCOVER1 also works if the original descriptor system realization
%       has E singular, but SYS2 is proper. In this case, the order 
%       reduction is performed working only with the proper part of SYS1. 
%       The improper part of SYS1 is included without modification in the  
%       resulting realization of SYSX. 
%
% See also GRMCOVER2. 


%  Author:      A. Varga, 30.11.2015.
%  Revision(s): A. Varga, 31.01.2016, 30.10.2016, 09.05.2017, 27.07.2019.
%
%  Method: The state feedback gain F is determined using the method  of [1]
%  to compute Type 1 minimum dynamic covers for standard systems (E = I) 
%  and the method of [2] for proper descriptor systems.   
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

nsys1 = size(A,1); 
isys1 = 1:m1; isys2 = m1+1:m;
sstype = isequal(E,eye(nsys1));

% handle simple cases
if nsys1 == 0
    sysx = sys1(:,isys1);
    info.stdim = 0;
    info.tcond = 1;
    info.fnorm = 0;
    if nargout == 3
        sysy = ss(zeros(m2,m1)); 
        sysy.Ts = Ts;
    end
    return
end

if  m1 == 0
    sysx = ss(zeros(p,m1)); 
    info.stdim = 0;
    info.tcond = 1;
    info.fnorm = 0;
    if nargout == 3
        sysy = ss(zeros(m2,m1)); 
        sysy.Ts = Ts;
    end
    return
end

% reduce to the special controllability form of [1] or [2]
[At,Et,Bt,Ct,ncont,tau] = sl_gstra(10,A,E,B,C,m1,tol); 

if ncont == 0
    sysx = ss(D(:,isys1));
    sysx.Ts = Ts;
    info.stdim = 0;
    info.tcond = 1;
    info.fnorm = 0;
    if nargout == 3
        sysy = ss(zeros(m2,m1)); 
        sysy.Ts = Ts;
    end
    return
end

if  m2 == 0
   if sstype
      sysx = ss(At,Bt,Ct,D,Ts);
   else
      sysx = dss(At,Bt,Ct,D,Et,Ts);
   end
   info.stdim = tau(tau > 0);
   info.tcond = 1;
   info.fnorm = 0;
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
   [sysxf,info,sysy] = grmcover1(sysf,m1,tol); 
   sysx = sysxf + sysi(:,1:m1); 
   return
end


% perform permutation for Type I cover
cind = length(tau)/2; 
ind = zeros(ncont,1);
iodd = 1:2:length(tau);
ind1 = tau(iodd); n1 = sum(ind1);
ind2 = tau(iodd+1); n2 = sum(ind2);
i2 = 1; i1 = 1; ioff1 = 0; ioff2 = 0; 
for i = 1:cind
   ioff1 =  ioff1 + ind1(i);
   ind((i2:i2+ind2(i)-1)+n1) =  (i2:i2+ind2(i)-1)+ioff1; 
   i2 = i2+ind2(i);
   ind(i1:i1+ind1(i)-1) =  (i1:i1+ind1(i)-1)+ioff2; 
   ioff2 =  ioff2 + ind2(i);
   i1 = i1+ind1(i);
end
At = At(ind,ind); Bt = Bt(ind,:); Ct = Ct(:,ind); 

if ~sstype
   Et = Et(ind,ind); 
   % check invertibility of the relevant part of Et 
   i1 = ind1(1)+1:n1;
   if ~isempty(i1) && (rcond(Et(i1,i1)) < eps || ...
      norm(Et(i1,i1),1) < numel(Et)*eps(norm(Et,1)))
      error('The system [SYS1 SYS2] is possibly non-proper')
   end
end

% anihilate lower left blocks of At and Et
tcond = ncont; % the Frobenius condition number of identity matrix
nlowc1 = n1;             nlowc2 = n1+n2-ind2(cind);
nlowr1 = n1+n2;          nlowr2 = nlowr1;
for i=cind:-1:2
   nlowc1 = nlowc1-ind1(i); nlowc2 = nlowc2-ind2(i-1); 
   nlowr1 = nlowr1-ind2(i); nlowr2 = nlowr2-ind2(i);
   ic1 = nlowc1+1:n1; ic2 = nlowc2+1:nlowc2+ind2(i-1);
   ir2 =  nlowr2+1:nlowr2+ind2(i); 
   x = -At(ir2,ic2)\At(ir2,ic1); 
   tcond = max(tcond,ncont+norm(x,'fro')^2);
   ir = 1:nlowr1;
   At(ir,ic1) = At(ir,ic2)*x+At(ir,ic1); 
   At(ir2,ic1) = zeros(ind2(i),n1-nlowc1);
   Ct(:,ic1) = Ct(:,ic2)*x+Ct(:,ic1); 
   if sstype
      x = -x; At(ic2,:) = At(ic2,:)+x*At(ic1,:);
      Bt(ic2,:) = Bt(ic2,:)+x*Bt(ic1,:);
   else
      Et(ir,ic1) = Et(ir,ic2)*x+Et(ir,ic1); 
      y = -Et(ic2,ic1)/Et(ic1,ic1); 
      tcond = max(tcond,ncont+norm(y,'fro')^2);
      Et(ic2,:) = Et(ic2,:)+y*Et(ic1,:);
      At(ic2,:) = At(ic2,:)+y*At(ic1,:);
      Bt(ic2,:) = Bt(ic2,:)+y*Bt(ic1,:);
   end
end
   

% form the reduced system
ic = 1:n1; i3 = n1+1:n1+ind2(1); 
f2 = -Bt(i3,isys2)\At(i3,ic);
if sstype
   sysx = ss(At(ic,ic)+Bt(ic,isys2)*f2,Bt(ic,isys1),Ct(:,ic)+D(:,isys2)*f2,D(:,isys1),Ts);
else
   sysx = dss(At(ic,ic)+Bt(ic,isys2)*f2,Bt(ic,isys1),Ct(:,ic)+D(:,isys2)*f2,D(:,isys1),Et(ic,ic),Ts);
end
info.fnorm = norm(f2,'fro');
info.tcond = tcond; 
info.stdim = ind1; 
if nargout == 3
   if sstype
      sysy = ss(At(ic,ic)+Bt(ic,isys2)*f2,Bt(ic,isys1),f2,zeros(m2,m1),Ts);
   else
      sysy = dss(At(ic,ic)+Bt(ic,isys2)*f2,Bt(ic,isys1),f2,zeros(m2,m1),Et(ic,ic),Ts);
   end
end

% end GRMCOVER1
end
    


