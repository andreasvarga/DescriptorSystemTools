function [nugapdist,fpeak] = gnugap(sys1,sys2,tol,freq,offset)
%GNUGAP  Computes the nu-gap distance between two LTI systems.
%
%    [NUGAPDIST,FPEAK] = GNUGAP(SYS1,SYS2,TOL) calculates NUGAPDIST, the  
%    nu-gap distance between two LTI systems SYS1 and SYS2 and the 
%    corresponding frequency FPEAK (in rad/TimeUnit), where the nu-gap 
%    distance achieves its peak value. NUGAPDIST satisfies 0 <= NUGAPDIST <= 1. 
%    The value NUGAPDIST = 1 results, if the winding number is different of 
%    zero in which case FPEAK = []. 
%    TOL is a nonnegative scalar to be used for rank determinations 
%    (Default: internally computed, if not specified). 
%
%    [NUGAPDIST,FPEAK] = GNUGAP(SYS1,SYS2,TOL,FREQ) calculates the 
%    pointwise nu-gap distances between two LTI systems SYS1 and SYS2 over 
%    a set of real frequency values FREQ. FPEAK (in rad/TimeUnit) is the 
%    frequency, for which the pointwise distance achieves its peak value.
%
%    [NUGAPDIST,FPEAK] = GNUGAP(SYS1,SYS2,TOL,FREQ,OFFSET) uses OFFSET, the 
%    stability boundary offset, to be used to assess the  
%    finite zeros which belong to the boundary of the stability domain 
%    as follows: in the continuous-time case, these are the finite 
%    zeros having real parts in the interval [-OFFSET, OFFSET], while 
%    in the discrete-time case, these are the finite zeros having moduli 
%    in the interval [1-OFFSET,1+OFFSET]. (Default: OFFSET = sqrt(eps)).

%   Copyright 2018 A. Varga.
%   Author:       A. Varga, 14.09.2018.
%   Revision(s):  A. Varga, 02.10.2018. 
%
%   Method: The evaluation of nu-gap uses the definition proposed in [1],
%   extended to generalized LTI systems. The computation of winding number
%   is based on enhancements covering zeros on the boundary of the 
%   stability domain and infinite zeros.
%
%   References:
%   [1] G. Vinnicombe. 
%       Uncertainty and feedback: H-infinity loop-shaping and the 
%       nu-gap metric. Imperial College Press, London, 2001. 

narginchk(2,5)
nargoutchk(0,3)

if ~isa(sys1,'ss') && ~isa(sys1,'tf') && ~isa(sys1,'zpk')
   error('The input system SYS1 must be an LTI object')
end 

if ~isa(sys2,'ss') && ~isa(sys2,'tf') && ~isa(sys2,'zpk')
   error('The input system SYS2 must be an LTI object')
end 
Ts = sys1.Ts; Ts2 = sys2.Ts;
disc = (Ts ~= 0);
if (disc && (Ts2 == 0)) || (~disc && (Ts2 ~= 0))
   error('The systems SYS1 and SYS2 must have the same type')
end
if disc && (Ts && Ts2 && Ts ~= Ts2) 
   error('The systems SYS1 and SYS2 must have the same sampling period')
end

% ensure SYS1 and SYS2 are in state-space form
if ~isa(sys1,'ss')
   sys1 = ss(sys1);
end

if ~isa(sys2,'ss')
   sys2 = ss(sys2);
end

if nargin == 2
   tol = 0;
else
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','TOL',3)  
end
tolinf = max(tol,sqrt(eps));

if nargin <= 3
   freq = [];
else
   if ~isempty(freq)
      validateattributes(freq, {'double'},{'real','vector','>=',0},'','FREQ',4)
   end
end
nf = length(freq);


if nargin <= 4
   offset = sqrt(eps);
else
   validateattributes(offset, {'double'},{'real','scalar','>',0,'<',1},'','OFFSET',5)  
end

[p,m] = size(sys1);

% % this version should also work
% syst = eye(m)+sys2'*sys1;
% % check invertibility 
% if m == 1
%    singular = (abs(evalfr(syst,rand)) <= tolinf);
% else
%    singular = (rank(evalfr(syst,rand),tolinf) < m);
% end
% if singular
%    nugap = 1; fpeak = [];
%    return
% end
% 
% tpol = gpole(syst,tol);
% tzer = gzero(syst,tol); 
% if disc
%    wno = sum(abs(tzer) > 1) - sum(abs(tpol) > 1);
% else
%    indpi = (abs(tpol) == inf);  
%    indzi = (abs(tzer) == inf);  
%    wno = sum(real(tzer(~indzi)) > 0) + sum(indzi) - sum(real(tpol(~indpi)) > 0) - sum(indpi);
% end
% % check if winding number is zero
% if wno ~= 0
%    % nonzero winding number
%    nugap = 1; fpeak = [];
%    return
% end

% compute the normalized right coprime factorizations
R1 = grange([sys1;eye(m)],struct('inner',true,'tol',tol)); % R1 = [N1;M1]
R2 = grange([sys2;eye(m)],struct('inner',true,'tol',tol)); % R2 = [N2;M2]

% check conditions on det(R2'*R1)
syst = gir(R2'*R1,tol);  % syst = R2'*R1 must also work
[~,infoz] = gzero(syst,tol,offset);
% check invertibility and presence of zeros on the boundary of stability
% domain
if infoz.nrank ~= order(syst)+m || infoz.nfsbz ;
   nugapdist = ones(max(1,nf),1); fpeak = [];
   return
end

% evaluate winding number 
[~,infop] = gpole(syst,tol,offset);
wno = infoz.nfuz - infop.nfuev + infoz.niz - infop.nip;
% check condition on winding number 
if wno ~= 0
   % nonzero winding number
   nugapdist = ones(max(1,nf),1); fpeak = [];
   return
end

% compute the normalized left coprime factorization s.t. L1 = [ N1t M1t]
L1 = gcrange([sys1 eye(p)],struct('coinner',true,'tol',tol)); 
% compute the underlying system to compute the nu-gap distance 
% using the definition of Vinnicombe
syst = L1*[zeros(m,p) -eye(m);eye(p,p+m)]*R2;
if isempty(freq)
   % compute the nu-gap using the definition of Vinnicombe
   [nugapdist,fpeak] = norm(syst,inf,tolinf);
else
   H = freqresp(syst,freq); nugapdist = zeros(nf,1);
   tmax = norm(H(:,:,1)); fpeak = freq(1);
   nugapdist(1) = tmax; 
   for i = 2:nf
       temp = norm(H(:,:,i)); 
       nugapdist(i) = temp; 
       if tmax < temp
          tmax = temp; fpeak = freq(i);
       end
   end          
end

% end GNUGAP

