% DSTOOLSDEMO    Demonstration of the Descriptor System Tools (DSTOOLS).
%                A. Varga.
clc
format compact
echo on
   
% DSTOOLS contains a collection of functions for the operation on and 
% manipulation of rational transfer-function matrices via their 
% descriptor system realizations. The DSTOOLS collection relies on the 
% Control System Toolbox and several mex-functions based on the 
% Systems and Control Library SLICOT. 

pause  % Press any key to continue after pauses.

% The current release of DSTOOLS is version 1.0.0.3, dated April 24, 2021.

% The available funtions in DSTOOLS are:
%
% System analysis.
%   gpole     - Poles of a LTI descriptor system.
%   gzero     - Zeros of a LTI descriptor system. 
%   gnrank    - Normal rank of the transfer function matrix of a LTI system.
%   ghanorm   - Hankel norm of a proper and stable LTI descriptor system.
%   gnugap    - Nu-gap distance between two LTI systems.
%
% System order reduction.
%   gir       - Irreducible realizations of LTI descriptor systems.
%   gminreal  - Minimal realization of a LTI descriptor system.
%   gbalmr    - Balancing-based model reduction of a LTI descriptor system.
%   gss2ss    - Conversions to SVD-like forms without non-dynamic modes.
%
% Operations on transfer function matrices.
%   grnull    - Right nullspace basis of a transfer function matrix.
%   glnull    - Left nullspace basis of a transfer function matrix.
%   grange    - Range space basis of a transfer function matrix. 
%   gcrange   - Coimage space basis of a transfer function matrix. 
%   grsol     - Solution of the linear rational matrix equation G*X = F.
%   glsol     - Solution of the linear rational matrix equation X*G = F.
%   ginv      - Generalized inverses.
%   gsdec     - Generalized additive spectral decompositions.
%   grmcover1 - Right minimum dynamic cover of Type 1 based order reduction.
%   glmcover1 - Left minimum dynamic cover of Type 1 based order reduction.
%   grmcover2 - Right minimum dynamic cover of Type 2 based order reduction.
%   glmcover2 - Left minimum dynamic cover of Type 2 based order reduction.
%   gbilin    - Generalized bilinear transformation.
%
% Factorizations of transfer function matrices.
%   grcf      - Right coprime factorization with proper and stable factors.
%   glcf      - Left coprime factorization with proper and stable factors.
%   grcfid    - Right coprime factorization with inner denominator.
%   glcfid    - Left coprime factorization with inner denominator.
%   gnrcf     - Normalized right coprime factorization.
%   gnlcf     - Normalized left coprime factorization.
%   giofac    - Inner-outer/QR-like factorization.
%   goifac    - Co-outer-co-inner/RQ-like factorization.
%   grsfg     - Right spectral factorization of gamma^2*I-G'*G.
%   glsfg     - Left spectral factorization of gamma^2*I-G*G'.
%
% Model-matching problem.
%   grasol    - Approximate solution of the linear rational matrix equation G*X = F.
%   glasol    - Approximate solution of the linear rational matrix equation X*G = F.
%   glinfldp  - Solution of the least distance problem min||G1-X G2||_inf.
%   gnehari   - Generalized Nehari approximation.
%
% Feedback stabilization. 
%   gsfstab    - Generalized state-feedback stabilization.
%
% Pencil similarity transformations
%   gklf      - Kronecker-like staircase forms of a linear matrix pencil. 
%   gsklf     - Special Kronecker-like form of a system matrix pencil. 
%   gsorsf    - Specially ordered generalized real Schur form.

%%
pause % Press any key to continue ...
clc
% DSTOOLS can manipulate improper rational matrices 

s = tf('s'); z = tf('z');     % define the complex variables s and z                     

pause % Press any key to continue ...

Gc = [s^2 s/(s+1); 0 1/s]     % define the 2-by-2 improper Gc(s)

pause % Press any key to continue ...

Gd = [z^2 z/(z-2); 0 1/z]     % define the 2-by-2 improper Gd(z)

pause % Press any key to continue ...
clc
% build LTI descriptor realizations of Gc(s) and Gd(z) using
% the Control System Toolbox
sysc = ss(Gc)      % build continuous-time descriptor system realization
sysd = ss(Gd);     % build discrete-time descriptor system realization

pause % Press any key to continue ...

% both realizations are actually minimal; unfortunately, checking minimality
% is not possible within the Control System Toolbox since MINREAL 
% fails for improper descriptor systems: 
try
   minreal(sysc);
catch ME
   ME.message
end

pause % Press any key to continue ...
clc
% we can use the functions GIR or GMINREAL from DSTOOLS instead

syscir = gir(sysc)       % computes an irreducible realization

syscmr = gminreal(sysc); % computes a minimal realization
order(syscmr)

% since no order reduction takes place, syscir = sysc and syscmr = sysc !

pause % Press any key to continue ...

%%
clc
% other functions of the Control System Toolbox have also their 
% limitations, when performing on improper descriptor systems

pole(sysc)   %  computes only the finite poles

tzero(sysc)  %  computes only the finite zeros
pause

% we can use instead the functions GPOLE and GZERO from DSTOOLS 

gpole(sysc)   %  computes all poles (finite and infinite)

gzero(sysc)   %  computes all zeros (finite and infinite)

pause % Press any key to continue ...
clc
% the Control Toolbox allows to build descriptor systems 
% (A-s*E,B,C,D) with singular pole pencil A-s*E
A = 0; E = 0; B = 1; C = 1; D = 0;
syst = dss(A,B,C,D,E);

% this may lead to misleading warnings or questionable results 

pole(syst)       % the pole must be NaN, similar to eig(A,E)

pause % Press any key to continue ...

tf(syst)         % the transfer function is infinite and not NaN!

pause % Press any key to continue ...

evalfr(syst,1)   % this evaluation of frequency response is correct

pause % Press any key to continue ...

isproper(syst)   % this test is wrong, because the system is not proper

pause % Press any key to continue ...

% GPOLE can be used to check regularity of the pole pencil
% and to compute the "correct" value of the pole, which is NaN 

gpole(syst)

pause % Press any key to continue ...

%%
clc
% let's try some decomposition and factorization functions

% to compute the separation of proper and polynomial parts of Gc(s), 
% as Gc(s) = Gcp(s) + Gci(s), we can use

[sysf,sysi] = gsdec(sysc);

% for checking the results, convert the terms to zeros/poles/gain form

zpk(sysf)   % proper part Gcp(s)

zpk(sysi)   % polynomial part Gci(s)

pause % Press any key to continue ...

% checking the decomposition
OK = gnrank(sysc-sysf-sysi,[],rand) == 0  % Gc(s)-Gcp(s)-Gci(s) = 0

pause % Press any key to continue ...

%%
clc
% to compute the stable and proper right coprime factorization of Gc(s), as
% Gc(s) = N(s)*inv(M(s)), with N(s) and M(s) stable and proper, 
% having a stability degree -1 and poles assigned to [-2,-3,-4], we can use

[sysn,sysm] = grcf(sysc,struct('poles',[-2,-3,-4],'sdeg',-1,...
    'mindeg',true,'mininf',true));

pause % Press any key to continue ...

% checking the factorization 

OKfact = gnrank(sysc*sysm-sysn,[],rand) == 0 % Gc(s)*M(s)-N(s) = 0

pause % Press any key to continue ...

% checking properness of factors
OKproper = isproper(sysm) & isproper(sysn)    

pause % Press any key to continue ...

% checking stability of poles
OKstab = max(real(gpole(sysm))) < 0 & max(real(gpole(sysn))) < 0    

pause % Press any key to continue ...

% checking coprimeness:  [N(s);M(s)] has no zeros
OKcoprime = isempty(gzero(gminreal([sysn;sysm])))

pause % Press any key to continue ...

%%
clc
% to compute the right coprime factorization with inner denominator of 
% Gd(z), as Gd(z) = N(z)*inv(M(z)), with N(z) and M(z) stable and proper,
% and M(z) inner, we can use

[sysni,sysmi] = grcfid(sysd);

pause % Press any key to continue ...

% checking the factorization

OKfact = gnrank(sysd*sysmi-sysni,[],rand) == 0      %  Gd(z)*M(z)-N(z) = 0

pause % Press any key to continue ...

% checking the innerness of M(z)

OKinner = gnrank(sysmi'*sysmi-eye(2),[],rand) == 0  % conj(M(z))*M(z)-I = 0

pause % Press any key to continue ...

% checking properness of factors
OKproper = isproper(sysmi) & isproper(sysni)    

pause % Press any key to continue ...

% checking stability of poles
OKstab = max(abs(gpole(sysmi))) < 1 & max(abs(gpole(sysni))) < 1    

pause % Press any key to continue ...

% checking coprimeness:  [N(s);M(s)] has no zeros
OKcoprime = isempty(gzero(gminreal([sysni;sysmi])))

pause % Press any key to continue ...

%%
clc
% to illustrate the computation of inner-outer factorization we
% consider a stable 3-by-3 transfer-function matrix G(s) with normal rank 2
G = [(s-1)/(s+2) s/(s+2) 1/(s+2);
    0 (s-2)/(s+1)^2 (s-2)/(s+1)^2;
    (s-1)/(s+2) (s^2+2*s-2)/(s+1)/(s+2) (2*s-1)/(s+1)/(s+2)]; 

% a minimal realization can be computed as

sys = minreal(ss(G)); 

pause % Press any key to continue ...

% analysis of some properties

gpole(sys)    % the system is stable

pause % Press any key to continue ...

gzero(sys)    % the system has 2 unstable zeros and an infinite zero

pause % Press any key to continue ...

gnrank(sys)    % the normal rank of G(s) is 2

pause % Press any key to continue ...
clc
% the inner-quasi-outer factorization of G(s) as G(s) = Gi(s)*[Go(s);0], 
% with Gi(s) inner and square, and Go(s) quasi-outer (i.e., full row rank, 
% and without zeros in the open right-half plane), can be computed by using

[sysi,syso] = giofac(sys,struct('tol',1.e-7));  % use a tolerance of 1.e-7  

pause % Press any key to continue ...

% checking the factorization: Gi(:,1:nr)(s)*Go(s)-G(s) = 0
nr = size(syso,1);    %  nr = 2 is also the normal rank of G(s) 
OKfact = gnrank(sysi(:,1:nr)*syso-sys,[],rand) == 0  

pause % Press any key to continue ...

% checking the innerness of Gi(s)
OKinner = gnrank(sysi'*sysi-eye(3),[],rand) == 0   %  conj(Gi(s))*Gi(s)-I = 0

pause % Press any key to continue ...

% if G(s) contains a so-called free inner factor, then this factor
% is included in Gi(s), but the realization of Go(s) is not minimal

syso = gir(syso);                  % a free inner factor is present in G(s) 

pause % Press any key to continue ...

% check that there are no zeros in the open right-half plane or at infinity 
zer = gzero(syso); 
OKzeros = max(real(zer(~isinf(zer)))) <= 0  

pause % Press any key to continue ...

%%
clc
% let's illustrate the solution of linear rational equations G(s)*X(s) = F(s)

% consider the Wang and Davison example (IEEE Trans. Autom. Contr.,1973)
% to determine a right inverse of G(s) by solving G(s)*X(s) = I.

G = [ s+1 s+2; s+3 s^2+2*s; s^2+3*s 0 ].'/(s^2+3*s+2);
F = eye(2);

sys = gir(ss(G));  % compute a minimal realization of order 3

pause % Press any key to continue ...

gzero(sys)         % G(s) has no zeros, thus stable right inverses exist

pause % Press any key to continue ...

% we can compute a stable right inverse with poles in [-1 -2 -3], by using
sysx = grsol(sys,ss(F),struct('poles',[-1 -2 -3])); 

pause % Press any key to continue ...

% checking the solution 
OK = gnrank(sys*sysx-F,[],rand) == 0   %  G(s)*X(s) - I = 0

pause % Press any key to continue ...

gpole(sysx)            % check assigned poles

pause % Press any key to continue ...
clc
% we can also compute a right inverse of least order = 2

sysxmin = grsol(sys,ss(F),struct('mindeg',true)); order(sysxmin)

pause % Press any key to continue ...

% checking the solution 
OK = gnrank(sys*sysxmin-F,[],rand) == 0   %  G(s)*X(s) - I = 0

pause % Press any key to continue ...

gpole(sysxmin)            % the least order inverse is unstable

pause % Press any key to continue ...

%%
clc
% we can try to explicitly determine a least order right inverse
% using stable generators (X0(s),XN(s)), where X0(s) is a particular
% stable right inverse satisfying G(s)*X0(s) = I and XN(s) is a stable
% basis of the right nullspace of G(s), satisfying G(s)*XN(s) = 0;
% all solutions of G(s)*X(s) = I are given by X(s) = X0(s)+XN(s)*Y(s),
% where Y(s) is arbitrary

% determine the generators [X0,XN] with poles assigned to [-1 -2 -3]
[~,~,sysgen] = grsol(sys,ss(F),struct('poles',[-1 -2 -3]));

pause % Press any key to continue ...

% compute a least order solution X2(s) = X0(s)+XN(s)*Y2(s), by using 
% order reduction based on minimal dynamic covers of Type 2

[sysx2,~,sysy2] = grmcover2(sysgen,2); order(gss2ss(sysx2))

pause % Press any key to continue ...

% checking the solution 
OK = gnrank(sys*sysx2-F,[],rand) == 0   %  G(s)*X2(s) - I = 0

pause % Press any key to continue ...

gpole(sysx2)            % the least order inverse is unstable

pause % Press any key to continue ...

% checking the minimal cover reduction results: X2(s) = X0(s)+XN(s)*Y2(s)
OK = gnrank(sysx2-sysgen(:,1:2)-sysgen(:,3)*sysy2,[],rand) == 0 

%%
pause % Press any key to continue ...
clc
% solution of a H-infinity model-matching problem  min||X(s)*G(s)-F(s)||_inf
% example taken from Francis' book (1987)
W = (s+1)/(10*s+1);    % weighting function
G = [ -(s-1)/(s^2+s+1); (s^2-2*s)/(s^2+s+1)]*W;
F = [ W; 0 ];

pause % Press any key to continue ...

% employ the gamma-iteration based solution proposed by Francis (1987)

% step 1: compute the extended outer-co-inner factorization to reduce the
%         problem to a least-distance problem

[Xopt,info] = grasol(ss(G),ss(F),struct('mindeg',true,'tol',1.e-7)); info
minreal(zpk(Xopt))

pause % Press any key to continue ...

norm(G*Xopt-F,inf)     % this fully agrees with info.mindist 

pause % Press any key to continue ...
clc
% the solution in Francis' book corresponds to a suboptimal solution

% step 2: compute the solution of the suboptimal least-distance problem
%         no gamma-iteration is performed

% compute the suboptimal solution for gamma = 0.2729
opts = struct('mindeg',true,'tol',1.e-7,'gamma',0.2729);
[Xsub,info] = grasol(ss(G),ss(F),opts); info
minreal(zpk(Xsub))

pause % Press any key to continue ...

norm(G*Xsub-F,inf)     % this fully agrees with suboptdist 

%%
clc
% let's illustrate the use of the bilinear transformation 

s = tf('s');              % define the complex variable s
G  = [s^2 s/(s+1); 0 1/s] % define the 2-by-2 improper G(s)
sys = ss(G);              % build continuous-time descriptor system realization
[p,m] = size(sys);        % get system dimensions

pause % Press any key to continue ...

% pole-zero analysis shows that G(s) has poles and zeros in the origin
% and at infinity
pol = gpole(sys), zer = gzero(sys)

pause % Press any key to continue ...

clc
% define a bilinear transformation with g(s) = (s+0.01)/(1+0.01*s) to make 
% all poles and zeros stable and finite using 
g = (s+0.01)/(1+0.01*s)

pause % Press any key to continue ...

% compute the transformed system 
syst = gbilin(sys,g,struct('tol',1.e-7,'minimal',true));
minreal(zpk(syst),1.e-5)

pause % Press any key to continue ...

% check poles and zeros (all are stable)
pol_new = gpole(syst), zer_new = gzero(syst)

pause % Press any key to continue ...

clc
% compute the nugap distance between models using Vinnicombe's formula
nugap = gnugap(sys,syst)

pause % Press any key to continue ...

% plot the Bode magnitude plot
bodemag(sys,syst,{0.01 100})

pause % Press any key to continue ...

echo off
format
