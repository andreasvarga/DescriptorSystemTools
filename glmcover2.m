function [sysx,info,sysy] = glmcover2(sys1,sys2,tol)
% GLMCOVER2  Left minimum dynamic cover of Type 2 based order reduction.
% [SYSX,INFO,SYSY] = GLMCOVER2(SYS1,SYS2,TOL) determines for the proper  
%   system SYS = [ SYS1; SYS2 ] with the descriptor system realization 
%   (A-lambda*E,B,[C1;C2],[D1;D2]) with E invertible, an output injection F 
%   and a feedforward gain G such that the descriptor system realizations 
%   SYSX = (A+F*C2-lambda*E,B+F*D2,C1+G*C2,D1+G*D2) and 
%   SYSY = (A+F*C2-lambda*E,F,C1+G*C2,G) have a maximum number of 
%   unobservable eigenvalues and the corresponding transfer-function 
%   matrices satisfy 
%       SYSX(lambda) = SYS1(lambda) + SYSY(lambda)*SYS2(lambda). 
%   A minimum dynamic cover of Type 2 is determined to obtain observable
%   realizations of SYSX and SYSY of the form 
%   SYSX = (Ar-lambda*Er,Br1,Cr,Dr1) and SYSY = (Ar-lambda*Er,Br2,Cr,Dr2),
%   with the pair (Ar-lambda*Er,Cr) in an observability staircase form. 
%   TOL is a tolerance for rank determinations 
%   (default: TOL = 0, i.e., internally computed).
%
%   The resulting INFO structure contains additional information:
%   INFO.stdim is a vector which contains the dimensions of the diagonal 
%        blocks of Ar-lambda*Er, which are the column dimensions of the  
%        full column rank subdiagonal blocks of the pencil 
%        [Ar-lambda*Er;Cr] in observability staircase form.
%   INFO.tcond is the maximum of the Frobenius-norm condition numbers of 
%        the employed non-orthogonal transformation matrices.
%   INFO.fnorm is the Frobenius-norm of the employed output-injection 
%        matrix F. 
%   INFO.gnorm is the Frobenius-norm of the employed feedforward gain G. 
%   Note: Large values of INFO.tcond, INFO.fnorm or INFO.gnorm indicate 
%         possible loss of numerical accuracy. 
%
% [SYSX,INFO,SYSY] = GLMCOVER2(SYS,M1,TOL) uses the compound realization 
%     of SYS = [SYS1; SYS2], where SYS1 has M1 outputs.
% 
% Note: GLMCOVER2 also works if the original descriptor system realization
%       has E singular, but SYS2 is proper. In this case, the order 
%       reduction is performed working only with the proper part of SYS1. 
%       The improper part of SYS1 is included in the resulting realization 
%       of SYSX. 
%
% See also GLMCOVER1. 

%  Author:      A. Varga, 30.11.2015.
%  Revision(s): A. Varga, 31.01.2016, 11.05.2017, 27.07.2019.
%
%  Method: The output-injection gain F and a feedforward gain G are
%  determined using the dual of the method  of [1] to compute Type 1 
%  minimum dynamic covers for standard systems (E = I) and the dual of 
%  the method of [2] for proper descriptor systems.   
%
%  Note: The resulting McMillan degree of SYSX is the least achievable one
%  provided the realization (A-lambda*E,B,C2,D2) is maximally controllable
%  (i.e., the pair (A+F*C2-lambda*E,B+F*D2) is controllable for any F). 
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
end

if ~isa(sys1,'ss')
   error('The input system SYS1 must be a state-space system object')
end

if ~isa(sys2,'ss')
   if ~isa(sys2,'double') 
      error('SYS2 must be an SS object or a nonnegative integer')
   else
      validateattributes(sys2,{'double'},{'integer','scalar','>=',0,'<=',size(sys1,1)},'','P1')
   end
else
   if size(sys1,2) ~= size(sys2,2) 
      error('The systems SYS1 and SYS2 must have the same number of inputs')
   end    
end  

if nargout <= 2
   [sysx,info] = grmcover2(sys1.',sys2.',tol);
elseif nargout == 3
    [sysx,info,sysy] = grmcover2(sys1.',sys2.',tol);
end
sysx = xperm(sysx,order(sysx):-1:1).';
if nargout == 3
   sysy = xperm(sysy,order(sysy):-1:1).';
end
info.stdim = flip(info.stdim); 

% end GLMCOVER2
end
    


