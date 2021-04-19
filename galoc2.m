function [f,u,v] = galoc2(a,e,b,poles,tola,tolb)
%GALOC2 Generalized pole assignment for second order systems.
%  [F,U,V] = GALOC2(A,E,B,P,TOLA,TOLB) calculates for the 
%  second order descriptor pair (A-lambda*E,B) with E invertible, the
%  state-feedback gain F such that the generalized eigenvalues of the pair
%  (A+B*F,E) are equal to the complex conjugate pair contained in POLES.
%  TOLA and TOLB are thresholds for nonzero elements in A and B,
%  respectively, and are used for controllability checks.
%  If the pair (A-lambda*E,B) is uncontrollable, then F = [] and U and V
%  contains orthogonal transformation matrices such that the transformed
%  pair (U'*A*V-lambda*U'*E*V,U'*B) is in the form
%
%                                  [ A11-lambda*E11       X        B1 ]
%   [ U'*A*V-lambda*U'*E*V U'*B] = [                                  ] ,
%                                  [     0          A22-lambda*E22 0  ]
% 
%  where the pair (A11-lambda*E11,B1) is controllable. If norm(B) < TOLB,
%  then U and V are the 2x2 identity matrices.


%  Author:    A. Varga 22.11.2015.
%  Revisions: 


% check controllability and determine rank of B
m = size(b,2);
[u,s,v1] = svd(b); s1 = s(1,1);
if m == 1
   s2 = 0;
else
   s2 = s(2,2);
end
rankB = (s1 > tolb) + (s2 > tolb);
if rankB == 0
   % return if norm(B) < TOLB
   f = []; u = eye(2); v = u;
   return
end
if isempty(e) || isequal(e,eye(2))
   at = u'*a*u;  
   if  s2 <= tolb && abs(at(2,1)) <= tola
      % return if rank(B) == 1 and (A-lambda*I,B) uncontrollable
      f = [ ]; v = u; return
   end
else
   at = u'*a; et = u'*e;
   % determine v such that et = u'*e*v is upper triangular
   [v,~] = qr([et(2,2); et(2,1)]); v = v([2,1],[2,1])';  
   at = at*v; 
   if s2 <= tolb && abs(at(2,1)) <= tola
      % return if rank(B) == 1 and (A-lambda*E,B) uncontrollable
      f = [ ]; return
   end
   % reduce to standard case
   [u,s,v1] = svd(e\b); s1 = s(1,1);
   at = u'*(e\a)*u;
end
if rankB == 2
   % try with a direct solution, without inversion
   gamma = [real(poles(1)) imag(poles(1)); imag(poles(2)) real(poles(2))]; 
   ftry = -b\(a-e*gamma);
end
sc = at(1,1)+at(2,2);
sp = poles(1)+poles(2); pp = poles(1)*poles(2);
k11 = (sc-sp)/s1; 
k12 = at(2,2)/at(2,1)*k11+(pp-at(1,1)*at(2,2)+at(1,2)*at(2,1))/at(2,1)/s1;
f = v1*[-k11 -k12; zeros(m-1,2)]*u'; v = eye(m); 
if rankB == 2 && norm(f) > norm(ftry)
   % choose the lower norm feedback 
   f = ftry;
end

% end GALOC2
end
