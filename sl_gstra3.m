%sl_gstra  MEX-function for reduction to RQ-coordinate form.
%
% [At,Et,Ct,Z] = sl_gstra(3,A,E,C) computes, for a descriptor pair 
% (A-lambda*E,C), an orthogonal transformation matrix Z such that 
% the transformed pair 
%         (At-lambda*Et, Ct) := (A*Z-lambda*E*Z, C*Z) 
% has Et upper triangular. 
%
% See also sl_gstra, sl_gstra2.



