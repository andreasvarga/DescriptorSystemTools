%sl_gstra  MEX-function for reduction to QR-coordinate form.
%
% [At,Et,Bt,Q] = sl_gstra(2,A,E,B) computes, for a descriptor pair 
% (A-lambda*E,B), an orthogonal transformation matrix Q such that 
% the transformed pair 
%       (At-lambda*Et, Bt) := (Q'*A-lambda*Q'*E, Q'*B) 
% has Et upper triangular. 
%
% See also sl_gstra, sl_gstra3.



