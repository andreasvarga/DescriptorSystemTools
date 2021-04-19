function [hanorm,hs] = ghanorm(sys)
%GHANORM  Hankel-norm of a proper and stable descriptor system.
%        HANORM = GHANORM(SYS)  returns the Hankel-norm of the
%        proper and stable descriptor system SYS.
%
%        [HANORM,HS] = GHANORM(SYS) also returns in HS the Hankel singular
%        values of SYS.

%  Author:      A. Varga, 17.01.2016.
%  Revision(s): A. Varga, 22.06.2017.

narginchk(1,1)
nargoutchk(0,2)

if ~isa(sys,'ss')
   error('The input system SYS must be an SS object')
end

% eliminate non-dynamic modes if possible
s2eps = sqrt(eps);
if rcond(sys.e) < s2eps
   sys = gss2ss(sys,s2eps,'triu');
   if rcond(sys.e) < s2eps
       error('The system SYS is not proper')
   end
end

[a,b,c,~,e,Ts] = dssdata(sys);
standsys = isequal(e,eye(size(a,1)));  
discr = (Ts ~= 0);

if size(a,1) == 0,
    % for a non-dynamic system, we set the Hankel norm to zero,
    % but the Hankel singular values are empty
    hanorm = 0; 
    hs = [];
    return
end

if  standsys
    % reduce the system to Schur coordinate form
    [Q,as] = schur(a,'real');
    ev = ordeig(as); 
    % check stability
    if (discr && max(abs(ev)) >= 1-sqrt(eps)) || ...
       (~discr && max(real(ev)) >= -sqrt(eps))
          error('The system SYS is unstable')
    end
    bs = Q'*b; cs = c*Q; es = e;
 else
    % reduce the system to generalized Schur coordinate form
    [as,es,Q,Z] = qz(a,e,'real');
    ev = ordeig(as,es); 
    if (discr && max(abs(ev)) >= 1-sqrt(eps)) || ...
       (~discr && max(real(ev)) >= -sqrt(eps) )
          error('The system SYS is unstable')
    end
    bs = Q*b; cs = c*Z; 
end

% compute the Cholesky factors of Gramians: P = S*S' and Q = R'*R.  
if discr
   flag = [ 1 1]; trans = 1;
   S = sl_glme(4,as,es,bs,flag,trans);
   trans = 0;
   R = sl_glme(4,as,es,cs,flag,trans);
else
   flag = [ 0 1]; trans = 1;
   S = sl_glme(4,as,es,bs,flag,trans);
   trans = 0;
   R = sl_glme(4,as,es,cs,flag,trans);
end

% compute the Hankel singular values
if standsys
   hs = svd(R*S); 
else
   hs = svd(R*es*S); 
end

hanorm = hs(1);  % set Hankel norm 

% end GHANORM
end
