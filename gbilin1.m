function [sys1,sysi1] = gbilin1(type,params)
%GBILIN1  Transfer functions of commonly used bilinear transformations. 
%    [SYS1,SYSI1] = GBILIN1(TYPE,PARAMS) generates the first degree 
%    rational transfer functions SYS1 and SYSI1 which define several 
%    commonly used bilinear transformations and their inverses, respectively.
%    The transfer function g(delta) of SYS1 has the form
%                         g(delta) = (a*delta+b)/(c*delta+d), 
%    where for a continuous-time SYS1, delta = s, the complex variable in 
%    the Laplace transform, while for a discrete-time SYS1, delta = z, 
%    the complex variable in the Z-transform. 
%    The transfer function ginv(lambda) of SYSI1 is the transfer function 
%    of the inverse of the bilinear transformation g(delta) and has the form
%                       ginv(lambda) = (d*lambda-b)/(-c*lambda+a) .
%    TYPE is a character string which specifies the type of the bilinear 
%    transformation and PARAMS is an optional vector of parameters.
%    The following options for TYPE can be used, with the corresponding 
%    values of parameters:
%    'c2d'    - continuous-to-discrete transformation
%               s = g(z) = (z-1)/(z+1); z = ginv(s) = (s+1)/(-s+1)
%               PARAMS = [ T ], where:
%               T - sampling time of SYS1 (default: not specified)
%    'd2c'    - discrete-to-continuous transformation
%               z = g(s) = (s+1)/(-s+1);  s = ginv(z) = (z-1)/(z+1);
%               PARAMS = [ T ], where:
%               T - sampling time of SYSI1 (default: not specified)
%    'Tustin' - Tustin Transform (also known as trapezoidal integration)
%               s = g(z) = (2*z-2)/(T*z+T); z = ginv(s) = (T*s+2)/(-T*s+2)
%               PARAMS = [ T ], where:
%               T - sampling time of SYS1 (default: T = 2)
%    'Euler'  - Euler integration (also known as forward Euler integration)
%               s = g(z) = s = (z-1)/T; z = ginv(s) = T*s+1
%               PARAMS = [ T ], where:
%               T - sampling time of SYS1 (default: T = 1)
%    'BEuler' - Backward Euler integration (also known as forward Euler)
%               s = g(z) = (z-1)/(T*z); z = ginv(s) = 1/(-T*s+1)
%               PARAMS = [ T ], where:
%               T - sampling time of SYSI1 (default: T = 1)
%    'gen'    - General bilinear transformation 
%               lambda = g(delta) =  (a*delta+b)/(c*delta+d); 
%               delta = ginv(lambda) = (d*lambda-b)/(-c*lambda+a)
%               PARAMS = [ a b c d T Ti], where:
%               [a b c d ] real parameters (default: [1 0 0 1])
%               T  - sampling time of SYS1  (default: T = 0)
%               Ti - sampling time of SYSI1 (default: T = 0)
%
%    See also  GBILIN.

%  Copyright 2018 A. Varga 
%  Author:    A. Varga, 21.08.2018.
%  Revisions: A. Varga, 06.12.2018.

narginchk(1,2)
nargoutchk(0,2)
if nargin == 1
    params = [];
end

validateattributes(type,{'char'},{'nonempty'},'','TYPE',1) 
if ~isempty(params)
   validateattributes(params, {'double'}, {'real','vector'},'','PARAMS',2)
end


s = tf('s'); z = tf('z');
switch lower(type)
   case 'c2d'
       sys1 = (z-1)/(z+1); sysi1 = (s+1)/(-s+1);
       if ~isempty(params)
           validateattributes(params, {'double'}, {'real','scalar','>',0},'','PARAMS',2)
           sys1.Ts = params(1);
       end
   case 'd2c'
       sys1 = (s+1)/(-s+1); sysi1 = (z-1)/(z+1);
       if ~isempty(params)
           validateattributes(params, {'double'}, {'real','scalar','>',0},'','PARAMS',2)
           sysi1.Ts = params(1);
       end
   case 'tustin'
       if isempty(params)
           T = 2;
       else
           validateattributes(params, {'double'}, {'real','scalar','>',0},'','PARAMS',2)
           T = params(1);
       end
       sys1 = (2*z-2)/(T*z+T); sysi1 = (T*s+2)/(-T*s+2);
       sys1.Ts = T;
   case 'euler'
       if isempty(params)
           T = 1;
       else
           validateattributes(params, {'double'}, {'real','scalar','>',0},'','PARAMS',2)
           T = params(1);
       end
       sys1 = (z-1)/T; sysi1 = T*s+1;
       sys1.Ts = T;
   case 'beuler'
       if isempty(params)
           T = 1;
       else
           validateattributes(params, {'double'}, {'real','scalar','>',0},'','PARAMS',2)
           T = params(1);
       end
       sys1 = (z-1)/(T*z); sysi1 = 1/(-T*s+1);
       sys1.Ts = T;
   case 'gen'
       if isempty(params)
          a = 1; b = 0; c = 0; d = 1; T = 0; Ti = -1;
       else
          if length(params) ~= 6
              error('The number of parameters must be 6')
          end
          a = params(1); b = params(2); c = params(3); d = params(4); 
          T = params(5); Ti = params(6);
       end
       if T == 0
           sys1 = minreal((a*s+b)/(c*s+d)); 
           if order(sys1) == 0
              error('The McMillan degree of SYS1 must be one')
           end
       else
           sys1 = minreal((a*z+b)/(c*z+d)); 
           if order(sys1) == 0
              error('The McMillan degree of SYS1 must be one')
           end
           sys1.Ts = T; 
       end
       if Ti == 0
           sysi1 = (d*s-b)/(-c*s+a);
       else
           sysi1 = (d*z-b)/(-c*z+a);
           sysi1.Ts = Ti;
       end
   otherwise
      error('Unknown type.')
end

% end GBILIN1
