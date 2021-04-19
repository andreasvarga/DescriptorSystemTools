%sl_gstra  MEX-function for descriptor system scaling.
%
% [At,Et,Bt,Ct,Dl,Dr] = sl_gstra(1,A,E,B,C,JOB,THRESH) scales, for a 
% descriptor system (A-lambda*E,B,C), the matrices of the system pencil
% 
%           S =  ( A  B ) - lambda * ( E  0 )                (1)
%                ( C  0 )            ( 0  0 )
% 
%  by applying diagonal similarity transformations Dl and Dr to determine
%  (At-lambda*Et,Bt,Ct) = (Dl*A*Dr - lambda*Dl*E*Dr, Dl*B, C*Dr) such that
%  the rows and columns of the matrices of the transformed system pencil
% 
%                diag(Dl,I) * S * diag(Dr,I)
% 
%  are as close in norm as possible. Scaling may reduce the 1-norms
%  of the matrices A, E, B, and C.
% 
%  JOB specifies the option for performing the balancing on the following
%  particular system pencils: 
%  JOB = 0: scale all matrices A, E, B, C using S in (1) (default);
%  JOB = 1: scale only A, E, B using S = ( A-lambda*E  B );
%  JOB = 2: scale only A, E, C using  S = ( A-lambda*E );
%                                         (     C      )
%  JOB = 3: scale using only A, E using S = A-lambda*E.
% 
% THRESH is an optional threshold on the size of the elements of Dl and Dr.
%
% See also sl_gstra.



