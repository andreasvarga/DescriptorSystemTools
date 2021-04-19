function [ev,poleupd] = eigselect1(poles,sdeg,polref)
% EIGSELECT1  Selection of a real eigenvalue to be assigned.
%  [EV,POLEUPD] = EIGSELECT1(POLES,SDEG,POLREF) selects, from a set of 
%  eigenvalues POLES, a real eigenvalue EV which is the nearest one to the 
%  reference value POLREF. The resulting POLEUPD contains the remaining 
%  eigenvalues from POLES. If POLES is empty, then the selected eigenvalue 
%  EV is set to the stability degree value SDEG. If both POLES and SDEG are 
%  empty, then EV is set to POLREF. If POLES does not contain real
%  eigenvalues, EV is set to an empty matrix. 
%
%  See also EIGSELECT2.

%  Author:    A. Varga 08.12.2016.
%  Revisions: 

if isempty(poles)
   if isempty(sdeg)
       ev = polref;
   else
       ev = sdeg;
   end
   poleupd = [];
else
   polr = poles(imag(poles)==0); poli = poles(imag(poles)~=0);
   if isempty(polr)
       ev = []; poleupd = poles;
   else
       [~,i] = min(abs(polr-polref));
       ev = polr(i);
       poleupd = [polr(1:i-1,1); polr(i+1:end,1); poli];
   end
end

% end EIGSELECT1
end
