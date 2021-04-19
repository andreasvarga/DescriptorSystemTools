function [ev,poleupd] = eigselect2(poles,sdeg,polref,systype)
% EIGSELECT2  Selection of a pair of eigenvalues to be assigned.
%  [EV,POLEUPD] = EIGSELECT2(POLES,SDEG,POLREF,SYSTYPE) selects,
%  from a set of eigenvalues POLES, a complex conjugate pair of eigenvalues 
%  (EV(1),EV(2)), such that EV(1) is the nearest one to the reference value 
%  POLREF. If no complex eigenvalues are available in POLES, a pair of real 
%  eigenvalues is selected, which are the nearest ones to POLREF. 
%  If POLES is empty, the selected eigenvalues are the nearest ones
%  on the boundary of the stability domain Cs, defined by the stability 
%  degreee parameter SDEG and system type flag SYSTYPE. If SYSTYPE = 0, Cs 
%  is the left half complex plane Re(s) <= SDEG, and if SYSTYPE > 0, Cs is 
%  the interior of unit circle |z| <= SDEG. The resulting POLEUPD contains 
%  the remaining eigenvalues. If both POLES and SDEG are empty, then EV is 
%  selected as [ POLREF; conj(POLREF) ]. If POLES contains only a real
%  value and if SDEG is nonempty, then EV = [ POLES; SDEG ], while if
%  SDEG is empty then EV = [ POLES; SDEGDEF ], where SDEGDEF = -0.2 for 
%  SYSTYPE = 0 and SDEGDEF = 0.8 if SYSTYPE > 0.
%
%  See also EIGSELECT1.

%  Author:    A. Varga 08.12.2016.
%  Revisions: 


polref = real(polref)+1i*abs(imag(polref));
if isempty(poles)
   if isempty(sdeg)
      ev = [polref; conj(polref)];
   else
     if systype
        ev = [polref; conj(polref)]; 
        ev = (sdeg/abs(polref))*ev; 
     else
        evi = imag(polref);
        ev = [sdeg+1i*evi; sdeg-1i*evi];
     end
   end
   poleupd = [];
else
   polr = poles(imag(poles)==0); poli = poles(imag(poles)~=0);
   if isempty(poli)
       % select two real eigenvalues
       if length(polr) < 2
          if isempty(sdeg)
             % use default values
             if systype, sdegdef = 0.8; else sdegdef = -0.2; end
             ev = [polr;sdegdef]; 
          else
             ev = [polr;sdeg]; 
          end
          poleupd = [];
       else
          [~,isort] = sort(abs(polr-polref));
          polr = polr(isort); 
          ev = [ polr(1); polr(2)];
          poleupd = [polr(3:end); poli];
       end
   else
       [~,i] = min(abs(poli-polref));
       ev = [poli(i);poli(i+1)];       
       poleupd = [polr; poli(1:i-1,1); poli(i+2:end,1)];
   end
end

% end EIGSELECT2
end

