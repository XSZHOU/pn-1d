%==========================================================================
function [out]=f_charge_e(mesh,kgrid,SB1,EFn,eDensity_target)
%==========================================================================
%
kB = 1.38D-23;      % Boltzmann constant, J/K
Q = 1.6021917D-19;  % electron charge
T = 300;            % K
kgrid = kgrid*1e10; % 1/m 
%==========================================================================
eDensity = 0;
for ib=1:mesh.ncb;
engy = SB1(ib,:); % eV
f = 1.0D0./(1.0D0 + exp((engy-EFn)/(kB*T/Q)));
%
eDensity = eDensity + 1/(2*pi)^2*(2*pi)*trapz(kgrid,f.*kgrid); end
out = abs((eDensity_target - (eDensity*1.0D-4))/eDensity_target);
%==========================================================================
