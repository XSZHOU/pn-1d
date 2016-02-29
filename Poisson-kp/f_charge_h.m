%==========================================================================
function [out]=f_charge_h(mesh,kgrid,SB2,EFh,hDensity_target)
%==========================================================================
%
kB = 1.38D-23;      % Boltzmann constant, J/K
Q = 1.6021917D-19;  % electron charge
T = 300;            % K
kgrid = kgrid*1e10; % 1/m
%==========================================================================
hDensity = 0;
for ib=1:mesh.nvb;
engy = SB2(ib,:); % eV
f = 1.0D0./(1.0D0 + exp(-(engy-EFh)/(kB*T/Q)));
%
hDensity = hDensity + 1/(2*pi)^2*(2*pi)*trapz(kgrid,f.*kgrid); end
% fprintf('HDENSITY,     (1/cm^2): %e\n',hDensity * 1.0D-4) 
out = abs((hDensity_target - (hDensity*1.0D-4))/hDensity_target);
%==========================================================================
