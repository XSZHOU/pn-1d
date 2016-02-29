%==========================================================================
function [eDensity,hDensity]=f_charge_evp(mesh,kgrid,SB1,SB2,XV1,XV2,EFn,EFh)
%==========================================================================
%
kB = 1.38D-23;      % Boltzmann constant, J/K
Q = 1.6021917D-19;  % electron charge
T = 300;            % K
kgrid = kgrid*1e10; % 1/m
%==========================================================================
% M0 = 9.1095D-31;    % electron mass, Kg
% HBAR = 1.05457D-34; % J*s
% EFn = 0; CBmin = 0;
% eMEFF = 0.2D0 * M0;
% eDensity_ana = eMEFF*(kB*T)/(pi*(HBAR^2));
% eDensity_ana = eDensity_ana * log(1.0D0 + exp((EFn-CBmin)/(kB*T/Q)));
% fprintf('EDENSITY_ANA, (1/cm^2): %e\n',eDensity_ana * 1.0D-4) 
% %
% num_kvectors = 401;               
% engy = CBmin + (HBAR^2)*(kgrid.^2)/(2.0D0*eMEFF)/Q;
% f = 1.0D0./(1.0D0 + exp((engy-EFn)/(kB*T/Q)));
% eDensity = 2/(2*pi)^2*(2*pi)*trapz(kgrid,f.*kgrid);
% fprintf('EDENSITY,     (1/cm^2): %e\n',eDensity * 1.0D-4) 
%==========================================================================
eDensity = 0;
for ib=1:mesh.ncb;
engy = SB1(ib,:); % eV
f = 1.0D0./(1.0D0 + exp((engy-EFn)/(kB*T/Q)));
%
fz2=zeros(mesh.nn,length(kgrid));
for ik2 = 1:length(kgrid);
fz=reshape(XV1(:,ib,ik2),8,mesh.nn);
fz2(:,ik2)=sum(conj(fz).*fz,1); end
%
f= bsxfun(@times,fz2,f.*kgrid);%1996Chuang_JQE,eq.(37)-eq.(38)
eDensity = eDensity + 1/(2*pi)^2*(2*pi)*trapz(kgrid,f,2); end
%eDensity = eDensity + 1/(2*pi)^2*(2*pi)*trapz(kgrid,f.*kgrid); end
%fprintf('EDENSITY,     (1/cm^2): %e\n',eDensity * 1.0D-4) 
%==========================================================================
hDensity = 0;
for ib=1:mesh.nvb;
engy = SB2(ib,:); % eV
f = 1.0D0./(1.0D0 + exp(-(engy-EFh)/(kB*T/Q)));
%
fz2=zeros(mesh.nn,length(kgrid));
for ik2 = 1:length(kgrid);
fz=reshape(XV2(:,ib,ik2),8,mesh.nn);
fz2(:,ik2)=sum(conj(fz).*fz,1); end
%
f= bsxfun(@times,fz2,f.*kgrid);
hDensity = hDensity + 1/(2*pi)^2*(2*pi)*trapz(kgrid,f,2); end
%hDensity = hDensity + 1/(2*pi)^2*(2*pi)*trapz(kgrid,f.*kgrid); end
%fprintf('HDENSITY,     (1/cm^2): %e\n',hDensity * 1.0D-4) 
%==========================================================================