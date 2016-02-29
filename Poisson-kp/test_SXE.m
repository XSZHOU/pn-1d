clear all, % close all
% only one band?
%
load subbands_AlGaN_201.mat % 201 points rgrid, 401 points kgrid
%
EFn = 3.5;
EFh = 0; 
% [eDensity,hDensity] = f_charge(mesh,kgrid,SB1,SB2,EFn,EFh);
OPTIONS = optimset('TolX',1e-9,'MaxIter',1000);
%
hDensity = 1e12; % 1/cm^2
[EFh,FVAL,EXITFLAG] = ...
    fminbnd(@(EFh) f_charge_h(mesh,kgrid,SB2,EFh,hDensity),-0.2,0.2,OPTIONS);
ib2 = 1; [SXE_1e12] = f_SXE(mesh,kgrid,SB1,SB2,XV1,XV2,EFn,EFh,ib2);
%
%
hDensity = 2e12; % 1/cm^2
[EFh,FVAL,EXITFLAG] = ...
    fminbnd(@(EFh) f_charge_h(mesh,kgrid,SB2,EFh,hDensity),-.5,.5,OPTIONS);
ib2 = 1; [SXE_2e12] = f_SXE(mesh,kgrid,SB1,SB2,XV1,XV2,EFn,EFh,ib2);
%
figure(2), hold on
plot(kgrid,abs(SXE_1e12),'r.-','linewidth',2)
plot(kgrid,abs(SXE_2e12),'b.-','linewidth',2)
set(gca,'FontSize',14,'FontName','Times','Box','on')
xlabel('kx, inplane kvectors')
ylabel('Energy, eV')
grid
legend('hDensity = 1e12 1/cm^3','hDensity = 2e12 1/cm^3')
% save test_SXE SXE_1e12 SXE_2e12 kgrid