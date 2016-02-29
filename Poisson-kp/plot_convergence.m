function [SB1,SB2]= plot_convergence(mesh,kx,ky,EFn,EFh,mate_sys)
% check the Fermi-levels, kp subbands and the targets used by kp model
kgrid  = sqrt(kx.^2 + ky.^2);
SB1 = zeros(mesh.ncb,mesh.nk);
SB2 = zeros(mesh.nvb,mesh.nk);
XV1 = zeros(8*mesh.nn,mesh.ncb,mesh.nk);
XV2 = zeros(8*mesh.nn,mesh.nvb,mesh.nk);
Q=  1.6021917e-19;
kB= 1.3807*1e-23;  
T = 300;
%---------------------------------------------------------------------SOLVE
for ik = 1:mesh.nk;
% 
fprintf('k, 1/A %e %e\n',[kx(ik) ky(ik)])
switch (mate_sys)
case{'AlGaN'}
[lmb1,lmb2,xv1,xv2] = solve_kp88_AlGaN(mesh,kx(ik),ky(ik));
case{'InGaN'}
[lmb1,lmb2,xv1,xv2] = solve_kp88_InGaN(mesh,kx(ik),ky(ik));
otherwise, error('specify correct material system!');
end
%
for ib = 1:mesh.ncb;
xv = xv1(:,ib);
f = reshape(xv,8,mesh.nn);
g = sum(abs(f).^2,1); 
xv = xv / sqrt(trapz(mesh.x,g));
XV1(:,ib,ik) = xv; 
SB1(ib,ik) = lmb1(ib); end 
%
for ib = 1:mesh.nvb;
xv = xv2(:,ib);
f = reshape(xv,8,mesh.nn);
g = sum(abs(f).^2,1); 
xv = xv / sqrt(trapz(mesh.x,g));
XV2(:,ib,ik) = xv; 
SB2(ib,ik) = lmb2(ib); end;
end
%---------------------------------------------------------------------SOLVE
%
%-----------------------------------------------------------------------OVP
for ik = 1:(mesh.nk-1);
    xv0 = XV1(:,:,ik);
    xv =  XV1(:,:,ik+1);
    ovp = abs(xv0'*xv);
    [ovpmax,ii] = max(ovp);
    XV1(:,ii,ik+1) = XV1(:,:,ik+1); 
    SB1(ii,ik) = SB1(:,ik); end
%
for ik = 1:(mesh.nk-1);
    xv0 = XV2(:,:,ik);
    xv =  XV2(:,:,ik+1);
    ovp = abs(xv0'*xv);
    [ovpmax,ii] = max(ovp);
    XV2(:,ii,ik+1) = XV2(:,:,ik+1);
    SB2(ii,ik) = SB2(:,ik); end
%-----------------------------------------------------------------------OVP
%
fcb1 = 1.0D0./(1.0D0 + exp( (SB1(1,:)-EFn)/(kB*T/Q)));                     %distrib of first C-subband
fvb1=  1.0D0./(1.0D0 + exp(-(SB2(1,:)-EFh)/(kB*T/Q)));                     %distrib of first V-subband
%fcb3 = 1.0D0./(1.0D0 + exp( (SB1(3,:)-EFn)/(kB*T/Q)));                     %distrib of second C-subband
%fvb3=  1.0D0./(1.0D0 + exp(-(SB2(3,:)-EFh)/(kB*T/Q)));                     %distrib of second V-subband
%
figure(2)
title('current kp results and Fermi level')
subplot(1,2,1),hold on
for ii= 1: mesh.ncb
plot(kgrid,SB1(ii,:), 'k.');end
plot(kgrid,mesh.target1,'*r')                                              %eigs target
plot(kgrid,EFn,'sb')                                                       %fermi level @ convergence
plot(fcb1/10,SB1(1,:),'db')                                                %distribution@ congergence
%plot(fcb3/10,SB1(3,:),'db')
set(gca,'FontSize',14,'FontName','Arial','Box','on')
Emin= mesh.target1-0.2;Emax=EFn+1;
axis([0 max(kgrid) Emin Emax])
xlabel('kx')
ylabel('Energy, eV')
%
subplot(1,2,2),hold on
for ii= 1: mesh.nvb
plot(kgrid,SB2(ii,:), 'r.');end
plot(kgrid,mesh.target2,'*k')
plot(kgrid,EFh,'sg')
plot(fvb1/10,SB2(1,:),'dg')
%plot(fvb3/10,SB2(3,:),'dg')
set(gca,'FontSize',14,'FontName','Arial','Box','on')
Emin=EFh-0.5;Emax=mesh.target2+0.1;
axis([0 max(kgrid) Emin Emax])
xlabel('kx')
ylabel('Energy, eV')
end