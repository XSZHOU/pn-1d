%==========================================================================
clear all, close all
clc
%load subbands_AlGaN_401.mat
%==========================================================================
% define FEM grid 
mesh.L  = 10e-9;                           % total length, m
mesh.nn = 201;                             % number of nodes 
mesh.ne = mesh.nn - 1;                     % number of elements
mesh.x  = linspace(0,mesh.L,mesh.nn);      % node coordinates, m
mesh.le = mesh.x(2:end) - mesh.x(1:end-1); % edge length, m
mesh.xc = (mesh.x(1:end-1) + mesh.x(2:end))/2; % center points, m
%
xmol_barrier= 0.25;
xmol_well = 0.0;
bow = 0.98;                                % conduction band bowing
bow_v = 0.33*bow;                          % valence band bowing
V0 = - xmol_barrier*0.8 + bow_v*xmol_barrier*(1-xmol_barrier);% VB offset
%
mesh.xmol = xmol_barrier*ones(1,mesh.ne);
mesh.evb = V0*ones(1,mesh.ne);
%
ii = ((mesh.xc>=4e-9)&(mesh.xc<=6e-9)); 
mesh.evb (ii) = 0; 
mesh.ecb = mesh.evb;
mesh.xmol(ii) = xmol_well;
mesh.target1 =  3.0 ;
mesh.target2 = -3e-3;
% 
num_kvectors = 401;              
kx = linspace(0,0.4,num_kvectors);
ky = zeros(1,num_kvectors);
kgrid  = sqrt(kx.^2 + ky.^2);
mesh.nk = 16; mesh.nvb = 10; mesh.ncb = 6;
%==========================================================================
%
SB1 = zeros(mesh.ncb,num_kvectors);
SB2 = zeros(mesh.nvb,num_kvectors);
XV1 = zeros(8*mesh.nn,mesh.ncb,num_kvectors);
XV2 = zeros(8*mesh.nn,mesh.nvb,num_kvectors);
%==========================================================================
%
for ik = 1:num_kvectors;
% 
fprintf('k, 1/A %e %e\n',[kx(ik) ky(ik)])
[lmb1,lmb2,xv1,xv2] = solve_kp88_AlGaN(mesh,kx(ik),ky(ik));
%
for ib = 1:mesh.ncb;
xv = xv1(:,ib);
f = reshape(xv,8,mesh.nn);
g = sum(abs(f).^2,1); 
xv = xv / sqrt(trapz(mesh.x,g)); %trapz(mesh.x,g)=sum(g)*mesh.L/mesh.ne=1
XV1(:,ib,ik) = xv; 
SB1(ib,ik) = lmb1(ib); end
%
for ib = 1:mesh.nvb;
xv = xv2(:,ib);
f = reshape(xv,8,mesh.nn);
g = sum(abs(f).^2,1); 
xv = xv / sqrt(trapz(mesh.x,g));
XV2(:,ib,ik) = xv; 
SB2(ib,ik) = lmb2(ib); end; end
%==========================================================================
for ik = 1:(num_kvectors-1);
    xv0 = XV1(:,:,ik);
    xv =  XV1(:,:,ik+1);
    ovp = abs(xv0'*xv);
    [ovpmax,ii] = max(ovp);
    XV1(:,ii,ik+1) = XV1(:,:,ik+1); 
    SB1(ii,ik) = SB1(:,ik); end
%
for ik = 1:(num_kvectors-1);
    xv0 = XV2(:,:,ik);
    xv =  XV2(:,:,ik+1);
    ovp = abs(xv0'*xv);
    [ovpmax,ii] = max(ovp);
    XV2(:,ii,ik+1) = XV2(:,:,ik+1);    
    SB2(ii,ik) = SB2(:,ik); end
%
% ovp = NaN*zeros(mesh.ncb,num_kvectors);
% xv0 = XV1(:,:,2);
% for ik = 2:num_kvectors;
%     xv =  XV1(:,:,ik);
%     ovp(:,ik) = diag(abs(xv0'*xv)); 
% end
% %
% ovp = NaN*zeros(mesh.nvb,num_kvectors);
% xv0 = XV2(:,:,2);
% for ik = 2:num_kvectors;
%     xv =  XV2(:,:,ik);
%     ovp(:,ik) = diag(abs(xv0'*xv)); 
% end
%==========================================================================
%
SB1=real(SB1);
SB2=real(SB2);
%
figure(1),   hold on
plot(kgrid,SB1(1,:), 'k.')
plot(kgrid,SB1(2,:), 'k.')
plot(kgrid,SB1(3,:), 'r.')
plot(kgrid,SB1(4,:), 'r.')
plot(kgrid,SB1(5,:), 'g.')
plot(kgrid,SB1(6,:), 'g.')
% plot(kgrid,SB1(7,:), 'b.')
% plot(kgrid,SB1(8,:), 'b.')
% plot(kgrid,SB1(9,:), 'm.')
% plot(kgrid,SB1(10,:),'m.')
% plot(kgrid,SB1(11,:),'k.')
% plot(kgrid,SB1(12,:),'k.')
% plot(kgrid,SB1(13,:),'r.')
% plot(kgrid,SB1(14,:),'r.')
% plot(kgrid,SB1(15,:),'g.')
% plot(kgrid,SB1(16,:),'g.')
% axis([-Kt Kt mesh.Vref mesh.Vref+1.5])
set(gca,'FontSize',14,'FontName','Arial','Box','on')
xlabel('kx, inplane kvectors')
ylabel('Energy, eV')
title('GaN/AlGaN/GaN:C-subbands')
%==========================================================================
figure(2),   hold on
plot(kgrid,SB2(1,:), 'k.')
plot(kgrid,SB2(2,:), 'k.')
plot(kgrid,SB2(3,:), 'r.')
plot(kgrid,SB2(4,:), 'r.')
plot(kgrid,SB2(5,:), 'g.')
plot(kgrid,SB2(6,:), 'g.')
plot(kgrid,SB2(7,:), 'b.')
plot(kgrid,SB2(8,:), 'b.')
plot(kgrid,SB2(9,:), 'm.')
plot(kgrid,SB2(10,:),'m.')
% plot(kgrid,SB2(11,:),'k.')
% plot(kgrid,SB2(12,:),'k.')
% plot(kgrid,SB2(13,:),'r.')
% plot(kgrid,SB2(14,:),'r.')
% plot(kgrid,SB2(15,:),'g.')
% plot(kgrid,SB2(16,:),'g.')
% axis([-Kt Kt Vv 0.01])
set(gca,'FontSize',14,'FontName','Arial','Box','on')
xlabel('kx, inplane kvectors')
ylabel('Energy, eV')
title('Al_{0.25}GaN/GaN/Al_{0.25}GaN:Valence subbands')
%==========================================================================
save subbands_AlGaN0_201 XV1 XV2 SB1 SB2 kgrid mesh
%==========================================================================