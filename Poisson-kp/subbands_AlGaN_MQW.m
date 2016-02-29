%==========================================================================
clear all, %close all
clc
%needs very dense k mesh to avoid points jumping
%==========================================================================
% define FEM grid 
mesh.L  = 40e-9;                               % total length, m
mesh.nn = 801;                                 % number of nodes 
mesh.ne = mesh.nn - 1;                         % number of elements
mesh.x  = linspace(-200e-10,200e-10,mesh.nn);  % node coordinates, m
mesh.le = mesh.x(2:end) - mesh.x(1:end-1);     % edge length, m
mesh.xc = (mesh.x(1:end-1) + mesh.x(2:end))/2; % center points, m
%
xmol_barrier= 0.25;
xmol_well = 0.0;
bow = 0.98;                                    % conduction band bowing
bow_v = 0.33*bow;                              % valence band bowing
V0 = - xmol_barrier*0.8 + bow_v*xmol_barrier*(1-xmol_barrier);% VB offset
%
mesh.xmol = xmol_barrier*ones(1,mesh.ne);
mesh.evb = V0*ones(1,mesh.ne);
%
ii = ((mesh.xc>=-7e-9)&(mesh.xc<=-5e-9))| ...
     ((mesh.xc>=-3e-9)&(mesh.xc<=-1e-9))| ...
     ((mesh.xc>= 1e-9)&(mesh.xc<= 3e-9))| ...
     ((mesh.xc>= 5e-9)&(mesh.xc<= 7e-9));
%
i_lower  = find(mesh.xc >=-7e-9,1,'first');
i_higher = find(mesh.xc <= 7e-9,1,'last');
vbias=-0.8;
mesh.bias=zeros(1,mesh.ne);
mesh.bias(i_lower:i_higher)=linspace(0,vbias,i_higher-i_lower+1);
mesh.bias(i_higher:end)=vbias;
%
mesh.evb (ii) = 0; mesh.evb=mesh.evb+mesh.bias;
mesh.ecb = mesh.evb;
mesh.xmol(ii) = xmol_well;
% mesh.target1 =  3.0 ;
mesh.target1 =  2.0 ;
mesh.target2 = -3e-3;
% 
num_kvectors = 201;
kx = linspace(0,0.2,num_kvectors);
ky = zeros(1,num_kvectors);
kgrid  = sqrt(kx.^2 + ky.^2);
mesh.nk = 32; mesh.nvb = 32; mesh.ncb = 32;
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
[lmb1,lmb2,xv1,xv2] = solve_kp88_AlGaN(mesh,kx(ik),ky(ik),1);
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
colarlib={'k.','k.','r.','r.','g.','g.','b.','b.','m.','m.'};
%colarlib={'ko','ko','ro','ro','go','go','bo','bo','mo','mo'};
colarlibs=repmat(colarlib,[1 5]);
figure(1),   hold on
for ii=1:mesh.ncb
     plot(kgrid,SB1(ii,:),colarlibs{ii})
end
set(gca,'FontSize',14,'FontName','Arial','Box','on')
xlabel('kx, inplane kvectors')
ylabel('Energy, eV')
title('AlGaN/GaN:C-subbands')
%==========================================================================
%
colarlib={'k.','k.','r.','r.','g.','g.','b.','b.','m.','m.'};
%colarlib={'ko','ko','ro','ro','go','go','bo','bo','mo','mo'};
colarlibs=repmat(colarlib,[1 5]);
figure(2),   hold on
for ii=1:mesh.nvb
     plot(kgrid,SB2(ii,:),colarlibs{ii})
end
axis([0 max(kx) -0.2 0.0])
set(gca,'FontSize',14,'FontName','Arial','Box','on')
xlabel('kx, inplane kvectors')
ylabel('Energy, eV')
title('AlGaN/GaN:V-subbands')
%
figure(3), hold on 
plot(mesh.xc*1e10,mesh.evb,'k-.')
%
esheet=4e13;hsheet=4e13;
OPTIONS = optimset('TolX',1e-9,'MaxIter',1000);
[EFn,FVAL,EXITFLAG] = ...
fminbnd(@(EFn) f_charge_e(mesh,kgrid,SB1,EFn,esheet), 2.5,3.9,OPTIONS);
[EFh,FVAL,EXITFLAG] = ...
fminbnd(@(EFh) f_charge_h(mesh,kgrid,SB2,EFh,hsheet),-0.5,0.5,OPTIONS);
fprintf('EFn,EFh,%e %e\n',[EFn EFh])
%==========================================================================
%save subbands_AlGaNSL_201 XV1 XV2 SB1 SB2 kgrid mesh
%==========================================================================