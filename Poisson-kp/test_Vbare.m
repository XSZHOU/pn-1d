%==========================================================================
clear all, %close all
%clc
%==========================================================================
% define FEM grid 
mesh.L  = 10e-9;                      % total length, m
mesh.nn = 401;                        % number of nodes 
mesh.ne = mesh.nn - 1;                % number of elements
mesh.x  = linspace(0,mesh.L,mesh.nn); % node coordinates, m
mesh.le = mesh.x(2:end) - mesh.x(1:end-1); % edge length, m
mesh.xc = (mesh.x(1:end-1) + mesh.x(2:end))/2; % center points, m
%
xmol_barrier = 0.0;              % molar fraction of the barrier 
xmol_well = 0.25;                % molar fraction of the well
bow = 1.4;                       % bowing factor, Moses
mesh.Vref = 3.4527*(1-xmol_well)+0.7577*xmol_well;
Eg_InGaN = 3.437*(1-xmol_well)+0.69*xmol_well; %-bow*xmol_well*(1-xmol_well)
DeltaE = 3.437-Eg_InGaN;
mesh.target1 = 3.4;
mesh.target2 = 3e-3;
%
Vc =  0.7*DeltaE;                % CB offset
Vv = -0.3*DeltaE;                % VB offset
mesh.ecb  = Vc*ones(1,mesh.ne);
mesh.evb  = Vv*ones(1,mesh.ne);     
mesh.xmol = xmol_barrier*ones(1,mesh.ne); 
%
ii = ((mesh.xc>=4e-9)&(mesh.xc<=6e-9)); 
mesh.ecb (ii) = 0;
mesh.evb (ii) = 0;
mesh.xmol(ii) = xmol_well;
mesh.ncb = 16;mesh.nvb = 16;
%==========================================================================
% test rotations
%==========================================================================
kx = 0.04; ky = 0.1; k = [kx ky]; knorm = norm(k);
phi=atan(ky/kx); k_rotated=[knorm 0];
ib = 16; %unconfined 
%
[lmb1,lmb2,xv1,xv2] = solve_kp88_InGaN(mesh,kx,ky);
f1 = reshape(xv2(:,ib),8,mesh.nn);
%
[lmb_1,lmb_2,xv_1,xv_2] = solve_kp88_InGaN(mesh,k_rotated(1),k_rotated(2));
f_1 = reshape(xv_2(:,ib),8,mesh.nn);
%
a=exp(1i*3/2*phi)*sqrt(1/2);
b=exp(1i*1/2*phi)*sqrt(1/2);
V =[b' b  0  0  0  0  0  0;...
    b' -b 0  0  0  0  0  0;...
    0  0  a' 0  0  a  0  0;...
    0  0  0  b  0  0  b' 0;...
    0  0  0  0  b' 0  0  b;...
    0  0  a' 0  0 -a  0  0;...
    0  0  0  b  0  0 -b' 0;...
    0  0  0  0  -b' 0  0 b].';    
%
f1_rotated = V*f_1;
%
figure(1), hold on
plot(mesh.x*1e9,abs(f1),'kd')
plot(mesh.x*1e9,abs(f_1),'b-')
plot(mesh.x*1e9,abs(f1_rotated),'ro')
% 
break
%==========================================================================
% test bare coulomb, non-ordered
%==========================================================================
[xx,yy] = meshgrid(mesh.x,mesh.x);
kx = 0.04; ky = 0.0; k = [kx ky]; knorm = norm(k);
ib1 = 1; ib2 = 1; 
[lmb1,lmb2,xv1,xv2] = solve_kp88_InGaN(mesh,kx,ky);   %InGaN
a=sqrt(1/2);b=sqrt(1/2);
V0=[b' b  0  0  0  0  0  0;...
    b' -b 0  0  0  0  0  0;...
    0  0  a' 0  0  a  0  0;...
    0  0  0  b  0  0  b' 0;...
    0  0  0  0  b' 0  0  b;...
    0  0  a' 0  0 -a  0  0;...
    0  0  0  b  0  0 -b' 0;...
    0  0  0  0  -b' 0  0 b].';
f2 = V0*reshape(xv2(:,ib2),8,mesh.nn); f4 = f2; 
%f2 = reshape(xv2(:,ib2),8,mesh.nn); f4 = f2;
%
phi = linspace(0,2*pi,31);
Vbare = zeros(1,length(phi));
Vbare_fft = zeros(1,length(phi));
%
for ip = 1:length(phi);
%    
a=exp(1i*3/2*phi(ip))*sqrt(1/2);
b=exp(1i*1/2*phi(ip))*sqrt(1/2);
%
V =[b' b  0  0  0  0  0  0;...
    b' -b 0  0  0  0  0  0;...
    0  0  a' 0  0  a  0  0;...
    0  0  0  b  0  0  b' 0;...
    0  0  0  0  b' 0  0  b;...
    0  0  a' 0  0 -a  0  0;...
    0  0  0  b  0  0 -b' 0;...
    0  0  0  0  -b' 0  0 b].';    
    
% U = f_U(phi(ip), mesh);
% U = trans_wf(phi(ip), mesh.nn, 8);
% xv1_rotated = U*xv1;
% xv2_rotated = U*xv2;
k_rotated = [cos(phi(ip)) sin(phi(ip))]*knorm;
%
 %[lmb1,lmb2,xv1_rotated,xv2_rotated] = solve_kp88_InGaN(mesh,k_rotated(1),k_rotated(2));
%
 %f1 = reshape(xv2_rotated(:,ib1),8,mesh.nn); f3 = f1; %@@@@@@@@@@@@@@@@@@@@@
 f1 = V*reshape(xv2(:,ib1),8,mesh.nn); f3 = f1; %@@@@@@@@@@@@@@@@@@@@@
f14 = sum(conj(f1).*f4,1); 
f23 = sum(conj(f2).*f3,1);
q = k - k_rotated; qnorm = norm(q)*1e10; % 1/m
MM = (f14.'*f23).*exp(-qnorm*abs(yy-xx))';
Vbare(ip) = trapz(mesh.x, trapz(mesh.x,MM,1), 2);
%
%
F14 = fftshift(fft(f14));  
F23 = fftshift(fft(f23));
L = mesh.L;
nh = floor(mesh.nn/2);
k0 = 2*pi/L; ss = - nh:nh;
tmp = 1./(qnorm^2 + ss.^2.*k0^2);
Vbare_fft(ip) = (1/mesh.nn)^2 * (2*L*qnorm*sum(tmp.*F14.*F23(end:-1:1)) + ...
     2*(exp(-qnorm*L)-1)*qnorm^2 *sum(tmp.*F14)*sum(tmp.*F23)  + ...
     2*(exp(-qnorm*L)-1)*   k0^2 *sum(ss.*tmp.*F14)*sum(ss.*tmp.*F23)); end 
%
%
figure(1), hold on
plot(phi,abs(Vbare),'ro')  
plot(phi(1:end-1),abs(Vbare_fft(1:end-1)),'r--') %singular point for fft
grid
%==========================================================================
% test bare coulomb, ordered
%==========================================================================
clear all
%load subbands_AlGaN_ordered.mat;
load subbands_AlGaN_201.mat;
[xx,yy] = meshgrid(mesh.x,mesh.x);
a=sqrt(1/2);b=sqrt(1/2);
V0=[b' b  0  0  0  0  0  0;...
    b' -b 0  0  0  0  0  0;...
    0  0  a' 0  0  a  0  0;...
    0  0  0  b  0  0  b' 0;...
    0  0  0  0  b' 0  0  b;...
    0  0  a' 0  0 -a  0  0;...
    0  0  0  b  0  0 -b' 0;...
    0  0  0  0  -b' 0  0 b].';
f2 = V0*reshape(XV2(:,1,41),8,mesh.nn); f4 = f2;
knorm = abs(kgrid(41)); k = [knorm 0];                   %0.04
%
phi = linspace(0,2*pi,31);
Vbare = zeros(1,length(phi));
for ip = 1:length(phi);
%    
a=exp(1i*3/2*phi(ip))*sqrt(1/2);
b=exp(1i*1/2*phi(ip))*sqrt(1/2);
%
V =[b' b  0  0  0  0  0  0;...
    b' -b 0  0  0  0  0  0;...
    0  0  a' 0  0  a  0  0;...
    0  0  0  b  0  0  b' 0;...
    0  0  0  0  b' 0  0  b;...
    0  0  a' 0  0 -a  0  0;...
    0  0  0  b  0  0 -b' 0;...
    0  0  0  0  -b' 0  0 b].';    
%
k_rotated = [cos(phi(ip)) sin(phi(ip))]*knorm;
f1 = V*reshape(XV2(:,1,41),8,mesh.nn); f3 = f1; 
f14 = sum(conj(f1).*f4,1); 
f23 = sum(conj(f2).*f3,1);
q = k - k_rotated; qnorm = norm(q)*1e10; % 1/m
MM = (f14.'*f23).*exp(-qnorm*abs(yy-xx))';
Vbare(ip) = trapz(mesh.x, trapz(mesh.x,MM,1), 2);
end
%
figure(1), hold on
plot(phi,abs(Vbare),'ro')
grid