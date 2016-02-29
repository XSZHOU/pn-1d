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
ic = 1; iv=4;
%
[lmb1,lmb2,xv1,xv2] = solve_kp88_InGaN(mesh,kx,ky);
fc = reshape(xv1(:,ic),8,mesh.nn);
fv = reshape(xv2(:,iv),8,mesh.nn);
%
[lmb_1,lmb_2,xv_1,xv_2] = solve_kp88_InGaN(mesh,k_rotated(1),k_rotated(2));
f_c = reshape(xv1(:,ic),8,mesh.nn);
f_v = reshape(xv2(:,iv),8,mesh.nn);
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
fc_rotated = V*f_c;
fv_rotated = V*f_v;
%
fc= reshape(fc,1,8*mesh.nn);
fv= reshape(fv,1,8*mesh.nn);
fc_rotated= reshape(fc,1,8*mesh.nn);
fv_rotated= reshape(fv,1,8*mesh.nn);
%
ovp=fc*fv'*mesh.L/mesh.ne
ovp_rotated=fc_rotated*fv_rotated'*mesh.L/mesh.ne
%