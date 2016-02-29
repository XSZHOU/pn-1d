clear all, % close all
clc
load subbands_AlGaN_401.mat   % 401 kpoints
%load subbands_InGaN_201.mat
%
Q=1.6022e-019;        % elementary charge, C
h=6.63e-34;           % Plank constant, J*s
hbar=h/(2*pi);        % reduced Plank constant, J*s
M0=9.10938188E-31;    % electron mass, kg
%
num_kvectors = length(kgrid);              
kx = linspace(0,max(kgrid),num_kvectors);
ky = zeros(1,num_kvectors);
%
 [Hb_te,~,~,Hb_tm] = assem_kp88_AlGaN(1,0,mesh.xmol(1),0);        %H_barrier
 [Hw_te,~,~,Hw_tm] = assem_kp88_AlGaN(1,0,mesh.xmol(mesh.ne/2),0);%H_well
%[Hb_te,~,~,Hb_tm] = assem_kp88_InGaN(1,0,mesh.xmol(1),0);        %H_barrier
%[Hw_te,~,~,Hw_tm] = assem_kp88_InGaN(1,0,mesh.xmol(mesh.ne/2),0);%H_well
ind =find(mesh.xmol==0);dmat=speye(mesh.nn);dmat(ind,ind)=0;
Pte=   (kron(speye(mesh.nn)-dmat,Hw_te)+kron(dmat,Hb_te));       %<|P|>te
Ptm=   (kron(speye(mesh.nn)-dmat,Hw_tm)+kron(dmat,Hb_tm));       %<|P|>tm
Pte=Pte/(hbar/M0)*(mesh.L/mesh.ne);                              % kg*m/s
Ptm=Ptm/(hbar/M0)*(mesh.L/mesh.ne);
M2_te=zeros(1,num_kvectors);M2_tm=zeros(1,num_kvectors);
%
Mp=zeros(1,num_kvectors); Mr=zeros(1,num_kvectors);
MD=spdiags((mesh.x)',0,length(mesh.x),length(mesh.x));
MD=kron(MD,speye(8))*mesh.L/mesh.ne;
%
for ik=1:num_kvectors
    for ic= 1:2%  CB1
        for iv=1:2 % VB1
    M2_te(ik) = M2_te(ik) + abs(XV1(:,ic,ik)'*Pte*XV2(:,iv,ik))^2;
    M2_tm(ik) = M2_tm(ik) + abs(XV1(:,ic,ik)'*Ptm*XV2(:,iv,ik))^2; %(kg*m/s)^2
%
    Mp(ik) = Mp(ik)+hbar/M0*abs(XV1(:,ic,ik)'*Ptm*XV2(:,iv,ik))/(SB1(ic,ik)-SB2(iv,ik));
    Mr(ik) = Mr(ik)+      Q*abs(XV1(:,ic,ik)'*MD *XV2(:,iv,ik)); % C*m
        end
    end
end
figure(1),hold on   % TE vs TM
plot(kgrid,M2_te,'r')
plot(kgrid,M2_tm,'b')
legend('|M|^2_{TE}','|M|^2_{TM}')
%
figure(2),hold on   % 4.93 vs 4.97
plot(kgrid,Mp,'r')
plot(kgrid,Mr,'b')
set(gca,'yscale','log','FontSize',14,'FontName','Times','Box','on')
set(gca,'yscale','lin')
xlabel('kx')
ylabel('oscillator strength')
grid
legend('Momentum matrix element','dipole matrix element')