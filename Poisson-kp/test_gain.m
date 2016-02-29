clear all, % close all
clc
load subbands_AlGaN_401.mat   % 401 kpoints
%
Q=1.6022e-019;        % elementary charge, C
h=6.63e-34;           % Plank constant, J*s
hbar=h/(2*pi);        % reduced Plank constant, J*s
M0=9.10938188e-31;    % electron mass, kg
eps0=8.854*1e-12;
c=2.99792458e8;
kB = 1.38e-23; 
T = 300;
%
n_GaN=sqrt(10.39);
gama=1/(0.1*1e-12);                             %(0.1ps)^-1
omega=linspace(3.55/hbar*Q,3.85/hbar*Q,101);      % eV
esheet=1e13; hsheet=1e13;                       % sheet density cm^-2    
%
num_kvectors = length(kgrid);              
kx = linspace(0,max(kgrid),num_kvectors);
ky = zeros(1,num_kvectors);
[Hb_te,~,~,Hb_tm] = assem_kp88_AlGaN(1,0,mesh.xmol(1),0);        %H_barrier
[Hw_te,~,~,Hw_tm] = assem_kp88_AlGaN(1,0,mesh.xmol(mesh.ne/2),0);%H_well
ind =find(mesh.xmol==0);dmat=speye(mesh.nn);dmat(ind,ind)=0;     %struc ind
Pte=abs(kron(speye(mesh.nn)-dmat,Hw_te)+kron(dmat,Hb_te));       %<|P|>te
Ptm=abs(kron(speye(mesh.nn)-dmat,Hw_tm)+kron(dmat,Hb_tm));       %<|P|>tm
Pte=Pte/(hbar/M0)*(mesh.L/mesh.ne);                              % kg*m/s
Ptm=Ptm/(hbar/M0)*(mesh.L/mesh.ne);                              % kg*m/s
M_te=zeros(1,num_kvectors);M_tm=zeros(1,num_kvectors);           % |M|^2        
guass=zeros(1,num_kvectors);                                 %stat factor
Gsp_te=zeros(1,length(omega));Gsp_tm=zeros(1,length(omega)); %sp emiss rate
Gm_te=zeros(1,length(omega));Gm_tm=zeros(1,length(omega));   %material gain
%
OPTIONS = optimset('TolX',1e-9,'MaxIter',1000);
[EFn,FVAL,EXITFLAG] = ...
fminbnd(@(EFn) f_charge_e(mesh,kgrid,SB1,EFn,esheet), 3,4,OPTIONS);
[EFh,FVAL,EXITFLAG] = ...
fminbnd(@(EFh) f_charge_h(mesh,kgrid,SB2,EFh,hsheet),-0.5,0.5,OPTIONS);
kgrid = kgrid*1e10;
%
for ie= 1:length(omega)
    ie
    for ic= 1:4%  CB
        for iv=1:8 % VB
            for ik=1:num_kvectors
              fc = 1.0D0./(1.0D0 + exp((SB1(ic,ik)-EFn)/(kB*T/Q)));
              fv = 1.0D0./(1.0D0 + exp((SB2(iv,ik)-EFh)/(kB*T/Q)));
              guass(ik)= fc*(1-fv)*(hbar*gama/pi)/...
              (((SB1(ic,ik)-SB2(iv,ik))*Q-hbar*omega(ie))^2+(hbar*gama)^2);
              M_te(ik) = abs(XV1(:,ic,ik)'*Pte*XV2(:,iv,ik))^2; %(kg*m/s)^2
              M_tm(ik) = abs(XV1(:,ic,ik)'*Ptm*XV2(:,iv,ik))^2; %(kg*m/s)^2 
            end
         Gsp_te(ie)=Gsp_te(ie)+trapz(kgrid,guass.*M_te.*kgrid/2/pi);
         Gsp_tm(ie)=Gsp_tm(ie)+trapz(kgrid,guass.*M_tm.*kgrid/2/pi);
        end
    end
end
%
mesh.Lz=mesh.L*length(ind)/mesh.ne;  
Gsp_te=Gsp_te*Q^2*pi/n_GaN/c/eps0/M0^2./omega/(mesh.Lz*1e2); %sp emiss rate
Gsp_tm=Gsp_tm*Q^2*pi/n_GaN/c/eps0/M0^2./omega/(mesh.Lz*1e2); %sp emiss rate         
Gte=Gsp_te.*(1-exp((hbar*omega/Q-(EFn-EFh))/(kB*T/Q)));      %material gain
Gtm=Gsp_tm.*(1-exp((hbar*omega/Q-(EFn-EFh))/(kB*T/Q)));      %material gain
%
figure(1),hold on   % material gain, TE vs TM
plot(hbar*omega/Q,Gte,'r')                                   %1/cm
plot(hbar*omega/Q,Gtm,'b')                                   %1/cm
ylim([-2000 10000])
%title('Gain: sheet density from 5e^{12} to 9e^{12} cm^{-2}');
%arrowobj=annotation('arrow',[0.4 0.8],[0.2 0.6],'Color','g')