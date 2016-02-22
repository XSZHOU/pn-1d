%1D Drift-diffusion solver, by S.Z. 03, 2015
%
%Run the code for a biased PN junction
%mesh=pn_1d(1.5e16,1.5e16,0.6);
%
%Poisson discretized by FEM, continuity Eqs dicretized by Finite Volume
%
function mesh = pn_1d(dop_n,dop_p,bias)
q=1.602e-19;T=300;
Kb=1.3807*1e-23;
vt=Kb*T/q;                                         %KbT
epsi=11.7*8.854*1e-14;                             %Si
Eg=1.08;                                           %Si_300K_1.12eV
ni=1.45*1e10;                                      %Si,intrinsic
tau=1e-10;                                         %->recombination const
N_V=1.04e19;
%
mesh.N=501;                                        %discretization
mesh.L=2e-6*100;                                   %2um->2e^-6m->cm
mesh.le=mesh.L/(mesh.N-1);
mesh.Le=mesh.le*ones(1,mesh.N-1);           
V_ref=vt*log(dop_n/ni);                            %V_ref
mesh.EF=vt*log(N_V/dop_p);
N_A=dop_p*[ones(1,(mesh.N-1)/2)  zeros(1,(mesh.N+1)/2)]';     
N_D=dop_n*[zeros(1,(mesh.N-1)/2) ones(1,(mesh.N+1)/2) ]';
%
%------------------------------------------------->Begin Equilibrium case
phi_0=[(V_ref-vt*log(dop_p/ni))*ones((mesh.N-1)/2,1);...
       (V_ref+vt*log(dop_n/ni))*ones((mesh.N+1)/2,1)];
phi=phi_0;
%
in1=1:(mesh.N-1);in2=2:mesh.N;
inr12=[in1 in1 in2 in2]; 
inc12=[in1 in2 in1 in2];
NxNx11=1./mesh.Le; NxNx12=-1./mesh.Le; NxNx22=1./mesh.Le;
Amat=sparse(inr12,inc12,epsi*[NxNx11 NxNx12 NxNx12 NxNx22],mesh.N,mesh.N);
Amat(1,:)=0; Amat(mesh.N,:)=0;                         
    %
    conver_eq=0;iter=0;
    while(~conver_eq)   
    iter=iter+1;
    rho=mesh.le*q*(N_D-N_A-2*ni*sinh((phi-V_ref)/vt));
    E=sparse(1:mesh.N,1:mesh.N, 2*mesh.le*q*ni/vt*cosh((phi-V_ref)/vt) );
    E(1,1)=1;E(mesh.N,mesh.N)=1;
    R_re=-(Amat*phi-rho);
    Jacob=Amat+E;                                      
    [LL,UU,PP,QQ]=lu(Jacob);
    delta_phi=QQ*(UU\(LL\(PP*R_re)));
    %   
    if(norm(delta_phi)/norm(phi)<=1e-9)
    conver_eq=1; end
    %
    phi=phi+delta_phi;
    end
%
mesh.phi=phi'; 
mesh.elecf=-diff(mesh.phi)/mesh.le;
mesh.nn=ni*exp(+(phi-V_ref)/vt);
mesh.np=ni*exp(-(phi-V_ref)/vt);
%rhoabs=abs(-2*ni*sinh((mesh.phi-V_ref)/vt));
rhot=(N_D'-N_A'-2*ni*sinh((mesh.phi-V_ref)/vt));
plot_equi(rhot,Eg,mesh);
%
%-------------------------------------------------->Begin non-Equilibrium
vstep=0.5*vt;ivm=ceil(bias/vstep);vbias=ones(1,ivm);
for iv = 1:ivm
    conver_neq=0;iter=0;                  
    mesh.phi(end)=mesh.phi(end)-vstep;    
    %
    while(~conver_neq)   
    iter=iter+1;
    fprintf('ibias: %d | iter: %d\n',iv,iter)
    mesh.elecf=-diff(mesh.phi)/mesh.le;   
    mesh=comput_mobility(mesh,vt);     
    
    %%%begin Continuity,FVM
    in0=1:mesh.N;in1=1:(mesh.N-1);in2=2:mesh.N;
    inr12=[in0 in1 in2]; 
    inc12=[in0 in2 in1];
    %
    nnw=q*mesh.Dn.*BER1(-diff(mesh.phi)/vt)/mesh.le/mesh.le;
    nne=q*mesh.Dn.*BER1( diff(mesh.phi)/vt)/mesh.le/mesh.le;
    Vw=[nnw(1:end-1) 0]; Ve=[0 nne(2:end)];
    nnc=-(nnw(2:end)+nne(1:end-1)); Vc=[1 nnc 1];
    Nmat=sparse(inr12,inc12,[Vc Ve Vw],mesh.N,mesh.N);
    %
    npw=q*mesh.Dp.*BER1( diff(mesh.phi)/vt)/mesh.le/mesh.le;
    npe=q*mesh.Dp.*BER1(-diff(mesh.phi)/vt)/mesh.le/mesh.le;
    Vw=[npw(1:end-1) 0]; Ve=[0 npe(2:end)];
    npc=-(npw(2:end)+npe(1:end-1)); Vc=[1 npc 1];
    Pmat=sparse(inr12,inc12,[Vc Ve Vw],mesh.N,mesh.N);
    %
    Gn= q*(mesh.nn.*mesh.np-ni*ni)./(tau*(mesh.nn+ni)+tau*(mesh.np+ni));
    Gp= q*(mesh.nn.*mesh.np-ni*ni)./(tau*(mesh.nn+ni)+tau*(mesh.np+ni));
    Gn(1)=mesh.nn(1);Gn(end)=mesh.nn(end);
    Gp(1)=mesh.np(1);Gp(end)=mesh.np(end);
    mesh.nn=Nmat\Gn;                                                      
    mesh.np=Pmat\Gp;                                                      
    
    %%%begin Poisson,FEM/FDM
    Apmat=Amat+ sparse(1:mesh.N,1:mesh.N,...
          q*mesh.le/vt*(mesh.nn+mesh.np) );
    bp=mesh.le*q*(N_D-N_A+mesh.np-mesh.nn)+...
          q*mesh.le/vt*(mesh.nn+mesh.np).*mesh.phi'; 
    mesh.oldphi=mesh.phi;
    Apmat(1,1)=1;Apmat(mesh.N,mesh.N)=1;
    bp(1)=mesh.phi(1);bp(end)=mesh.phi(end);
    phi=Apmat\bp;delta_phi=phi'-mesh.oldphi;                              
    %
    if(norm(delta_phi)/norm(mesh.phi)<=1e-6)
    conver_neq=1;vbias(iv)=vstep*iv;
    fprintf('converged! bias=%f V\n',vbias(iv));
    end
    mesh.phi=phi';
   end
    
    %%%elemental current
    Jn1=q/mesh.le*mesh.Dn.*(BER1( diff(mesh.phi)/vt).*mesh.nn(2:end)'-...
                            BER1(-diff(mesh.phi)/vt).*mesh.nn(1:end-1)');
    Jp1=q/mesh.le*mesh.Dp.*(BER1(-diff(mesh.phi)/vt).*mesh.np(2:end)'-...
                            BER1( diff(mesh.phi)/vt).*mesh.np(1:end-1)');    
    mesh.Jnx=Jn1;
    mesh.Jpx=Jp1;
    mesh.Jtot=mesh.Jnx-mesh.Jpx; 
    mesh.Jv(iv)=mesh.Jtot(1);
    mesh.Jn(iv)=mesh.Jnx(1);
    mesh.Jp(iv)=mesh.Jpx(1);

end   %------------------------------------------->End Non-Equilibrium

%--------------------------------------------------------->End main job!

mesh.EFn=-mesh.phi+vt*log(mesh.nn'/ni);
mesh.EFp=-mesh.phi-vt*log(mesh.np'/ni);
rhot=N_D-N_A+mesh.np-mesh.nn;
plot_nonequi(vbias,rhot,Eg,mesh);
% figure(3),hold on   %I-V
% plot(mesh.Jv,'-g')
% plot(mesh.Jn,'-.r')
% plot(mesh.Jp,'-.b')
figure(3),hold on
plot(mesh.Jnx,'-r') %I-x
plot(-mesh.Jpx,'-b')
plot(mesh.Jtot,'-.g')
legend('Jn','Jp','Jtot')
title('Converged current,A/cm^2','fontweight','bold')

%---------------------------------------------------------------------%

function mesh = comput_mobility(mesh,vt)
% Mobility: Thomas model | ATLAS: VSAT = (2.4*10^7)/(1+0.8*exp(T/600))
vsatn=1.07e7;vsatp=8.37e6;betan=2; betap=1; %cm/s
mobn0=960;mobp0=435;% ~2e16 doping level dependent,cm^2/V/sec
%mobn0=1500;mobp0=1000;
mesh.mobn=mobn0./((1+(mobn0.*abs(mesh.elecf)/vsatn).^betan).^(1/betan));
mesh.mobp=mobp0./((1+(mobp0.*abs(mesh.elecf)/vsatp).^betap).^(1/betap));
%mesh.mobn(1)=mesh.mobn(2);mesh.mobn(end)=mesh.mobn(end-1);
%mesh.mobp(1)=mesh.mobp(2);mesh.mobp(end)=mesh.mobp(end-1);
mesh.Dn=mesh.mobn*vt;
mesh.Dp=mesh.mobp*vt;

%---------------------------------------------------------------------%

function B=BER1(x)
% Cf. book by Selberherr(1984), p.169.
% The parameters x1,...,x5 are evaluated for MATLAB.
% Defines the nodes for the approximation
x1=-36.25; x2=-7.63e-6; x3=-x2; x4=32.92; x5=36.5;
%
B=zeros(size(x));
%
B1 = (x<=x1); 
B(B1) = -x(B1);
%
B2= (x>x1)&(x<x2);
B(B2) = +  x(B2)./(exp(x(B2))-1+1.e-99);
%
B3 = (x>=x2)&(x<=x3);
B(B3) = 1-x(B3)/2;
%
B4 = (x>x3)&(x<x4);
B(B4) = x(B4).*exp(-x(B4))./(1-exp(-x(B4))+1.e-99);
%
B5= (x>=x4)&(x<x5);
B(B5) =  x(B5).*exp(-x(B5));
%
B6= (x>=x5);
B(B6) = 0.0;

%---------------------------------------------------------------------%

function [] = plot_nonequi(vbias,rhot,Eg,mesh)
x_ratio=1e-6/(mesh.le/100);   %1e-6->um|le/100:cm->um
N=mesh.N;
%
figure(2),
%
subplot(2,2,1),
plot(vbias,mesh.Jv,'b','linewidth',1.5);
set(gca,'yscale','log')
xlabel('bias voltage, V','fontweight','bold')
ylabel('total current, A/cm^2','fontweight','bold')
title('Forward biased PN: IV curve','fontweight','bold')
%
subplot(2,2,2)
plot(mesh.elecf,'b','linewidth',1.5)
xlabel('position in 1D,um','fontweight','bold')
ylabel('intensity V/cm','fontweight','bold')
title('electric field,V/cm','fontweight','bold')
set(gca,'xlim',[0,N],'Xtick',0:(N-1)/4:N-1);
set(gca,'xticklabel',str2double(get(gca,'XTickLabel'))/x_ratio);
%
subplot(2,2,3),hold on
plot(mesh.nn,'g','linewidth',1.5)
plot(mesh.np,'b','linewidth',1.5)
plot(rhot,'--r','linewidth',1.5)
legend('electron','hole','total charge')
xlabel('position in 1D,um','fontweight','bold')
ylabel('concerntration cm_-3 (log)','fontweight','bold')
title('carrier density','fontweight','bold')
%set(gca,'yscale','log');ylim([1e3 1e18]);
set(gca,'xlim',[0,N],'Xtick',0:(N-1)/4:N-1);
set(gca,'xticklabel',str2double(get(gca,'XTickLabel'))/x_ratio);
%
subplot(2,2,4),hold on
p1=plot(-mesh.phi+Eg,'m','linewidth',1.5);
p2=plot(-mesh.phi   ,'m','linewidth',1.5);
p3=plot(-mesh.phi+Eg/2,'--k','linewidth',1.5);%approximately meffe!=meffh
p4=plot( mesh.EFn+Eg/2,'r','linewidth',1.5);
p5=plot( mesh.EFp+Eg/2,'b','linewidth',1.5);
legend([p1 p2 p3 p4 p5],'Ec','Ev','EFi','EFn','EFp')
xlabel('position in 1D,um','fontweight','bold')
ylabel('Energy level eV','fontweight','bold')
title('band diagram','fontweight','bold')
ylim([-.5 1.5]);
set(gca,'xlim',[0,N],'Xtick',0:(N-1)/4:N-1);
set(gca,'xticklabel',str2double(get(gca,'XTickLabel'))/x_ratio);

%---------------------------------------------------------------------%

function [] = plot_equi(rhot,Eg,mesh)
x_ratio=1e-6/(mesh.le/100);   %1e-6->um|le/100:cm->um
N=mesh.N;
%
figure(1),
%
subplot(2,2,1)
plot(mesh.phi,'b','linewidth',1.5)
xlabel('position in 1D,um','fontweight','bold')
ylabel('potential V','fontweight','bold')
title('Equilibrium PN: Potential','fontweight','bold')
set(gca,'xlim',[0,N],'Xtick',0:(N-1)/4:N-1);
set(gca,'xticklabel',str2double(get(gca,'XTickLabel'))/x_ratio);
%
subplot(2,2,2)
plot(mesh.elecf,'b','linewidth',1.5)
xlabel('position in 1D,um','fontweight','bold')
ylabel('intensity V/cm','fontweight','bold')
title('electric field,V/cm','fontweight','bold')
set(gca,'xlim',[0,N],'Xtick',0:(N-1)/4:N-1);
set(gca,'xticklabel',str2double(get(gca,'XTickLabel'))/x_ratio);
%
subplot(2,2,3),hold on
%plot(rhoabs,'b','linewidth',1.5)
plot(mesh.nn,'g','linewidth',1.5)
plot(mesh.np,'b','linewidth',1.5)
plot(rhot,'--r','linewidth',1)
legend('electron','hole','total charge')
xlabel('position in 1D,um','fontweight','bold')
ylabel('concerntration cm_-3 (log)','fontweight','bold')
title('carrier density','fontweight','bold')
set(gca,'xlim',[0,N],'Xtick',0:(N-1)/4:N-1);
set(gca,'xticklabel',str2double(get(gca,'XTickLabel'))/x_ratio);
%
subplot(2,2,4),hold on
p1=plot(-mesh.phi,'m','linewidth',1.5);
p2=plot(-mesh.phi+Eg,'b','linewidth',1.5);
p3=plot([0 N],[mesh.EF mesh.EF],'--r','linewidth',1.5);
legend([p1 p2 p3],'Ec','Ev','EF')
xlabel('position in 1D,um','fontweight','bold')
ylabel('Energy level eV','fontweight','bold')
title('band diagram','fontweight','bold')
ylim([-1 1.5]);
set(gca,'xlim',[0,N],'Xtick',0:(N-1)/4:N-1);
set(gca,'xticklabel',str2double(get(gca,'XTickLabel'))/x_ratio);