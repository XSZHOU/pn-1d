%==========================================================================
% 1D Poisson_kp solver,  @Polito S.Z 05/2014
%- coupled with 8-band kp, using full set of WFs to compute charge
%- underrelaxation mixing used to achieve convergence
%- Broyden's second algorithm
%==========================================================================
clear all; clc
%
mate_sys='InGaN';mate_sys='AlGaN';
mesh.nn = 301; mesh.L = 300*1e-10;
mesh.ne=mesh.nn-1;
mesh.x  = linspace(0,mesh.L,mesh.nn);
mesh.le = mesh.x(2:end) - mesh.x(1:end-1);
mesh.xc =(mesh.x(1:end-1) + mesh.x(2:end))/2;
%--------------------------------------------------------------------------specify composition
%mol_barrier=0.54; mol_well=0.34;
%mol_barrier=0.0; mol_well=0.25;
mol_barrier=0.25; mol_well=0.0;
mesh.ii = ((mesh.xc>=10e-9)&(mesh.xc<=12e-9))| ...
          ((mesh.xc>=14e-9)&(mesh.xc<=16e-9))| ...
          ((mesh.xc>=18e-9)&(mesh.xc<=20e-9));
mesh.xmol= mol_barrier*ones(1,mesh.ne);
mesh.xmol(mesh.ii) = mol_well;
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------load mate_params to mesh
mesh=load_params(mesh,mate_sys);
%--------------------------------------------------------------------------
mesh.nk = 41;       
mesh.ncb=8*4;mesh.nvb=8*4;
esheet=2e13;hsheet=2e13;                                                 
kx = linspace(0,0.3,mesh.nk);
ky = zeros(1,mesh.nk);
%
mesh.evb0 = mesh.evb;                                                      %initial guess for input
unrelax=0.5;                                                               %under-relaxation factor
%--------------------------------------------------------------------------iteration
max_iter=60;res_err=1;iter=1;
while res_err>1e-5
%   
display(['iteration: ' num2str(iter)])           
%                                                                          
%[eDensity, hDensity, EFn, EFh, mesh, phi] = ...                            
%solve_PS_linear(mesh,kx,ky,esheet,hsheet,unrelax,mate_sys);               %1e13\2e13\5e13cm^-2   %0.75\0.5\0.15 
%
[eDensity, hDensity, EFn, EFh, mesh, phi] = ...                          
solve_PS_Broyden(mesh,kx,ky,esheet,hsheet,iter,unrelax,mate_sys);          %unrelax here refer to initial guess of inverse J
%
mesh.target1= min(mesh.evb + mesh.gap)+0.5;                                %update target,0.5 for c-subband
mesh.target2= max(mesh.evb);
%
res_err=norm(mesh.evb-mesh.evb0-phi',inf);
%res_err=norm(mesh.evb-mesh.evb0-phi(2:end)')/norm(mesh.evb-mesh.evb0);
%
if(res_err<=1e-5)
   fprintf('convergence achieved with %d iterations !', iter);
   close all;figure(1),hold on;
end
pause(.1)
plot(mesh.xc,mesh.evb,'o-g'), hold on
plot(mesh.xc,mesh.evb+ mesh.gap,'o-g')
plot(mesh.xc,phi'+ mesh.evb0,'r.-','linewidth',1.5)
plot(mesh.xc,phi'+ mesh.evb0 + mesh.gap,'k.-','linewidth',1.5)
plot(mesh.x ,eDensity/max(hDensity)+EFn, 'r-')
plot(mesh.x ,hDensity/max(hDensity)+EFh, 'k-')
axis([0 mesh.L -1 6])
pause(.1)
   if(iter==max_iter)
      fprintf('fail to converge with %d iterations !', max_iter);
      break
   end
iter=iter+1;  
end
%title('In_{0.25}GaN, charge density: 5e^{13} cm^{-2}')
%--------------------------------------------------------------------------check (k,E), targets and Fermi levels
[SB1,SB2]=plot_convergence(mesh,kx,ky,EFn,EFh,mate_sys);
%--------------------------------------------------------------------------