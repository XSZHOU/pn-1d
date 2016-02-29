%==========================================================================
% 1D Poisson_kp solver,  @Polito S.Z 05/2014
%- coupled with 8-band kp, using full set of WFs to compute charge
%- underrelaxation mixing used to achieve convergence
%- Broyden's second algorithm
%==========================================================================
clear all; clc
%
mate_sys='InGaN';
mesh.nn = 401; mesh.L = 400*1e-10;
mesh.ne=mesh.nn-1;
mesh.x  = linspace(0,mesh.L,mesh.nn);
mesh.le = mesh.x(2:end) - mesh.x(1:end-1);
mesh.xc =(mesh.x(1:end-1) + mesh.x(2:end))/2;
%--------------------------------------------------------------------------specify composition
%mol_barrier=0.54; mol_well=0.34;
mol_barrier=0.0; mol_well=0.25;
mesh.ii = ((mesh.xc>=13e-9)&(mesh.xc<=15e-9))| ...
          ((mesh.xc>=17e-9)&(mesh.xc<=19e-9))| ...
          ((mesh.xc>=21e-9)&(mesh.xc<=23e-9))| ...
          ((mesh.xc>=25e-9)&(mesh.xc<=27e-9));
mesh.xmol= mol_barrier*ones(1,mesh.ne);
mesh.xmol(mesh.ii) = mol_well;
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------load mate_params to mesh
mesh=load_params(mesh,mate_sys);
%--------------------------------------------------------------------------
mesh.nk = 41;       
mesh.ncb=10*4;mesh.nvb=12*4;
esheet=4e13;hsheet=4e13;                                                   
kx = linspace(0,0.3,mesh.nk);
ky = zeros(1,mesh.nk);
%
mesh.evb0 = mesh.evb;                                                      %initial guess for input
unrelax=0.50;                                                              %under-relaxation factor
%--------------------------------------------------------------------------iteration
max_iter=50;res_err=1;iter=1;
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
mesh.target1= min(mesh.evb + mesh.gap);                                    %update target for eigs also
mesh.target2= max(mesh.evb);
%
res_err=norm(mesh.evb-mesh.evb0-phi')/norm(mesh.evb-mesh.evb0);
%
if(iter==max_iter)
   close all;
   figure(1),hold on; 
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
%--------------------------------------------------------------------------check (k,E), targets and Fermi levels
[SB1,SB2]=plot_convergence(mesh,kx,ky,EFn,EFh,mate_sys);
%--------------------------------------------------------------------------