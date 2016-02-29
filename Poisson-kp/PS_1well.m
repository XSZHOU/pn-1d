%==========================================================================
% 1D Poisson_kp solver,  @Polito S.Z 05/2014
%- coupled with 8-band kp, using full set of WFs to compute charge
%- underrelaxation mixing used to achieve convergence
%- Broyden's second algorithm
%==========================================================================
clear all; clc
%
mate_sys='AlGaN';
mesh.nn = 201; mesh.L = 200*1e-10;
mesh.ne=mesh.nn-1;
mesh.x  = linspace(0,mesh.L,mesh.nn);
mesh.le = mesh.x(2:end) - mesh.x(1:end-1);
mesh.xc =(mesh.x(1:end-1) + mesh.x(2:end))/2;
%--------------------------------------------------------------------------specify composition
mol_barrier=0.25; mol_well=0.00;%AlGaN
%mol_barrier=0.0; mol_well=0.18;%InGaN
mesh.ii = ((mesh.xc>=8.5e-9)&(mesh.xc<=11.5e-9));
mesh.xmol= mol_barrier*ones(1,mesh.ne);
mesh.xmol(mesh.ii) = mol_well;
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------load mate_params to mesh
mesh=load_params(mesh,mate_sys);
%--------------------------------------------------------------------------
mesh.nk = 41;       
mesh.ncb=18;mesh.nvb=30;
esheet=1e12;hsheet=1e12;                                                   %*(2e-7)^-1cm^-3, well width-active                                                  
kx = linspace(0,0.3,mesh.nk);
ky = zeros(1,mesh.nk);
%
mesh.evb0 = mesh.evb;                                                      %initial guess for input
unrelax=0.5;                                                               %under-relaxation factor
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
%--------------------------------------------------------------------------check (k,E), targets and Fermi levels
[SB1,SB2]=plot_convergence(mesh,kx,ky,EFn,EFh,mate_sys);
%--------------------------------------------------------------------------
%fft
% d_h = abs(min(mesh.evb0));
% W_well=2e-9;dkz = 2*pi/mesh.L;
% evb0_fft=fftshift(fft(mesh.evb0))/mesh.nn;
% kpoint  = round(mesh.ne/2); qz= -kpoint:kpoint; qz=qz*dkz;
% evb0_ana= d_h * W_well/mesh.L * sin(qz*W_well/2)./(qz*W_well/2);
% figure(2), hold on
% p1=plot(abs(evb0_fft),'ob');
% p2=plot((evb0_ana),'.r');
% %
% % fft, self-consistent evb
% evb_fft=fftshift(fft(mesh.evb))/mesh.nn;
% p3=plot(real(evb_fft),'om');
% legend([p1 p2 p3],'square-fft','square-ana','true-fft')
% title('confining potential: fourier components ')
%