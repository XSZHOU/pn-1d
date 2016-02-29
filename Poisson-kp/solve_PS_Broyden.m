function [eDensity, hDensity, EFn, EFh, mesh, phi_xc] = ...
 solve_PS_Broyden(mesh,kx,ky,esheet,hsheet,iter,unrelax,mate_sys)
%
SB1 = zeros(mesh.ncb,mesh.nk);
SB2 = zeros(mesh.nvb,mesh.nk);
XV1 = zeros(8*mesh.nn,mesh.ncb,mesh.nk);
XV2 = zeros(8*mesh.nn,mesh.nvb,mesh.nk);
%
kgrid = sqrt(kx.^2 + ky.^2);
num_kvectors=length(kgrid);
Q =1.6021917e-19;
%-------------------------------------------------------------------solvekp
for ik = 1:num_kvectors;
fprintf('k, 1/A %e %e\n',[kx(ik) ky(ik)])
%
switch (mate_sys)
case{'AlGaN'}
[lmb1,lmb2,xv1,xv2] = solve_kp88_AlGaN(mesh,kx(ik),ky(ik));
case{'InGaN'}
[lmb1,lmb2,xv1,xv2] = solve_kp88_InGaN(mesh,kx(ik),ky(ik));
otherwise, error('specify correct material system!');
end
%
for ib = 1:mesh.ncb;
xv = xv1(:,ib);
f = reshape(xv,8,mesh.nn);
g = sum(abs(f).^2,1); 
xv = xv / sqrt(trapz(mesh.x,g));
XV1(:,ib,ik) = xv; 
SB1(ib,ik) = lmb1(ib); end 
%
for ib = 1:mesh.nvb;
xv = xv2(:,ib);
f = reshape(xv,8,mesh.nn);
g = sum(abs(f).^2,1); 
xv = xv / sqrt(trapz(mesh.x,g));
XV2(:,ib,ik) = xv; 
SB2(ib,ik) = lmb2(ib); end;
end
%-------------------------------------------------------------------solvekp
%
%------------------------------------OVPcheck,do not use for compute charge
%for
% for ik = 1:(num_kvectors-1);
%     xv0 = XV1(:,:,ik);
%     xv =  XV1(:,:,ik+1);
%     ovp = abs(xv0'*xv);
%     [ovpmax,ii] = max(ovp);
%     XV1(:,ii,ik+1) = XV1(:,:,ik+1); 
%     SB1(ii,ik) = SB1(:,ik); end
% %
% for ik = 1:(num_kvectors-1);
%     xv0 = XV2(:,:,ik);
%     xv =  XV2(:,:,ik+1);
%     ovp = abs(xv0'*xv);
%     [ovpmax,ii] = max(ovp);
%     XV2(:,ii,ik+1) = XV2(:,:,ik+1);
%     SB2(ii,ik) = SB2(:,ik); end
%------------------------------------OVPcheck,do not use for compute charge
%
%--------------------------------------------compute fermi level and charge
OPTIONS = optimset('TolX',1e-9,'MaxIter',1000);
[EFn,FVAL,EXITFLAG] = ...
fminbnd(@(EFn) f_charge_e(mesh,kgrid,SB1,EFn,esheet), 1.5,5.0,OPTIONS);
[EFh,FVAL,EXITFLAG] = ...
fminbnd(@(EFh) f_charge_h(mesh,kgrid,SB2,EFh,hsheet),-1.0,1.0,OPTIONS);
display(['EFn: ', num2str(EFn)]);
display(['EFh: ', num2str(EFh)]);
if(iter==1)
[eDensity,hDensity]=f_charge_evp(mesh,kgrid,SB1,SB2,XV1,XV2,EFn,EFh);
else
[eDensity,hDensity]=f_charge_evp(mesh,kgrid,SB1,SB2,XV1,XV2,EFn,EFh);
end
%--------------------------------------------compute fermi level and charge
%
%-------------------------------------------------------------solve poisson
ind1=1:(mesh.nn-1);ind2=2:mesh.nn;
inr12=[ind1 ind1 ind2 ind2]; inc12=[ind1 ind2 ind1 ind2];
NxNx11=1./mesh.le.*mesh.eps; NxNx12=-1./mesh.le.*mesh.eps; 
NxNx22=1./mesh.le.*mesh.eps;
Amat=sparse(inr12,inc12,[NxNx11 NxNx12 NxNx12 NxNx22],mesh.nn,mesh.nn);
%
% compute e/h
le= mesh.L/mesh.ne;
rho=Q*le*(hDensity-eDensity)+ mesh.piechar;                                %1/le implicit in mesh.piechar due to use of diff
rho(1)=1/2*rho(1);rho(end)=1/2*rho(end);
%rho(1)=0;rho(end)=0;
Amat(1,:)=0;Amat(end,:)=0;                                                 %D.C. in this sense, follows the B.C. of kp solver
Amat(1,1)=1;Amat(end,end)=1;
%phi=Amat\rho;
[LL,UU,PP,QQ]=lu(Amat);
phi= -QQ*(UU\(LL\(PP*rho)));                                               %interpretation -phi! not  - rho! (band-diagram use e-energy)
phi_xc=(phi(1:end-1) + phi(2:end))/2;
%-------------------------------------------------------------solve poisson
%
%------------------------------------------------------------Broyden update
if (iter>1)
    Fn0=mesh.Fn;
end
%mesh.Fn = phi(2:end)-(mesh.evb'-mesh.evb0');                              % 2:end vs 1:end-1
mesh.Fn =phi_xc-(mesh.evb'-mesh.evb0');                                    % electric potential defined on the middle of an element
if(iter==1)
    mesh.G = -unrelax*eye(length(mesh.Fn));                                % approximated inverse Jacobian, "initial guess"
else
    dF=mesh.Fn-Fn0;
    mesh.G = mesh.G + 1/(dF'*dF)* ( mesh.update*dF' -mesh.G*(dF*dF'));
end
mesh.update = -mesh.G * mesh.Fn;
mesh.evb= mesh.evb + mesh.update';
mesh.ecb= mesh.evb;                                                        %update "ecb"
%------------------------------------------------------------Broyden update
end