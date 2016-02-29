%==========================================================================
function [SXE]=f_SXE(mesh,kgrid,SB1,SB2,XV1,XV2,EFn,EFh,ib2)
%==========================================================================
%
kB = 1.38D-23;      % Boltzmann constant, J/K
Q = 1.6021917D-19;  % electron charge
T = 300;            % K 
E0 = 8.854D-12;     % F/m
ER = 10.38;         % 1999Chow_PRB
%
kgrid = kgrid*1e10; % 1/m
[xx,yy] = meshgrid(mesh.x,mesh.x);
nk = length(kgrid); np = 101;
phi = linspace(0,pi,np);
Vp  = zeros(1,np);
Vk  = zeros(1,nk);
nk1 = zeros(1,nk);
Qk = zeros(1,nk);
SXE = NaN*zeros(1,nk);
%
a=sqrt(1/2);b=sqrt(1/2);
U0=[b' b  0  0  0  0  0  0;...
    b' -b 0  0  0  0  0  0;...
    0  0  a' 0  0  a  0  0;...
    0  0  0  b  0  0  b' 0;...
    0  0  0  0  b' 0  0  b;...
    0  0  a' 0  0 -a  0  0;...
    0  0  0  b  0  0 -b' 0;...
    0  0  0  0  -b' 0  0 b].'; 
%
for ik2 = 1:10:401;   % loop on k,inital
k2 = [kgrid(ik2) 0];  % 1/m
f2 = U0*reshape(XV2(:,ib2,ik2),8,mesh.nn); f4 = f2;
%==========================================================================
%
  for ik = 1:nk;      % loop on k+q,final
      ib1 = ib2; k1 = [kgrid(ik) 0]; k1norm = norm(k1); %  1/m
      f1_tmp = reshape(XV2(:,ib1,ik),8,mesh.nn); 
      E1 =  SB2(ib1,ik);
      nk1(ik) = 1.0D0./(1.0D0 + exp(-(E1-EFh)/(kB*T/Q)));
%
       for ip = 1:np; % loop on phi
%    
a=exp(1i*3/2*phi(ip))*sqrt(1/2);
b=exp(1i*1/2*phi(ip))*sqrt(1/2);
U =[b' b  0  0  0  0  0  0;...
    b' -b 0  0  0  0  0  0;...
    0  0  a' 0  0  a  0  0;...
    0  0  0  b  0  0  b' 0;...
    0  0  0  0  b' 0  0  b;...
    0  0  a' 0  0 -a  0  0;...
    0  0  0  b  0  0 -b' 0;...
    0  0  0  0  -b' 0  0 b].';    
%
k1_rotated = [cos(phi(ip)) sin(phi(ip))]*k1norm;
%
f1 = U*f1_tmp; f3 = f1;
f14 = sum(conj(f1).*f4,1); 
f23 = sum(conj(f2).*f3,1);
q = k1_rotated - k2; qnorm = norm(q);% 1/m
%MM = (f14.'*f23).*exp(-qnorm*abs(yy-xx))';
%Vp(ip) = Q^2/(2*E0*ER*qnorm) * trapz(mesh.x, trapz(mesh.x,MM,1), 2);
Vp(ip)=Q^2/(2*E0*ER*qnorm)*(mesh.L/mesh.ne)^2*f14*exp(-qnorm*abs(yy-xx))'*f23.';
if(qnorm<1), Vp(ip) = 0; end  
        end
%
   Vk(ik) = trapz(phi,Vp);  Qk(ik)=qnorm;
   end
%          2 >   [0,pi]
SXE(ik2) = 2 * 1/(2*pi)^2*trapz(kgrid,Vk.*nk1.*kgrid)/Q; % eV
fprintf('SXE (eV): %e\n',SXE(ik2))
end
%==========================================================================