%==========================================================================
%clear all%, close all
clc
%
xmol = 0.25;                                                 
num_kvectors = 101;
%
Kx = [zeros(1,num_kvectors) , linspace(0,0.15,num_kvectors)];
Ky = [zeros(1,num_kvectors), zeros(1,num_kvectors)];
Kz = [linspace(0.15,0,num_kvectors), zeros(1,num_kvectors)];
k =  [-Kz(1:num_kvectors) Kx(num_kvectors+(1:num_kvectors))];
%
EDiagram = zeros(8,2*num_kvectors);
%
for ik = 1:2*num_kvectors;
kx = Kx(ik)*1e10;
ky = Ky(ik)*1e10;
kz = Kz(ik)*1e10;
%
[H0,H1L,H1R,H2] = assem_kp88_AlGaN(kx,ky,xmol);
 Hv = H0 + (H1L+H1R)*kz + H2*kz^2;
[WF, E] = eig(Hv);
[E,Ind] = sort(real(diag(E)),'descend'); WF = WF(:, Ind);
%
EDiagram(:,ik) = E;
end                
%
EDiagram = EDiagram - max(max(EDiagram(3:8,:)));
%
figure(1), hold on
plot(k,EDiagram(1,:),'b.')
plot(k,EDiagram(2,:),'b.')
plot(k,EDiagram(3,:),'g.')
plot(k,EDiagram(4,:),'g.')
plot(k,EDiagram(5,:),'r.')
plot(k,EDiagram(6,:),'r.')
plot(k,EDiagram(7,:),'m.')
plot(k,EDiagram(8,:),'m.')
%
set(gca,'FontSize',14,'FontName','Arial','Box','on')
xlabel('k, 1/A')
ylabel('Energy, eV')
title('Al_{0.25}GaN')
axis([-0.15 0.15 -0.6 5.6])