%==========================================================================
function [lmb1,lmb2,xv1,xv2]=solve_kp88_AlGaN(mesh,kx,ky,flag)
%==========================================================================
%
if (nargin<4); p_asy=1;
elseif strcmp(flag,'symm'); p_asy=0.5; 
else   p_asy=flag;
end
%
kx = kx*1e10; % 1/m              
ky = ky*1e10; % 1/m
nn = mesh.nn;
NPlaneWaves = 8; % number of basis functions
%
I = eye(NPlaneWaves);
I(1:2,1:2) = 0;
Iv = I;
I = eye(NPlaneWaves);
I(3:8,3:8)=0;
Ic = I;
I = eye(NPlaneWaves);  
%
Amat=spalloc(nn*NPlaneWaves,nn*NPlaneWaves,NPlaneWaves^2*(3*nn-2));
Bmat=spalloc(nn*NPlaneWaves,nn*NPlaneWaves,NPlaneWaves*(3*nn-2));
%
for ie=1:mesh.ne;
%
ecb = mesh.ecb(ie);%12/05/2014
evb = mesh.evb(ie);
xmol= mesh.xmol(ie);
le =  mesh.le(ie);
%
[H0,H1l,H1r,H2] = assem_kp88_AlGaN(kx,ky,xmol);
%[H0,H1l,H1r,H2] = assem_kpdirect88_AlGaN(kx,ky,xmol);
H1L=H1l*p_asy+H1r*(1-p_asy);
H1R=H1r*p_asy+H1l*(1-p_asy);
%
NN11=1/3*le;  NN12=1/6*le;   NN22=1/3*le;
NxNx11=1./le; NxNx12=-1./le; NxNx22=1./le;
NxN = 1i/2*[-1 -1; 1  1]; 
NNx = NxN';
% node indices
i1 = (1:NPlaneWaves) + (ie-1)*NPlaneWaves;
i2 = (1:NPlaneWaves) + ie*NPlaneWaves;
%   
Amat(i1,i1) = Amat(i1,i1) + H2*NxNx11;
Amat(i1,i2) = Amat(i1,i2) + H2*NxNx12;
Amat(i2,i1) = Amat(i2,i1) + H2*NxNx12;
Amat(i2,i2) = Amat(i2,i2) + H2*NxNx22;
%
Amat(i1,i1) = Amat(i1,i1) + H1L*NxN(1,1);
Amat(i1,i2) = Amat(i1,i2) + H1L*NxN(1,2);
Amat(i2,i1) = Amat(i2,i1) + H1L*NxN(2,1);
Amat(i2,i2) = Amat(i2,i2) + H1L*NxN(2,2);
%   
Amat(i1,i1) = Amat(i1,i1) + H1R*NNx(1,1);
Amat(i1,i2) = Amat(i1,i2) + H1R*NNx(1,2);
Amat(i2,i1) = Amat(i2,i1) + H1R*NNx(2,1);
Amat(i2,i2) = Amat(i2,i2) + H1R*NNx(2,2);
%   
Amat(i1,i1) = Amat(i1,i1) + (H0+Iv*evb+Ic*ecb)*NN11;
Amat(i1,i2) = Amat(i1,i2) + (H0+Iv*evb+Ic*ecb)*NN12;
Amat(i2,i1) = Amat(i2,i1) + (H0+Iv*evb+Ic*ecb)*NN12;
Amat(i2,i2) = Amat(i2,i2) + (H0+Iv*evb+Ic*ecb)*NN22;
%   
Bmat(i1,i1) = Bmat(i1,i1) + I*NN11;
Bmat(i1,i2) = Bmat(i1,i2) + I*NN12;
Bmat(i2,i1) = Bmat(i2,i1) + I*NN12;
Bmat(i2,i2) = Bmat(i2,i2) + I*NN22;
end
%
opts.disp=0; % verbose mode
target1 = mesh.target1;                    
target2 = mesh.target2;           
%
[xv1,lmb1,flag] = eigs(Amat,Bmat,mesh.ncb,target1,opts);
[lmb1,ind] = sort(diag(real(lmb1)),'ascend'); xv1 = xv1(:,ind);
%
[xv2,lmb2,flag] = eigs(Amat,Bmat,mesh.nvb,target2,opts);
[lmb2,ind] = sort(diag(real(lmb2)),'descend'); xv2 = xv2(:,ind);