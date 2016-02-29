function [H0,H1L,H1R,H2] = assem_kp88_AlGaN( kx,ky,xmol,flag)%flag for H_te
%==========================================================================
Q=1.6022e-019;        % elementary charge, C 
h=6.63e-34;           % Plank constant, J*s
hbar=h/(2*pi);        % reduced Plank constant, J*s
M0=9.10938188E-31;    % electron mass, kg
H2M0 = hbar^2/(2*M0); % (J*s)^2/Kg/m^2
s2=sqrt(2);
%==========================================================================
% GaN
%============================================================8by8parameters
Delta1 = 0.01*Q;
Delta2 = 0.017/3*Q;
Delta3 = 0.017/3*Q;
mz=0.186*M0;mt=0.209*M0;                                           
Eg1=3.437*Q;                                                               %Eg1->fundamental gap, for GaN/InN
Eg=Eg1+Delta1+Delta2;                                                      %Eg ->cb+Eg - hh = fundamental gap  
P11=H2M0*(M0/mz-1)/(Eg1+2*Delta2)*((Eg1+Delta1+Delta2)*...
    (Eg1+2*Delta2)-2*Delta3^2);
P22=H2M0*(M0/mt-1)*Eg1*((Eg1+Delta1+Delta2)*(Eg1+2*Delta2)-2*Delta3^2)/...
    ((Eg1+Delta1+Delta2)*(Eg1+Delta2)-Delta3^2);
P1=sqrt(P11);
P2=sqrt(P22);
%effectively decoupled
%P1=0;P11=0; 
%P2=0;P22=0;
ax=H2M0 * ((1/mt)*M0-P22/H2M0/Eg1);
ay=ax;
az=H2M0 * ((1/mz)*M0-P11/H2M0/Eg);                              %25/10/2013
%============================================================1996Chuang_PRB
D1 = -3.6*Q;
D2 =  1.7*Q; 
D3 =  5.2*Q;
D4 = -2.7*Q;
D5 = -2.8*Q; 
D6 = -4.3*Q;
%
A1 = -7.21+P11/H2M0/Eg;
A2 = -0.44; 
A3 =  6.68-P11/H2M0/Eg;
A4 = -3.46+P22/H2M0/Eg1/2;
A5 = -3.40+P22/H2M0/Eg1/2;
A6 = -4.90+P1*P2/H2M0/sqrt(Eg1*Eg)/s2;
A7 = 0.0937*1e-10*Q;%A7=0;
%
a0 = 3.189;
a  = a0;
%
c0 = 5.185;
c  = c0;
%
C11 = 390;
C12 = 145;
C13 = 106;
C33 = 398;
d13 = -1;
%
exx = (a0-a)/a;       %a0=a_sub, -0.01 compressive
eyy = exx;
ezz = -2*C13/C33*exx;
exy = 0;
exz = 0;
eyz = 0;
%
Delta1 = 0.01*Q;
Delta2 = 0.017/3*Q;
Delta3 = 0.017/3*Q;
%
alphaxx=A2+A4;
alphayy=alphaxx;
alphazz=A1+A3;
alphaxz=0;
alphazx=alphaxz;
betaxx=-A5;                              
betayy=-betaxx;
betaxy=2*A5;
betayz=0;
betazy=betayz;
betazz=0;
gamaxz=-A6;                                %Veprek's form: A5, A6 have more 
gamazx=0;                                   %24/11/2012
gamayz= A6;                                 %02/04/2013
gamazy=0;                                   %02/04/2013
gamaxx=0;
gamayy=0;
gamazz=0;
gamaxy=0;
gamayx=0;
muxx=A2;
muyy=muxx;
muzz=A1;
muxz=0;
muzx=muxz;
%
l1=D2+D4+D5;
l2=D1;
m1=D2+D4-D5;
m2=D1+D3;
m3=D2;
n1=-2*D5;
n2=sqrt(2)*D6;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=alphaxx*kx^2+alphayy*ky^2;
h12=betaxx* kx^2+betayy *ky^2+1i*betaxy*kx*ky; 
h13=gamaxx* kx^2+gamayy *ky^2+1i*gamaxy*kx*ky;
h33=muxx  * kx^2+muyy   *ky^2;
h22= h11';
h23=-h13';
h21=h12';
h31=h13';
h32=h23';
%
H0_GaN = [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
%kane's band edge
H0_GaNE= [Delta1+Delta2 0 0  0      0           0; ...
          0 Delta1-Delta2 0  0      0   s2*Delta3; ...
          0       0       0  0  s2*Delta3       0; ...
          0       0       0  Delta1+Delta2  0   0; ...
          0   0   s2*Delta3  0  Delta1-Delta2   0; ...
          0   s2*Delta3   0  0        0         0];
%strain induced
H0_GaND= [l1*exx+m1*eyy+m2*ezz  n1* exy  n2*exz  0  0  0; ...
          n1*exy m1*exx+l1*eyy+m2*ezz    n2*eyz  0  0  0; ...
          n2*exz  n2*eyz  m3*exx+m3*eyy+l2*ezz   0  0  0; ...
          0  0  0  l1*exx+m1*eyy+m2*ezz  n1* exy  n2*exz; ...
          0  0  0  n1*exy m1*exx+l1*eyy+m2*ezz    n2*eyz; ...
          0  0  0  n2*exz  n2*eyz  m3*exx+m3*eyy+l2*ezz ];
%include the A7 parameter
H0_GaNA7= [0  0     1i*A7*(kx-1i*ky)    0   0   0;...
           0  0    -1i*A7*(kx+1i*ky)    0   0   0;...
          -1i*A7*(kx+1i*ky) 1i*A7*(kx-1i*ky) 0 0 0 0;...
          0 0 0    0  0   -1i*A7*(kx+1i*ky);...
          0 0 0    0  0    1i*A7*(kx-1i*ky);...
          0 0 0    1i*A7*(kx-1i*ky)  -1i*A7*(kx+1i*ky) 0]; 
H0_GaN =H0_GaN+H0_GaND+H0_GaNE+H0_GaNA7;
Hcc=[kx^2*ax+ky^2*ay+Eg 0;0 kx^2*ax+ky^2*ay+Eg];
Hcv=[-P2*(kx+1i*ky)/s2 P2*(kx-1i*ky)/s2 0 0 0 0;...
      0 0 0 P2*(kx-1i*ky)/s2 -P2*(kx+1i*ky)/s2 0];
Hvc=Hcv';
H0_GaN=[Hcc Hcv; Hvc H0_GaN];
% assemble only when computing the polarization
if (nargin>3)
H0_GaNte = [zeros(2) Hcv; Hvc H0_GaNA7]/sqrt(kx^2+ky^2);
Hcv=[0 0 P1 0 0 0;0 0 0 0 0 P1]; Hvc=Hcv'; 
H0_GaNtm = [zeros(2) Hcv; Hvc zeros(6)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=alphazx*kx;   
h12=1i*betazy*ky;
h13=0;            
h33=muzx*kx;      
h22= h11' ;       
h23=0;           
h21=-1i*ky*betayz;          
h31= gamaxz'*kx-1i*ky*gamayz';
h32=-gamaxz'*kx-1i*ky*gamayz';
%                                                       kz*xxx
H1_GaNL= [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
%                                                     
h11=alphaxz*kx;  
h12=1i*betayz*ky; 
h13=gamaxz*kx+1i*gamayz*ky;
h33=muxz*kx;      
h22=h11' ;       
h23=-gamaxz*kx+1i*gamayz*ky;
h21=-1i*ky*betazy;
h31=0;
h32=0;      
%                                                       xxx*kz
H1_GaNR= [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
%
Hcc=[0 0;0 0];
Hcv=[0 0 P1 0 0 0;...
     0 0 0 0 0 P1]*0.5;    
H1_GaNL=[Hcc Hcv; Hcv' H1_GaNL];
H1_GaNR=[Hcc Hcv; Hcv' H1_GaNR];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=alphazz;
h12=betazz; 
h13=gamazz; 
h33=muzz;
h22=h11';
h23=-h13';  
h21=h12';   
h31=h13';   
h32=h23';   
%
H2_GaN = [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
Hcc=[az 0;0 az];
Hcv=[0 0 0 0 0 0;...
     0 0 0 0 0 0];
Hvc=Hcv';
H2_GaN=[Hcc Hcv; Hvc H2_GaN];
%==========================================================================
% AlN
%============================================================8by8parameters
Delta1 = -0.227*Q;
Delta2 = 0.036/3*Q;
Delta3 = 0.036/3*Q;
mz=0.30*M0;mt=0.32*M0;
Eg1=6.00*Q;                                                                %Eg1->fundamental gap, for AlN
Eg=Eg1+Delta1+Delta2;                                                      %Eg ->cb+Eg -so = fundamental gap
P11=H2M0*(M0/mz-1)/(Eg1+2*Delta2)*((Eg1+Delta1+Delta2)*...
    (Eg1+2*Delta2)-2*Delta3^2);
P22=H2M0*(M0/mt-1)*Eg1*((Eg1+Delta1+Delta2)*(Eg1+2*Delta2)-2*Delta3^2)/...
    ((Eg1+Delta1+Delta2)*(Eg1+Delta2)-Delta3^2);
P1=sqrt(P11);
P2=sqrt(P22);
%effectively decoupled
%P1=0;P11=0;
%P2=0;P22=0;
ax=H2M0 * ((1/mt)*M0-P22/H2M0/Eg1);
ay=ax;
az=H2M0 * ((1/mz)*M0-P11/H2M0/Eg);
%============================================================1996Chuang_PRB
D1 = -2.9*Q;
D2 =  4.9*Q; 
D3 =  9.4*Q;
D4 = -4.0*Q;
D5 = -3.3*Q; 
D6 = -2.7*Q;
%
A1 = -3.86+P11/H2M0/Eg;
A2 = -0.25; 
A3 =  3.58-P11/H2M0/Eg;
A4 = -1.32+P22/H2M0/Eg1/2;
A5 = -1.47+P22/H2M0/Eg1/2;
A6 = -1.64+P1*P2/H2M0/sqrt(Eg1*Eg)/s2;
A7 =  0*1e-10*Q;
%
a0 = 3.112;
a  = a0;
%
c0 = 4.982;
c  = c0;
%
C13 = 108;
C33 = 373; 
%
exx = (a0-a)/a;
eyy = exx;
ezz = -2*C13/C33*exx;
exy = 0;
exz = 0;
eyz = 0;
%
alphaxx=A2+A4;
alphayy=alphaxx;
alphazz=A1+A3;
alphaxz=0;
alphazx=alphaxz;
betaxx=-A5;                       
betayy=-betaxx;
betaxy=2*A5;
betayz=0;
betazy=betayz;
betazz=0;
gamaxz=-A6;                                                  %Veprek's form
gamazx=0;                                                       %24/11/2012
gamayz= A6;                                                     %02/04/2013
gamazy=0;                                                       %02/04/2013
gamaxx=0;
gamayy=0;
gamazz=0;
gamaxy=0;
gamayx=0;
muxx=A2;
muyy=muxx;
muzz=A1;
muxz=0;
muzx=muxz;
%
l1=D2+D4+D5;
l2=D1;
m1=D2+D4-D5;
m2=D1+D3;
m3=D2;
n1=-2*D5;
n2=sqrt(2)*D6;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=alphaxx*kx^2+alphayy*ky^2;
h12=betaxx* kx^2+betayy *ky^2+1i*betaxy*kx*ky; 
h13=gamaxx* kx^2+gamayy *ky^2+1i*gamaxy*kx*ky;
h33=muxx  * kx^2+muyy   *ky^2;
h22= h11';
h23=-h13';
h21=h12';
h31=h13';
h32=h23';
%
H0_AlN = [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
%kane's band edge
H0_AlNE= [Delta1+Delta2 0 0  0      0           0; ...
          0 Delta1-Delta2 0  0      0   s2*Delta3; ...
          0       0       0  0  s2*Delta3       0; ...
          0       0       0  Delta1+Delta2  0   0; ...
          0   0   s2*Delta3  0  Delta1-Delta2   0; ...
          0   s2*Delta3   0  0        0         0];
%strain induced
H0_AlND= [l1*exx+m1*eyy+m2*ezz  n1* exy  n2*exz  0  0  0; ...
          n1*exy m1*exx+l1*eyy+m2*ezz    n2*eyz  0  0  0; ...
          n2*exz  n2*eyz  m3*exx+m3*eyy+l2*ezz   0  0  0; ...
          0  0  0  l1*exx+m1*eyy+m2*ezz  n1* exy  n2*exz; ...
          0  0  0  n1*exy m1*exx+l1*eyy+m2*ezz    n2*eyz; ...
          0  0  0  n2*exz  n2*eyz  m3*exx+m3*eyy+l2*ezz ];
%include the A7 parameter
H0_AlNA7= [0  0     1i*A7*(kx-1i*ky)    0   0   0;...
           0  0    -1i*A7*(kx+1i*ky)    0   0   0;...
          -1i*A7*(kx+1i*ky) 1i*A7*(kx-1i*ky) 0 0 0 0;...
          0 0 0    0  0   -1i*A7*(kx+1i*ky);...
          0 0 0    0  0    1i*A7*(kx-1i*ky);...
          0 0 0    1i*A7*(kx-1i*ky)  -1i*A7*(kx+1i*ky) 0]; 
H0_AlN =H0_AlN+H0_AlND+H0_AlNE+H0_AlNA7;
Hcc=[kx^2*ax+ky^2*ay+Eg 0;0 kx^2*ax+ky^2*ay+Eg];
Hcv=[-P2*(kx+1i*ky)/s2 P2*(kx-1i*ky)/s2 0 0 0 0;...
      0 0 0 P2*(kx-1i*ky)/s2 -P2*(kx+1i*ky)/s2 0];
Hvc=Hcv';
H0_AlN=[Hcc Hcv; Hvc H0_AlN];
% assemble only when computing the polarization
if (nargin>3)
H0_AlNte = [zeros(2) Hcv; Hvc H0_AlNA7]/sqrt(kx^2+ky^2);
Hcv=[0 0 P1 0 0 0;0 0 0 0 0 P1]; Hvc=Hcv'; 
H0_AlNtm = [zeros(2) Hcv; Hvc zeros(6)];
H0 =  (1-xmol)*H0_GaNte  + xmol*H0_AlNte;
H2 =  (1-xmol)*H0_GaNtm  + xmol*H0_AlNtm;
H1L=0;H1R=0; return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=alphazx*kx;   
h12=1i*betazy*ky;
h13=0;            
h33=muzx*kx;      
h22= h11' ;       
h23=0;           
h21=-1i*ky*betayz;          
h31= gamaxz'*kx-1i*ky*gamayz';
h32=-gamaxz'*kx-1i*ky*gamayz';
%                                                       kz*xxx
H1_AlNL= [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
%                                                     
h11=alphaxz*kx;  
h12=1i*betayz*ky; 
h13=gamaxz*kx+1i*gamayz*ky;
h33=muxz*kx;      
h22=h11' ;       
h23=-gamaxz*kx+1i*gamayz*ky;
h21=-1i*ky*betazy;
h31=0;
h32=0;      
%                                                       xxx*kz
H1_AlNR= [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
%
Hcc=[0 0;0 0];
Hcv=[0 0 P1 0 0 0;...
     0 0 0 0 0 P1]*0.5;    
H1_AlNL=[Hcc Hcv; Hcv' H1_AlNL];
H1_AlNR=[Hcc Hcv; Hcv' H1_AlNR];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=alphazz;
h12=betazz; 
h13=gamazz; 
h33=muzz;
h22=h11';
h23=-h13';  
h21=h12';   
h31=h13';   
h32=h23';   
%
H2_AlN = [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
Hcc=[az 0;0 az];
Hcv=[0 0 0 0 0 0;...
     0 0 0 0 0 0];
Hvc=Hcv';
H2_AlN=[Hcc Hcv; Hvc H2_AlN];
%==========================================================================
% AlGaN paramerts are linearly interploated between the binaries 
%==========================================================================
H0 =  (1-xmol)*H0_GaN  + xmol*H0_AlN;
H1L = (1-xmol)*H1_GaNL + xmol*H1_AlNL;
H1R = (1-xmol)*H1_GaNR + xmol*H1_AlNR;
H2 =  (1-xmol)*H2_GaN  + xmol*H2_AlN;
%
H0  = H0/Q;
H1L = H1L/Q;
H1R = H1R/Q;
H2  = H2/Q;
%==========================================================================
end