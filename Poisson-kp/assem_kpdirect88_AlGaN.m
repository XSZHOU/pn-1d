function [H0,H1L,H1R,H2] = assem_kpdirect88_AlGaN(kx,ky,xmol)
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
Delta1 = 0.0427*Q;
%
Delta2 = 0.000006/3*Q;
Delta3 = 0.000006/3*Q;
%
mz=1.1703*M0;   
mt=1.2815*M0;
%
Eg1=3.52076*Q;
Eg=Eg1+Delta1+Delta2;                                           %19/02/2013
%
P11=60.77*Q*Q*1e-20;
P22=47.22*Q*Q*1e-20;
P1=sqrt(P11);
P2=sqrt(P22);
ax=H2M0*(1/mt)*M0;
ay=ax;
az=H2M0*(1/mz)*M0;
%============================================================1996Chuang_PRB
A1 = -2.0248;
A2 = -0.8279;
A3 =  1.0861;
A4 = -0.53;
A5 = -0.4544;
A6 =  0.175;
A7 = 0.408*1e-10*Q;
%
D1 = -3.6*Q;
D2 =  1.7*Q; 
D3 =  5.2*Q;
D4 = -2.7*Q;
D5 = -2.8*Q; 
D6 = -4.3*Q;
%
a0 = 3.189;
a  = a0;
%
c0 = 5.185;
c  = c0;
%
C13 = 106;
C33 = 398; 
%
exx = (a0-a)/a;
%
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
gamaxz=-A6;                            %Veprek's form: A5, A6 have more 
gamazx=0;                              %24/11/2012
gamayz= A6;                            %02/04/2013
gamazy=0;                              %02/04/2013
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
h12=betaxx* kx^2+betayy *ky^2+1i*betaxy*kx*ky; %%err in mireles1999
h13=gamaxx* kx^2+gamayy *ky^2+1i*gamaxy*kx*ky; % 0 @ [0001]
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
          0   s2*Delta2   0  0        0         0];
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
          0 0 0    1i*A7*(kx-1i*ky)  -1i*A7*(kx+1i*ky) 0]; %Dugdale in |u>
      %
H0_GaN =H0_GaN+H0_GaND+H0_GaNE+H0_GaNA7;
Hcc=[kx^2*ax+ky^2*ay+Eg 0;0 kx^2*ax+ky^2*ay+Eg];
Hcv=[-P2*(kx+1i*ky)/s2 P2*(kx-1i*ky)/s2 0 0 0 0;...
      0 0 0 P2*(kx-1i*ky)/s2 -P2*(kx+1i*ky)/s2 0];
Hvc=Hcv';
H0_GaN=[Hcc Hcv; Hvc H0_GaN];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=alphazx*kx;   % 0 @ [0001]
h12=1i*betazy*ky; % 0 @ [0001]
h13=0;            % 
h33=muzx*kx;      % 0 @ [0001]
h22= h11' ;       % 0 @ [0001]
h23=0;            %
h21=-1i*ky*betayz;          %  0 @ [0001]
h31= gamaxz*kx-1i*ky*gamayz;%  2000Mireles(10)--veprek (3.90 ita')
h32=-gamaxz*kx-1i*ky*gamayz;%  2000Mireles(10)--veprek (3.90 -ita)
%                                                       kz*xxx
H1_GaNL= [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
%                                                     
h11=alphaxz*kx;   % 0 @ [0001]
h12=1i*betayz*ky; % 0 @ [0001]
h13=gamaxz*kx+1i*gamayz*ky;
h33=muxz*kx;      % 0 @ [0001]
h22=h11' ;        % 0 @ [0001]
h23=-gamaxz*kx+1i*gamayz*ky;
h21=-1i*ky*betazy;% 0 @ [0001]
h31=0;
h32=0;      
%                                                      xxx*kz
H1_GaNR= [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
%
Hcc=[0 0;0 0];
Hcv=[0 0 P1 0 0 0;...
     0 0 0 0 0 P1]/2;    
Hvc=Hcv';
H1_GaNL=[Hcc Hcv; Hvc H1_GaNL];
H1_GaNR=[Hcc Hcv; Hvc H1_GaNR];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=alphazz;
h12=betazz; %  0 @ [0001]
h13=gamazz; %  0 @ [0001]
h33=muzz;
h22=h11';
h23=-h13';  %  0 @ [0001]
h21=h12';   %  0 @ [0001]
h31=h13';   %  0 @ [0001]
h32=h23';   %  0 @ [0001]
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
%
%==========================================================================
% AlN
%============================================================8by8parameters
Delta1 = -0.20167*Q;
%
Delta2 = 0.00006/3*Q;   % 0.036 in standard reference
Delta3 = 0.00006/3*Q;
%
mz=1.1518*M0; 
mt=1.1658*M0;
%
Eg1=5.8504*Q - Delta1;  % Eg1= |S> -|hh> in kp sense
Eg=Eg1+Delta1+Delta2;    
%
P11=54.3971*Q*Q*1e-20;
P22=49.5617*Q*Q*1e-20;
%
P1=sqrt(P11);
P2=sqrt(P22);
ax=H2M0*(1/mt)*M0;
ay=ax;
az=H2M0*(1/mz)*M0;
%============================================================1996Chuang_PRB
%
A1 = -1.54;
A2 = -0.61677;
A3 =  0.966;
A4 = -0.255;
A5 = -0.3124;
A6 = -6.45e-2;
A7 = 0.138*1e-10*Q;
%
D1 = -3.6*Q;
D2 =  1.7*Q; 
D3 =  5.2*Q;
D4 = -2.7*Q;
D5 = -2.8*Q; 
D6 = -4.3*Q;
%
a0 = 3.189;
a  = a0;
%
c0 = 5.185;
c  = c0;
%
C13 = 106;
C33 = 398; 
%
exx = (a0-a)/a;
%
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
gamaxz=-A6;                            %Veprek's form: A5, A6 have more 
gamazx=0;                              %24/11/2012
gamayz= A6;                            %02/04/2013
gamazy=0;                              %02/04/2013
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
h12=betaxx* kx^2+betayy *ky^2+1i*betaxy*kx*ky; %
h13=gamaxx* kx^2+gamayy *ky^2+1i*gamaxy*kx*ky; %
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
          0   s2*Delta2   0  0        0         0];
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
          0 0 0    1i*A7*(kx-1i*ky)  -1i*A7*(kx+1i*ky) 0]; %Dugdale in |u>
      %
H0_AlN =H0_AlN+H0_AlND+H0_AlNE+H0_AlNA7;
Hcc=[kx^2*ax+ky^2*ay+Eg 0;0 kx^2*ax+ky^2*ay+Eg];
Hcv=[-P2*(kx+1i*ky)/s2 P2*(kx-1i*ky)/s2 0 0 0 0;...
      0 0 0 P2*(kx-1i*ky)/s2 -P2*(kx+1i*ky)/s2 0];
Hvc=Hcv';
H0_AlN=[Hcc Hcv; Hvc H0_AlN];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=alphazx*kx;   % 0 @ [0001]
h12=1i*betazy*ky; % 0 @ [0001]
h13=0;           
h33=muzx*kx;      % 0 @ [0001]
h22= h11' ;       % 0 @ [0001]
h23=0;           
h21=-1i*ky*betayz;          %  0 @ [0001]
h31= gamaxz*kx-1i*ky*gamayz;
h32=-gamaxz*kx-1i*ky*gamayz;
%                                                       kz*xxx
H1_AlNL= [h11     h12     h13      0    0    0  ; ...
          h21     h22     h23      0    0    0  ; ...
          h31     h32     h33      0    0    0  ; ...
          0       0       0       h11  h21  h23 ; ...
          0       0       0       h12  h22  h13 ; ...
          0       0       0       h32  h31  h33]*H2M0;
%                                                     
h11=alphaxz*kx;   % 0 @ [0001]
h12=1i*betayz*ky; % 0 @ [0001]
h13=gamaxz*kx+1i*gamayz*ky;
h33=muxz*kx;      % 0 @ [0001]
h22=h11' ;        % 0 @ [0001]
h23=-gamaxz*kx+1i*gamayz*ky;
h21=-1i*ky*betazy;% 0 @ [0001]
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
     0 0 0 0 0 P1]/2;    
Hvc=Hcv';
H1_AlNL=[Hcc Hcv; Hvc H1_AlNL];
H1_AlNR=[Hcc Hcv; Hvc H1_AlNR];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=alphazz;
h12=betazz; %  0 @ [0001]
h13=gamazz; %  0 @ [0001]
h33=muzz;
h22=h11';
h23=-h13';  %  0 @ [0001]
h21=h12';   %  0 @ [0001]
h31=h13';   %  0 @ [0001]
h32=h23';   %  0 @ [0001]
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