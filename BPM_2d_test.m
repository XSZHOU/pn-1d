clear all;clc;
lambda=0.5;
W=100;
L=100;
%2d YZ
DY=lambda/2;
DZ=lambda/2;
n0=1.6+0.01j;

IndexFun=@(z,y) n0+ 1*((abs(y)<1));
FieldFun=@(y) exp(-(y/20).^2);

k=2*pi/lambda;
PointsY=round(W/DY)+1;

if rem(PointsY,2)==0
    PointsY=PointsY+1;
end

Y=linspace(-W/2,W/2,PointsY).';
DY=Y(2)-Y(1);

field=FieldFun(Y);
%field=double((abs(Y)<100));
%plot(Y,field)

PointsZ=round(L/DZ)+1;
Z=linspace(0,L,PointsZ);
DZ=Z(2)-Z(1);
indexfield=IndexFun(Z,Y);

Aconst=sparse(PointsY-2,PointsY-2);

for m=1:PointsY-2
    Aconst(m,m)=2+DY^2*k*4j*n0/DZ;
end

for m=2:PointsY-2
    Aconst(m,m-1)=-1;
end

for m=1:PointsY-3
    Aconst(m,m+1)=-1;
end

Bconst=sparse(PointsY-2,PointsY);

for m=1:PointsY-2
    Bconst(m,m+1)=-2+DY^2*k*4j*n0/DZ;
    Bconst(m,m)=1;
    Bconst(m,m+2)=1;
end

FieldMatrix=zeros(length(Y), length(Z));

for kz=0:PointsZ-1;                      %main longitudinal loop 
    indexlocal=IndexFun(Z(kz+1)+DZ/2,Y);
    %if((Z(kz+1)>0.15*L)&&(Z(kz+1)<0.25*L))
    %    indexlocal=indexlocal+1;
    %end
    AverageIndex=indexlocal;
    dnSquared=AverageIndex.^2-n0^2;
    C=sparse(PointsY-2,PointsY-2);
    for m=1:PointsY-2
        C(m,m)=-(DY*k)^2*dnSquared(m+1);
    end
    D=sparse(PointsY-2,PointsY);
    for m=1:PointsY-2
        D(m,m+1)=(DY*k)^2*dnSquared(m+1);
    end
    A=Aconst+C;
    B=Bconst+D;
    if field(3)~=0
        ExpL=field(2)/field(3);
        if imag(log(ExpL))>0
           ExpL=conj(ExpL);
        end
    else
        ExpL=1;
    end
    field(1)=field(2)*ExpL;
    A(1,1)=A(1,1)-ExpL;
    if field(PointsY-2)~=0 %transparent boundary contitions 
       ExpR=field(PointsY-1)/field(PointsY-2); 
       if imag(log(ExpR))>0 %handlenonphysicalsolutions
       ExpR=conj(ExpR);
       end
    else
       ExpR=1; 
    end
    field(PointsY)=field(PointsY-1)*ExpR;
    A(PointsY-2,PointsY-2)=A(PointsY-2,PointsY-2)-ExpR;
    field(2:PointsY-1)=A\(B*field);
    field(1)= field(2)*ExpL;
    field(PointsY)= field(PointsY-1)*ExpR;
    FieldMatrix(:,kz+1)=abs(field).^2; 
end
figure;
pcolor(Z,Y,FieldMatrix); xlabel('Z [um]'); ylabel('Y [um]'); 
shading interp;
