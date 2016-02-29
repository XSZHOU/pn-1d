function mesh=load_params(mesh,mate_sys)
%
Q=  1.6021917e-19;
KB= 1.3807*1e-23;                  
Eps=8.854*1e-12;
%
Eg_AlN=6.00; a_AlN=3.112;            %%%%%%%%%%%% 
Eg_GaN=3.437;a_GaN=3.189;            a_sub=a_GaN;
Eg_InN=0.69; a_InN=3.545;            %%%%%%%%%%%%
bow_AlGaN= 0.0; bow_InGaN= 0.0;
VBO_AlGaN= 0.3; VBO_InGaN= 0.3;
AlN_eps=Eps*10.1;
GaN_eps=Eps*10.4;                                                          %Piprek,pp.215
InN_eps=Eps*15.3;                                                          
%
discont=diff(mesh.ii);                                                     %discontinuity capturing
ind = find(discont==1);                                                    %find corresponding indices
ind = ind(end);                                                              %
%==========================================================================
switch (mate_sys)
case{'AlGaN'}
bow=bow_AlGaN;
ratio =VBO_AlGaN;
mesh.gap =Eg_AlN*mesh.xmol+Eg_GaN*(1-mesh.xmol)-bow*mesh.xmol.*(1-mesh.xmol);
mesh.vbo =- ratio*(mesh.gap(ind) - mesh.gap(ind+1));
mesh.evb  = mesh.vbo*ones(1,length(mesh.xc));mesh.evb(mesh.ii)=0;
%-------------------------------------------------------------$Polarization
mesh.Psp =mesh.xmol*(-0.0898)+(1-mesh.xmol)*(-0.0339) ...                  %Piprek, pp.53,59, note ! opposite sign eq.2.9~eq.3.14
          +0.019* mesh.xmol.*(1-mesh.xmol);
mesh.exx =(a_sub-mesh.xmol*a_AlN-a_GaN*(1-mesh.xmol)) ...                  %Piprek, pp.53,59
         ./(mesh.xmol*a_AlN+a_GaN*(1-mesh.xmol));
mesh.Ppz =(-1.808*mesh.exx).*mesh.xmol+ ...                                %Piprek, pp.60, [$approximated$] eq.3.16
          (-0.918*mesh.exx+9.541*mesh.exx.*mesh.exx).*(1-mesh.xmol);     
mesh.Pol = mesh.Ppz + mesh.Psp;
%-------------------------------------------------------------$Polarization
mesh.eps = GaN_eps*(1-mesh.xmol)+AlN_eps*mesh.xmol;
mesh.piechar = -[0 diff(mesh.Pol) 0]';                                     %Piprek, pp.55, note ! -div(P),also Veprek pp.28
mesh.target1 =  mesh.gap(ind+1);                                           %initial guess for target, flat band
mesh.target2 =  3e-1;                                                      %initial guess for target
mesh.ecb=mesh.evb;
%==========================================================================
case{'InGaN'}
bow =  bow_InGaN;                                                          %note! 2010Moses_APL -> 1.4
ratio =VBO_InGaN;
mesh.gap =Eg_InN*mesh.xmol+Eg_GaN*(1-mesh.xmol)-bow*mesh.xmol.*(1-mesh.xmol);
mesh.vbo  = -ratio*(mesh.gap(ind) - mesh.gap(ind+1));
mesh.evb  = mesh.vbo*ones(1,length(mesh.xc));mesh.evb(mesh.ii)=0;
%-------------------------------------------------------------$Polarization
mesh.Psp =mesh.xmol*(-0.0413)+(1-mesh.xmol)*(-0.0339) ...                  %Piprek, pp.53,59, note opposite sign eq.2.9~eq.3.14
          +0.03* mesh.xmol.*(1-mesh.xmol);                                 %[$approximated$], table 3.4
mesh.exx =(a_sub-mesh.xmol*a_InN-a_GaN*(1-mesh.xmol)) ...                  %Piprek, pp.53,59
         ./(mesh.xmol*a_InN+a_GaN*(1-mesh.xmol));                            
mesh.Ppz =(-1.373*mesh.exx+7.559*mesh.exx.*mesh.exx).*mesh.xmol+ ...       %Piprek, pp.60
          (-0.918*mesh.exx+9.541*mesh.exx.*mesh.exx).*(1-mesh.xmol);
mesh.Pol = mesh.Ppz + mesh.Psp;
%-------------------------------------------------------------$Polarization
mesh.eps = GaN_eps*(1-mesh.xmol)+InN_eps*mesh.xmol;
mesh.piechar = -[0 diff(mesh.Pol) 0]';                                     %Piprek, pp.55, note ! -div(P),also Veprek pp.28
mesh.target1 = mesh.gap(ind+1);                                            %initial guess for target, flat band
mesh.target2 = 3e-1 ;                                                      %initial guess for target
mesh.ecb=mesh.evb;
%
otherwise, error('specify correct material system!s');
end