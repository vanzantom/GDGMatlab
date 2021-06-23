%-------------------------------------------------------------------------
% Example_Stokes_Darcy is an example of script using the GDGMatlab package to
% solve a Stokes-Darcy system.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
clear all;
addpath('lib')
addpath('solvers')
% The structure geo contains a function describing the geometry of the problem and the meshsize. 
geo.regular=0; % flag to have a regular triangular mesh.( geo.regular needs to set geo.fun='square').
geo.funDarcy='@square_translated'; %geometry of Darcy problem
geo.funStokes='@square'; %geometry of Stokes problem.
geo.label=-3; % label for the interface
%===============
%Parameters
%===============
paraDarcy.order_Gauss=2; %Order Gauss quadrature.
paraDarcy.basis_type=202; %P2 FE space 202
paraDarcy.Robin=0;
paraStokes.order_Gauss=2; %Order Gauss quadrature.
paraStokes.mu=1;
paraStokes.basis_type_ve=202; %P2 FE space 202
paraStokes.basis_type_pre=201; %P2 FE space 202
paraStokes.index_fix_pressure=0;
paraStokes.Robin=0;

%=================================
%   Darcy
%=================================
qex=@(x,y) -(2-pi*sin(pi*x))*(y+1);
dataDarcy.g1=@(x,y) pi^2*cos(pi*x)*(y + 1); %force term 
dataDarcy.g2=@(x,y)0 ;
dataDarcy.c=@(x,y,El) 1; % diffusion term
dataDarcy.Dirichlet_fun=@(x,y) qex(x,y); % boundary condition
dataDarcy.Neumann_fun=@(x,y)  0; % -(2-pi*sin(pi*x));%;
dataDarcy.label=[-3,-1,-1,-1];% Robin boundary condition on the top edge.
%=================================
%   Stokes
%=================================
uex=@(x,y) x^2*y^2+exp(-y);
vex=@(x,y) -2/3*x*y^3+2-pi*sin(pi*x);
pex=@(x,y) -(2-pi*sin(pi*x))*cos(2*pi*y);
dataStokes.f1=@(x,y) -2*paraStokes.mu*x^2-2*paraStokes.mu*y^2-paraStokes.mu*exp(-y)+pi^2*cos(pi*x)*cos(2*pi*y); %force term 
dataStokes.f2=@(x,y) 4*paraStokes.mu*x*y-paraStokes.mu*pi^3*sin(pi*x)+2*pi*(2-pi*sin(pi*x))*sin(2*pi*y); %force term 
dataStokes.mu=@(x,y,El) paraStokes.mu; % diffusion term
dataStokes.Dirichlet_u=@(x,y) uex(x,y); % boundary condition
dataStokes.Dirichlet_v=@(x,y) vex(x,y); % boundary condition
dataStokes.Normal_stress=@(x,y) 0;
dataStokes.labelu=[-1,-1,-1,-1];% -1 Dirichlet, -2 Neumann,-3 Robin.
dataStokes.labelv=[-1,-1,-3,-1];% 

plt=0;
hh=[1/4;1/8;1/16;1/32];
for i=1:length(hh);
    geo.h=hh(i);
    [result,Atot,C,matrixS,FES,matrixD,FED]=Stokes_Darcy_solver_2D(geo,dataStokes,paraStokes,dataDarcy,paraDarcy);
    for k=1:size(FES.Pb,2)
            x=FES.Pb(1,k); y=FES.Pb(2,k);
            uexvec(k)=uex(x,y);
            vexvec(k)=vex(x,y);
        end
        for k=1:size(FES.P,2)
            x=FES.P(1,k); y=FES.P(2,k);
            pexvec(k)=pex(x,y);
       end
       for k=1:size(FED.Pb,2)%evaluate exact solution on the mesh nodes.
            x=FED.Pb(1,k); y=FED.Pb(2,k);
            qexvec(k)=qex(x,y);
       end
    u=result(1:size(FES.Pb,2));
    v=result(size(FES.Pb,2)+1:2*size(FES.Pb,2));
    p=result(2*size(FES.Pb,2)+1:2*size(FES.Pb,2)+size(FES.P,2));
    q=result(2*size(FES.Pb,2)+size(FES.P,2)+1:end);
       
    if(plt==1)
        
        figure(1)
        subplot(1,2,1)
        pdeplot(FES.Pb,[],FES.Tb,'xydata',u,'xystyle','flat','zdata',u,'zstyle','continuous','colorbar','on');
        title('Velocity field u');xlabel('x');ylabel('y');
        subplot(1,2,2)
        pdeplot(FES.Pb,[],FES.Tb,'xydata',uexvec,'xystyle','flat','zdata',uexvec,'zstyle','continuous','colorbar','on');
        title('Velocity field u');xlabel('x');ylabel('y');
        figure(2)
        subplot(1,2,1)
        pdeplot(FES.Pb,[],FES.Tb,'xydata',v,'xystyle','flat','zdata',v,'zstyle','continuous','colorbar','on');
        title('Velocity field v');xlabel('x');ylabel('y');
        subplot(1,2,2)
        pdeplot(FES.Pb,[],FES.Tb,'xydata',vexvec,'xystyle','flat','zdata',vexvec,'zstyle','continuous','colorbar','on');
        title('Velocity field v');xlabel('x');ylabel('y');
        figure(3)
        subplot(1,2,1)
        pdeplot(FES.P,[],[FES.T],'xydata',p,'xystyle','flat','zdata',p,'zstyle','continuous','colorbar','on');
        title('Pressure field p');xlabel('x');ylabel('y');
        subplot(1,2,2)
        pdeplot(FES.P,[],[FES.T],'xydata',pexvec,'xystyle','flat','zdata',pexvec,'zstyle','continuous','colorbar','on');
        title('Pressure field p');xlabel('x');ylabel('y');
%         figure(4)
%         subplot(1,2,1)
%         pdeplot(Pbve,[],Tbve,'Flowdata',[u,v],'mesh','on');
%         title('Velocity field [u,v]');xlabel('x');ylabel('y');
%         subplot(1,2,2)
%          pdeplot(Pbve,[],Tbve,'Flowdata',[uexvec,vexvec],'mesh','on');
%         title('Velocity field [u,v]');xlabel('x');ylabel('y');
%       
       figure(5)
        subplot(1,2,1)
        pdeplot(FED.Pb,[],FED.Tb,'xydata',q,'xystyle','flat','zdata',q,'zstyle','continuous','colorbar','on');
        title('Pressure Darcy');xlabel('x');ylabel('y');
        subplot(1,2,2)
        pdeplot(FED.Pb,[],FED.Tb,'xydata',qexvec,'xystyle','flat','zdata',qexvec,'zstyle','continuous','colorbar','on');
        title('Pressure Darcy');xlabel('x');ylabel('y');
       
        pause
    end
     [error_uL2(i),error_uH1(i)]=compute_errors_Stokes(u,uexvec,FES.Pb,FES.Tb,FES.Tb((1:3),:),FES.Eb,paraStokes.order_Gauss,paraStokes.basis_type_ve); 
    [error_vL2(i),error_vH1(i)]=compute_errors_Stokes(v,vexvec,FES.Pb,FES.Tb,FES.Tb((1:3),:),FES.Eb,paraStokes.order_Gauss,paraStokes.basis_type_ve);
    [error_pL2(i),~]=compute_errors_Stokes(p,pexvec,FES.P,FES.T,FES.T((1:3),:),FES.E,paraStokes.order_Gauss,paraStokes.basis_type_pre);
     [error_qL2(i),~]=compute_errors_Stokes(q,qexvec,FED.Pb,FED.Tb,FED.Tb((1:3),:),FED.E,paraDarcy.order_Gauss,paraDarcy.basis_type);
   error(i)=error_uH1(i)+error_vH1(i)+error_pL2(i)+error_qL2(i);

end
loglog(hh,error,'k-o','Linewidth',1.3);
hold on
loglog(hh,hh.^2,'k--','Linewidth',1.3)
grid on
xlabel('Mesh size')
legend({'$Err(\mathbf{u}^h_s,\mathbf{p}_s^h,\mathbf{p}_d^h)$','$h^2$'},'interpreter','latex')
set(gca,'FontSize',16)