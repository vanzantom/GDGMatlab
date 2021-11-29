%-------------------------------------------------------------------------
% example_Stokes is an example of script using the GDGMatlab package to
% solve a Stokes problem.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
addpath('lib');
addpath('solvers');
clear all;
% The structure geo contains a function describing the geometry of the problem and the meshsize. 
geo.fun='@square';% decomposition of the square (0,1)^2 
geo.regular=1; % flag to have a regular triangular mesh.( geo.regular needs to set geo.fun='square').
% The structure para contains parameters of the problem such as the order
% of Gauss quadrature.
para.order_Gauss=4; %Order Gauss quadrature.
para.mu=1;% viscosity
para.basis_type_ve=2011; 
para.basis_type_pre=201;
para.index_fix_pressure=1;
para.Robin=1;

% The structure data contains information about the data of the problem,
% e.g. force term, diffusion coefficient, boundary conditions
uex=@(x,y) 20*x*y^3;
vex=@(x,y) 5*x^4-5*y^4;
pex=@(x,y) 60*x^2*y-20*y^3;
data.mu=@(x,y,El) 1; % viscosity
data.f1=@(x,y) 0;
data.f2=@(x,y) 0; %force term 
data.pset=@(x,y) pex(x,y);%support function to fix pressure
data.Dirichlet_u=@(x,y) uex(x,y); 
data.Dirichlet_v=@(x,y) vex(x,y); 
data.Normal_stress=@(x,y) 2*para.mu*(-20*y^3) -pex(x,y);
data.Robin_stress=@(x,y) (2*para.mu*(-2*x*y^2) -pex(x,y))+para.Robin*vex(x,y); %(n'*T(u,p)*n + p n'*(u,v))
data.labelu=[-1,-1,-1,-1];% -1 Dirichlet, -2 Neumann,-3 Robin. [Top,right,bottom,left]
data.labelv=[-1,-1,-1,-1];

plt=1;
hh=[1/4;1/8;1/16;1/32;1/64];
for i=1:length(hh)
    geo.h=hh(i);
    [u,v,p,T,Eb,Pbve,Tbve,Pbpre,Tbpre]=Stokes_solver_2D(geo,data,para);
    %Evaluate exact solution on the mesh
    for k=1:size(Pbve,2)
        x=Pbve(1,k); y=Pbve(2,k);
        uexvec(k)=uex(x,y);
        vexvec(k)=vex(x,y);
    end
    for k=1:size(Pbpre,2)
        x=Pbpre(1,k); y=Pbpre(2,k);
        pexvec(k)=pex(x,y);
    end
    if(plt==1)
        figure(1)
        pdeplot(Pbve,[],Tbve,'xydata',u,'xystyle','flat','zdata',u,'zstyle','continuous','colorbar','on');
        title('Velocity field u');xlabel('x');ylabel('y');
        figure(2)
        pdeplot(Pbve,[],Tbve,'xydata',v,'xystyle','flat','zdata',v,'zstyle','continuous','colorbar','on');
        title('Velocity field v');xlabel('x');ylabel('y');
        figure(3)
        pdeplot(Pbpre,[],[Tbpre;T(4,:)],'xydata',p,'xystyle','flat','zdata',p,'zstyle','continuous','colorbar','on');
        title('Pressure field p');xlabel('x');ylabel('y');
        figure(4)
        pdeplot(Pbve,[],Tbve,'Flowdata',[u,v],'mesh','on');
        title('Velocity field [u,v]');xlabel('x');ylabel('y');
    end
    
    [error_uL2(i),error_uH1(i)]=compute_errors_Stokes(u,uexvec,Pbve,Tbve,Tbve((1:3),:),Eb,para.order_Gauss,para.basis_type_ve); 
    [error_vL2(i),error_vH1(i)]=compute_errors_Stokes(v,vexvec,Pbve,Tbve,Tbve((1:3),:),Eb,para.order_Gauss,para.basis_type_ve);
    [error_pL2(i),~]=compute_errors_Stokes(p,pexvec,Pbpre,Tbpre,Tbpre((1:3),:),Eb,para.order_Gauss,para.basis_type_pre);
   
end
loglog(hh,error_uH1+error_vH1+error_pL2,'b-x');
    hold on
loglog(hh,hh.^2,'r')
loglog(hh,hh,'k')
grid on
