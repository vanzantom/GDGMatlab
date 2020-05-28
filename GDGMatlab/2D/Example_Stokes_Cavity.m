%-------------------------------------------------------------------------
% example_Lid_cavity is an example of script using the GDGMatlab package to
% solve the Lid cavity problem

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
clear all;

% The structure geo contains a function describing the geometry of the problem and the meshsize. 
geo.fun='@square';% decomposition of the square (0,1)^2 
geo.regular=1; % flag to have a regular triangular mesh.( geo.regular needs to set geo.fun='square').
geo.h=1/8;
% The structure data contains information about the data of the problem,
% e.g. force term, diffusion coefficient, boundary conditions
data.f1=@(x,y) 0; %force term f=(f1,f2);
data.f2=@(x,y) 0; %force term 
data.mu=@(x,y) 1; % viscosity
data.Dirichlet_u=@(x,y) (y==1)*1; % boundary condition_on top boundary I impose u=1
data.Dirichlet_v=@(x,y) 0; % boundary condition_on top boundary I impose u=1
data.pset=@(x,y) 1;
data.labelu=[-1,-1,-1,-1];
data.labelv=[-1,-1,-1,-1];


% The structure para contains parameters of the problem such as the order
% of Gauss quadrature.
para.order_Gauss=2; %Order Gauss quadrature.
para.mu=1;% viscosity
para.basis_type_ve=202; 
para.basis_type_pre=201;
para.index_fix_pressure=2;
% Solver
[u,v,p,T,Eb,Pbve,Tbve,Pbpre,Tbpre]=Stokes_solver_2D(geo,data,para);

figure(1)
pdeplot(Pbve,[],Tbve,'xydata',u,'xystyle','flat','zdata',u,'zstyle','continuous','colorbar','on');
title('Velocity field u')
xlabel('x');ylabel('y');
figure(2)
pdeplot(Pbve,[],Tbve,'xydata',v,'xystyle','flat','zdata',v,'zstyle','continuous','colorbar','on');
title('Velocity field v')
xlabel('x');ylabel('y');
figure(3)
pdeplot(Pbpre,[],[Tbpre;T(4,:)],'xydata',p,'xystyle','flat','zdata',p,'zstyle','continuous','colorbar','on');
title('Pressure field p')
xlabel('x');ylabel('y');
figure(4)
pdeplot(Pbve,[],Tbve,'Flowdata',[u,v],'mesh','on');
title('Velocity field [u,v]')
xlabel('x');ylabel('y');
