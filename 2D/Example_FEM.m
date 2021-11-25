%-------------------------------------------------------------------------
% example_FEM is an example of script using Poisson_solver_2D function to solve
% a Poisson problem on the geometry specified and P1 or P2 FE
% space. It is possible to specify Dirichlet/Neumann and Robin boundary
% condition
% Further, it checks the order of convergence.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
clear all;
addpath('lib');
addpath('solvers');

% The structure geo contains a function describing the geometry of the problem and the meshsize. 
geo.fun='@square';% decomposition of the square (0,1)^2 into 2 subdomains
geo.regular=1; % flag to have a regular triangular mesh.( geo.regular needs to set geo.fun='square').
% The structure para contains parameters of the problem such as the order
% of Gauss quadrature.
para.order_Gauss=2; %Order Gauss quadrature.
para.basis_type=202; %P1 FE space 201, P2 FE space 202
para.Robin=1;
% The structure data contains information about the data of the problem,
% e.g. force term, diffusion coefficient, boundary conditions
uex=@(x,y) sin(pi/2*x).*sin(pi/2*y); % Exact solution
data.f=@(x,y) (1/2)*pi^2*sin(pi/2*x)*sin(pi/2*y); %force term 
data.c=@(x,y,El) 1; % diffusion term
data.Dirichlet_fun=@(x,y) uex(x,y); % boundary condition
data.Neumann_fun=@(x,y) (pi/2)*cos(pi/2*x)*sin(pi/2*y);
data.Robin_fun=@(x,y) (pi/2)*cos(pi/2*x)*sin(pi/2*y)+para.Robin*sin(pi/2*x)*sin(pi/2*y);
data.label=[-1,-3,-1,-1];% Edge -1 Dirichlet,-2 Neumann, -3 Robin. [Top,right,bottom,left]


plt=0; % Set equal to 0 to disable plot, equal to 1 to plot solution
%=================================
%   Check Convergence
%=================================
hh=[1/2;1/8;1/16;1/32]; % Vector of mesh sizes
for i=1:length(hh)
geo.h=hh(i);
[P,E,T,Pb,Tb,Eb,u]=Poisson_solver_2D_FEM(geo,data,para);%Solve Poisson Problem
for k=1:size(Pb,2)%evaluate exact solution on the mesh nodes.
	x=Pb(1,k); y=Pb(2,k);
	uexvec(k)=uex(x,y);
end
[error_L2(i),error_H1(i)]=compute_errors(u,uexvec,P,T,Tb((1:3),:),E,para,'FEM');%Compute L2 and H1 errors.
if plt==1
    figure(1)
    visualize_solution_mesh(T,Pb,Tb,u,1,1)
%   figure(2)
%   pdemesh(P,E,T,u, 'mesh', 'on');% Plot the mesh
%   pdemesh(P,E,T,uexvec,); % Plot the exact solution.
    pause
end
end
loglog(hh,hh,'m')
hold on
loglog(hh,hh.^2,'r')
loglog(hh,hh.^3,'c')
loglog(hh,error_L2,'b*-')
loglog(hh,error_H1,'k*-')
grid on
xlabel('h')
ylabel('errors')
legend({'$$h$$','$$h^{2}$$','$$h^{3}$$','$$\|u_h-u\|_{L^2}$$','$$\|u_h-u\|_{H_1}$$',},'interpreter','latex')
set(gca,'FontSize',12);