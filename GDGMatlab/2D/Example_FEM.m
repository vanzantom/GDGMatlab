%-------------------------------------------------------------------------
% example_FEM is an example of script using Poisson_solver_2D function to solve
% a Poisson problem on the geometry specified with Dirichlet BC on the boundary and P1 or P2 FE
% space.
% Further, it checks the order of convergence.

% Attention: Currently, the function generate_boundarynodes works only if
% the geometry is a square (0,1)^2
% author: Tommaso Vanzan
%-------------------------------------------------------------------------
clear all;
% The structure geo contains a function describing the geometry of the problem and the meshsize. 
%geo.fun='@square_foursub';% decomposition of the square (0,1)^2 into 4 subdomains
geo.fun='@square_twosub';% decomposition of the square (0,1)^2 into 2 subdomains
geo.regular=0; % flag to have a regular triangular mesh.( geo.regular needs to set geo.fun='square').

% The structure data contains information about the data of the problem,
% e.g. force term, diffusion coefficient, boundary conditions
data.f=@(x,y) 2*pi^2*sin(pi*x)*sin(pi*y); %force term 
data.c=@(x,y) 1; % diffusion term
data.Dirichlet_fun=@(x,y) 0; % boundary condition

% The structure para contains parameters of the problem such as the order
% of Gauss quadrature.
para.order_Gauss=2; %Order Gauss quadrature.
para.basis_type=201; %P1 FE space 201, P2 FE space 202

% Variables for Postprocessing
uex=@(x,y) sin(pi*x).*sin(pi*y); % Exact solution
hh=[1/2;1/8;1/16;1/32;1/64]; % Vector of mesh sizes
plt=0; % Set equal to 0 to disable plot, equal to 1 to plot solution
%=================================
%   Check Convergence
%=================================

for i=1:length(hh)
geo.h=hh(i);
[P,E,T,Pb,Tb,Eb,u]=Poisson_solver_2D_FEM(geo,data,para);%Solve Poisson Problem
for k=1:size(Pb,2)%evaluate exact solution on the mesh nodes.
	x=Pb(1,k); y=Pb(2,k);
	uexvec(k)=uex(x,y);
end
[error_L2(i),error_H1(i)]=compute_errors(u,uexvec,P,T,T((1:3),:),E,para,'FEM');%Compute L2 and H1 errors.
if plt==1
    figure(1)
    pdemesh(Pb,Eb,Tb,u); % Plot the numerical solution.
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