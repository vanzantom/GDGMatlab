%-------------------------------------------------------------------------
% example_DG is an example of script to solve
% a Poisson problem with Discontinous Galerkin method on the geometry specified with Dirichlet boundary condition on the
% boundary.
% The user can use either Interior Penalty method (IP) or Hybridizable
% Interior penalty (IPH) method by choosing the appropriare
% 'Poisson_solver'

% Further, it checks the order of convergence.

% Attention: Currently, the function generate_boundarynodes works only if
% the geometry is a square (0,1)^2

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
clear all;close all;
addpath('lib');
addpath('solvers');
% The structure geo contains a function describing the geometry of the problem and the meshsize. 
%geo.fun='@square_foursub';% decomposition of the square (0,1)^2 into 4 subdomains
geo.fun='@square_twosub';% decomposition of the square (0,1)^2 into 2 subdomains
geo.regular=0; % flag to have a regular triangular mesh.( geo.regular needs to set geo.fun='square').

% The structure data contains information about the data of the problem,
% e.g. force term, diffusion coefficient, boundary conditions
data.f=@(x,y) 2*pi^2*sin(pi*x)*sin(pi*y); %force term
data.c=@(x,y,El) 1; % diffusion term
data.Dirichlet_fun=@(x,y) 0; % boundary condition

% The structure para contains parameters of the problem such as the order
% of Gauss quadrature.
para.order_Gauss=2;%order of quadrature
para.alpha_coef=10; %DG penalization coefficient
para.basis_type=2010;%2010 is the code for P1 discontinous Finite element

% Variables for Postprocessing
uex=@(x,y) sin(pi*x).*sin(pi*y);%exact solution
hh=[1/2;1/8;1/16;1/32];% set of meshsize to check convergence
plt=1; % Set equal to 0 to disable plot, equal to 1 to plot solution, equal to 2 to plot mesh elements.
       % Set equal to 3 to plot decomposition into subdomains
Nsub=1;

%=================================
%   Check Order Convergence
%=================================

for i=1:length(hh)
geo.h=hh(i);
[P,E,T,Pb,Tb,Eb,hmax,u]=Poisson_solver_2D_DG_IP(geo,data,para);
%[P,E,T,Pb,Tb,Eb,hmax,u]=Poisson_solver_2D_DG_IPH(geo,data,para);

hvec(i)=hmax;%maximum length of an edge
for k=1:size(Pb,2)%computing the exact solution in the degrees of freedom
	x=Pb(1,k); y=Pb(2,k);
	uexvec(k)=uex(x,y);
end
if plt~=0
    visualize_solution_mesh(T,Pb,Tb,u,Nsub,plt)%plot
    pause
    close;
end
 [error_L2(i),error_H1(i),error_DG(i)]=compute_errors(u,uexvec,P,T,Tb,Eb,para,'DG');%compute L2,H1 and DG errors.
end
loglog(hh,error_L2,'b*-')
hold on
loglog(hh,hh.^2,'r')
loglog(hh,error_H1,'k*-')
loglog(hh,error_DG,'c*-')
loglog(hh,hh,'m')
grid on
xlabel('h')
ylabel('errors')
legend({'$$\|u_h-u\|_{L^2}$$','$$h^{2}$$','$$\|u_h-u\|_{H_1}$$','$$\|u_h-u\|_{DG}$$','$$h$$'},'interpreter','latex')
set(gca,'FontSize',12);