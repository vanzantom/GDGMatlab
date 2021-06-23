%-------------------------------------------------------------------------
% example_HDG is an example of script using Poisson_solver_2D_HDG function to solve
% a Poisson problem on the geometry specified with Dirichlet BC on the boundary and P1 FE
% space.
% Further, it checks the order of convergence.
% The bilinear form explicity includes a substructured variable `\lambda',
% see eq(2.43) Soheil's thesis.

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

% The structure data contains information about the data of the problem,
% e.g. force term, diffusion coefficient, boundary conditions
data.f=@(x,y) 2*pi^2*sin(pi*x)*sin(pi*y); %force term
data.c=@(x,y,El) 1; % diffusion term
data.Dirichlet_fun=@(x,y) 0; % boundary condition
data.Neumann_fun=@(x) cos(1)-sin(1);

% The structure para contains parameters of the problem such as the order
% of Gauss quadrature.
para.order_Gauss=2;
para.alpha_coef=10;
para.basis_type=2010;

% Variables for Postprocessing
uex=@(x,y) sin(pi*x).*sin(pi*y);
hh=[1/2;1/8;1/16;1/32];%128
plt=0; % Set equal to 0 to disable plot, equal to 1 to plot solution, equal to 2 to plot mesh elements.
       % Set equal to 3 to plot decomposition into subdomains
%=======convergence plot
%Linear

for i=1:length(hh)
    geo.h=hh(i);
    [P,E,T,Pb,Tb,Eb,hmax,AHDG,AII,Aii,AIGAM,AGAM,b,Nsub,ngamma,indexGamma,result]=Poisson_solver_2D_HDG(geo,data,para);
    hvec(i)=hmax;
    for k=1:size(Pb,2)-ngamma %compute exact solution on the mesh nodes
        x=Pb(1,k); y=Pb(2,k);
        uexvec(k)=uex(x,y);
    end
    u=result(1:Nsub(end));
    if plt~=0
        visualize_solution_mesh(T,Pb,Tb,u,Nsub,0)
        pause
        close;
    end
     [error_L2(i),error_H1(i),error_DG(i)]=compute_errors(u,uexvec,P,T,Tb,Eb,para,'DG');% compute different error norms
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