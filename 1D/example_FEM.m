%-------------------------------------------------------------------------
% example_FEM is an example of script using Poisson_solver_1D function to solve
% a Poisson problem on the interval (a,b) with Dirichlet/Neumann BC with P1
% or P2 FE spaces.
% Further, it checks the order of convergence.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------

clear all; close all;
addpath('lib1D')
addpath('solvers')
% The structure geo contains the extrema of the interval Omega=(a,b)
geo.a=0;
geo.b=1;
% The structure data contains information about the data of the problem,
% e.g. force term, diffusion coefficient, boundary conditions
data.f=@(x) -exp(x)*(cos(x)-2*sin(x)-x*cos(x)-x*sin(x)); 
data.c=@(x) exp(x); % diffusion term
data.Dirichlet_fun=@(x) (1-x)*0 + x*cos(1) ; % boundary condition. % Label Dirichlet=-1
data.Neumann_fun=@(x) cos(1)-sin(1); %Label Neumann=-2
data.left=-1;%label variable for left node. If data.left=-1, then we have Dirichlet B.C. on x=a. If data.left=-1, then we have Neumann B.C. on x=a
data.right=-1;%label variable for right node. 
% Variables for Postprocessing
uex=@(x) x.*cos(x);%exact solution
plt=0;%plot variable
hh=[1/4;1/8;1/16;1/32;1/64;1/128];%array of mesh size
%=================================
%   Linear
%=================================
basis_type=101;% Linear P1
for i=1:length(hh)
geo.h=hh(i);%mesh size
u=Poisson_solver_1D(geo,basis_type,data);%solve the Problem
x=(0:geo.h:1)';
uexvec=uex(x);%evaluate exact solution
normvec(i)=norm(uexvec-u,2)*sqrt(geo.h);
if plt==1
    plot(x,u)%plot approximation
    hold on
    plot(x,uexvec,'r')%plot exact solution
    pause
    close
end
end
loglog(hh,normvec,'b*-')
hold on
loglog(hh,hh.^2,'r')
%=================================
%   Quadratic
%=================================
basis_type=102;%Quadratic Finite element
for i=1:length(hh)
geo.h=hh(i);%mesh size
u2order=Poisson_solver_1D(geo,basis_type,data);%solve problem
x2order=(0:geo.h/2:1)';
uexvec=uex(x2order);
normvec(i)=norm(uexvec-u2order,2)*sqrt(geo.h);
if plt==1
    plot(x2order,u2order)
    hold on
    plot(x2order,uexvec,'r')
    pause
    close
end
end
loglog(hh,normvec,'k-*')
hold on
loglog(hh,hh.^4,'m')
grid on
xlabel('h')
ylabel('errors')
legend({'$$\|u_h-u\|_{L^2}-P1$$','$$h^2$$','$$\|u_h-u\|_{L^2}-P2$$','$$h^{4}$$'},'interpreter','latex')
set(gca,'FontSize',12);