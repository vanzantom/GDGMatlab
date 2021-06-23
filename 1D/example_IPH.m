%-------------------------------------------------------------------------
% example_IPH is an example of script using Poisson_solver_1DIPH() function to solve
% a Poisson problem on the interval (a,b) with Dirichlet/Neumann BC with
% linear or quadratic DG spaces using HYBRIDAZABLE Interior penalty method.
% Further, it checks the order of convergence.

% Attention: Neumann boundary condition not implemented jet!

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
clear all;
addpath('lib1D')
addpath('solvers')
% The structure geo contains the extrema of the interval Omega=(a,b).
geo.a=0;
geo.b=1;
% The structure data contains information about the data of the problem,
% e.g. force term, diffusion coefficient, boundary conditions
data.f=@(x) 2*pi*exp(x).*(2*pi*sin(2*pi*x)-cos(2*pi*x)); %force term
data.c=@(x) exp(x); % diffusion term
data.Dirichlet_fun=@(x) 0 ; % boundary condition
data.Neumann_fun=@(x) cos(1)-sin(1); %DO NOT USE NEUMANN
data.left=-1;%label variable for left node. If data.left=-1, then we have Dirichlet B.C. on x=a. If data.left=-1, then we have Neumann B.C. on x=a
data.right=-1;%label variable for right node. 
% Penalization Constant
alpha=10; %costant for penalization
% Variables for Postprocessing
uex=@(x) sin(2*pi*x);
uex_der=@(x) 2*pi*cos(2*pi*x);
hh=[1/4;1/8;1/16;1/32;1/64;1/128];
plt=0;
%=================================
%   IPH method with linear FE space
%=================================
basis_type=101;
for i=1:length(hh)
geo.h=hh(i); %mesh size
penalty=alpha/geo.h; %penalty term 
[u,P,T,Pb,Tb]=Poisson_solver_1DIPH(geo,basis_type,data,penalty);
x=(0:geo.h:1)'; %discrete mesh
uexvec=uex(x); %comput exact solution
if plt==1
    figure(1); hold on;
    j=1;
 for k=1:length(P)-1  
     v=P(k:k+1);
 plot(v,u(j:j+1),'b') %plot numerical solution in the current element
 j=j+2;
 end
 plot(x,uexvec,'r')
 pause
 close
end
err(i)=errornormaDG(P,T,Tb,uex_der,u,basis_type,penalty); %compute error in DG norm
end
figure(2)
loglog(hh,err,'b*-')
hold on
loglog(hh,hh.^1,'r')
%=================================
%   IPH method with quadratic FE space
%=================================
basis_type=102;
for i=1:length(hh)
geo.h=hh(i);
penalty=alpha/geo.h;
[u2order,P,T,Pb,Tb]=Poisson_solver_1DIPH(geo,basis_type,data,penalty);
x2order=(0:geo.h/2:1)'; %discrete mesh
uexvec=uex(x2order); %exact solution
if plt==1
    figure(1); hold on;
    j=1;
 for k=1:length(P)-1
     v=Pb(j:j+2);
 plot(v,u2order(j:j+2),'b') %plot solution in the current element
 j=j+3;
 end
 plot(x2order,uexvec,'r')
 pause
close
end
err(i)=errornormaDG(P,T,Tb,uex_der,u2order,basis_type,penalty); %compute error in DG norm
end
figure(2)
loglog(hh,err,'m*-')
hold on
loglog(hh,hh.^2,'c')
grid on;
xlabel('h')
ylabel('errors')
legend({'$$\|u_h-u\|_{DG}-P1$$','$$h$$','$$\|u_h-u\|_{DG}-P2$$','$$h^2$$'},'interpreter','latex')
set(gca,'FontSize',12);
