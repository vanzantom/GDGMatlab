%% Example script to use the 1D FEM code
clear all;
a=0;% geometry
b=1;
f=@(x) -exp(x)*(cos(x)-2*sin(x)-x*cos(x)-x*sin(x)); %force term
c=@(x) exp(x); % diffusion term
Dirichlet_fun=@(x) (1-x)*0 + x*cos(1) ; % boundary condition. Label Dirichlet-1
Neumann_fun=@(x) cos(1)-sin(1); %Label Dirichlet-2
uex=@(x) x.*cos(x);%exact solution
plt=0;%plot variable
hh=[1/4;1/8;1/16;1/32;1/64;1/128];%array of mesh size
%=======convergence plot
%Linear
basis_type=101;% Linear P1
for i=1:length(hh)
h=hh(i);%mesh size
u=Poisson_solver_1D(a,b,h,basis_type,c,f,-1,-1,Dirichlet_fun,Neumann_fun);%solve the Problem
x=(0:h:1)';
uexvec=uex(x);%evaluate exact solution
normvec(i)=norm(uexvec-u,2)*sqrt(h);
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
% Quadratic
basis_type=102;%Quadratic Finite element
for i=1:length(hh)
h=hh(i);%mesh size
u2order=Poisson_solver_1D(a,b,h,basis_type,c,f,-1,-1,Dirichlet_fun,Neumann_fun);%solve problem
x2order=(0:h/2:1)';
uexvec=uex(x2order);
normvec(i)=norm(uexvec-u2order,2)*sqrt(h);
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
legend({'$$\|u_h-u\|_{L^2}-P1$$','$$h^2$$','$$\|u_h-u\|_{L^2}-P2$$','$$h^4}$$'},'interpreter','latex')
set(gca,'FontSize',12);