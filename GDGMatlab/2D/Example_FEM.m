%% Example script to use the 2D FEM code.
clear all;
f=@(x,y) 2*pi^2*sin(pi*x)*sin(pi*y); %force term 
c=@(x,y) 1; % diffusion term
Dirichlet_fun=@(x,y) 0; % boundary condition
uex=@(x,y) sin(pi*x).*sin(pi*y); % Exact solution
hh=[1/4;1/8;1/16;1/32;1/64]; % Vector of mesh sizes
order_quad=2; %Order of quadrature.
plt=0;
%=======convergence plot
%Linear
basis_type=201; %P1 finite element.
for i=1:length(hh)
h=hh(i);
[P,E,T,u]=Poisson_solver_2D_FEM(h,basis_type,c,f,Dirichlet_fun,order_quad);%Solve Poisson Problem
for k=1:size(P,2)%evaluate exact solution on the mesh nodes.
	x=P(1,k); y=P(2,k);
	uexvec(k)=uex(x,y);
end
[error_L2(i),error_H1(i)]=compute_errors(u,uexvec,P,T,T((1:3),:),E,basis_type,order_quad,0,'FEM');%Compute L2 and H1 errors.
if plt==1
    figure(1)
    pdemesh(P,E,T,u); % Plot the numerical solution.
%   figure(2)
%   pdemesh(P,E,T,u, 'mesh', 'on');% Plot the mesh
%   pdemesh(P,E,T,uexvec,); % Plot the exact solution.
    pause
end
end
loglog(hh,error_L2,'b*-')
hold on
loglog(hh,hh.^2,'r')
loglog(hh,error_H1,'k*-')
loglog(hh,hh,'m')
grid on
xlabel('h')
ylabel('errors')
legend({'$$\|u_h-u\|_{L^2}$$','$$h^{2}$$','$$\|u_h-u\|_{H_1}$$','$$h$$'},'interpreter','latex')
set(gca,'FontSize',12);