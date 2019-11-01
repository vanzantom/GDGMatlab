%% Example script to use the HDG ( based on lectures by Xiaoming HE)
clear all;close all;
f=@(x,y) 2*pi^2*sin(pi*x)*sin(pi*y); %force term
c=@(x,y) 1; % diffusion term
Dirichlet_fun=@(x,y) 0; % boundary condition
Neumann_fun=@(x) cos(1)-sin(1);
uex=@(x,y) sin(pi*x).*sin(pi*y);
hh=[1/2;1/8;1/16;1/32];%128
order_Gauss=2;
plt=0;
alpha_coef=10;
%=======convergence plot
%Linear
basis_type=2010;
for i=1:length(hh)

h=hh(i);
[P,E,T,Pb,Tb,Eb,hmax,AHDG,A1,A2,A1GAMMA,A2GAMMA,AGAMMA,b,nsub1,nsub2,ngamma,result]=Poisson_solver_2D_HDG(h,basis_type,c,f,Dirichlet_fun,Neumann_fun,order_Gauss,alpha_coef);
hvec(i)=hmax;
for k=1:size(Pb,2)-ngamma %compute exact solution on the mesh nodes
	x=Pb(1,k); y=Pb(2,k);
	uexvec(k)=uex(x,y);
end
u=result(1:nsub2);
if plt==1
    visualize_solution_mesh(T,Pb,Tb,u,0)
    pause
    close;
end
 [error_L2(i),error_H1(i),error_DG(i)]=compute_errors(u,uexvec,P,T,Tb,Eb,basis_type,order_Gauss,alpha_coef,'DG');% compute different error norms
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