%% Example script to use the DG IP AND IPH
%clear all;close all;
f=@(x,y) 2*pi^2*sin(pi*x)*sin(pi*y); %force term
c=@(x,y) 1; % diffusion term
Dirichlet_fun=@(x,y) 0; % boundary condition
uex=@(x,y) sin(pi*x).*sin(pi*y);%exact solution
hh=[1/2;1/8;1/16;1/32;1/64];%128
order_Gauss=2;%order of quadrature
plt=0;%variable to control the plotting
alpha_coef=10;
%=======convergence plot
%Linear
basis_type=2010; %2010 is the code for P1 discontinous Finite element
for i=1:length(hh)

h=hh(i);
[P,E,T,Pb,Tb,Eb,hmax,u]=Poisson_solver_2D_DG_IP(h,basis_type,c,f,Dirichlet_fun,order_Gauss,alpha_coef);
%[P,E,T,Pb,Tb,Eb,hmax,u]=Poisson_solver_2D_DG_IPH(h,basis_type,c,f,Dirichlet_fun,order_Gauss,alpha_coef);
hvec(i)=hmax;%maximum length of an edge
for k=1:size(Pb,2)%computing the exact solution in the degrees of freedom
	x=Pb(1,k); y=Pb(2,k);
	uexvec(k)=uex(x,y);
end
if plt==1
    visualize_solution_mesh(T,Pb,Tb,u,0)%plot
    pause
    close;
end
 [error_L2(i),error_H1(i),error_DG(i)]=compute_errors(u,uexvec,P,T,Tb,Eb,basis_type,order_Gauss,alpha_coef,'DG');%compute L2,H1 and DG errors.
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