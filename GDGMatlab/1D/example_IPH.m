%% Example script to use the 1D DG Symmetric HYBRID internal penalty method
clear all;
a=0;% geometry
b=1;
f=@(x) 2*pi*exp(x).*(2*pi*sin(2*pi*x)-cos(2*pi*x)); %force term
c=@(x) exp(x); % diffusion term
Dirichlet_fun=@(x) 0 ; % boundary condition
Neumann_fun=@(x) cos(1)-sin(1); %DO NOT USE NEUMANN
uex=@(x) sin(2*pi*x);
uex_der=@(x) 2*pi*cos(2*pi*x);
hh=[1/4;1/8;1/16;1/32;1/64;1/128];
alpha=10; %costant for penalization
plt=0;
%=======convergence plot
%Linear
basis_type=101;
for i=1:length(hh)
h=hh(i); %mesh size
penalty=alpha/h; %penalty term 
[u,P,T,Pb,Tb]=Poisson_solver_1DIPH(a,b,h,basis_type,c,f,-1,-1,Dirichlet_fun,Neumann_fun,penalty);
x=(0:h:1)'; %discrete mesh
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
% Quadratic
basis_type=102;
for i=1:length(hh)
h=hh(i);
penalty=alpha/h;
[u2order,P,T,Pb,Tb]=Poisson_solver_1DIPH(a,b,h,basis_type,c,f,-1,-1,Dirichlet_fun,Neumann_fun,penalty);
x2order=(0:h/2:1)'; %discrete mesh
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
