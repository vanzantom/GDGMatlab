%% Example script to implement OSM for HDG
clear all;close all;
f=@(x,y) 2*pi^2*sin(pi*x)*sin(pi*y); %force term
c=@(x,y) 1; % diffusion term
Dirichlet_fun=@(x,y) 0; % boundary condition
Neumann_fun=@(x) cos(1)-sin(1);
uex=@(x,y) sin(pi*x).*sin(pi*y);
hh=[1/2;1/8;1/16;1/32;1/64];%128
index= 4;%choose mesh size
order_Gauss=2;
plt=0;
alpha_coef=10;
basis_type=2010;
gamma=(1/2)*(1+sqrt(hh(index)));
NC=5; %number of functions in coarse level.
%======= Assembly matrices and compute solution
[P,E,T,Pb,Tb,Eb,hmax,AHDG,A1,A2,A1GAMMA,A2GAMMA,AGAMMA,b,nsub1,nsub2,ngamma,result]=Poisson_solver_2D_HDG(hh(index),basis_type,c,f,Dirichlet_fun,Neumann_fun,order_Gauss,alpha_coef);
uI=result(1:nsub2);
ulambda=result(nsub2+1:end);
B1=A1GAMMA'*(A1\A1GAMMA);
B2=A2GAMMA'*(A2\A2GAMMA);
S=AGAMMA-B1-B2;
Shat=[AGAMMA-B1,-B2;-B1, AGAMMA-B2];
%norm(sort(abs(eig(full(Shat))))-sort(abs([eig(full(S));eig(full(AGAMMA))]))) %check Theorem 1 of the notes
G1=(AGAMMA-B1)\B2;
G2=(AGAMMA-B2)\B1;
G=[0*AGAMMA,G1;G2,0*AGAMMA];
G1bar=B1*inv(AGAMMA-B1);
G2bar=B2*inv(AGAMMA-B2);
Gbar=[0*AGAMMA,G1bar;G2bar,0*AGAMMA];
Sbar=eye(2*ngamma)-Gbar;
G1tilde=(gamma*AGAMMA-B1)\(B2-(1-gamma)*AGAMMA);
G2tilde=(gamma*AGAMMA-B2)\(B1-(1-gamma)*AGAMMA);
Gtilde=[0*AGAMMA,G1tilde;G2tilde,0*AGAMMA];
Stilde=eye(2*ngamma)-Gtilde;
%======== convergence one level methods.
vg=eig(full(G));
max(abs(vg))% check spectral radius OSM is smaller than one.
[Vbar,Dbar]=eig(full(Gbar));
[vgbar,ibar]=sort(abs(diag(Dbar)),'descend');
max(vgbar)% check spectral radius OSM on Robin traces is smaller than one
[Vtilde,Dtilde]=eig(full(Gtilde));
[vgtilde,itilde]=sort(abs(diag(Dtilde)),'descend');
max(vgtilde)
%======= Two-level methods
Pbar=Vbar(:,ibar(1:NC));
Rbar=Pbar';
G2lbar=(eye(2*ngamma)-Pbar*((Rbar*Sbar*Pbar)\(Rbar*Sbar)))*Gbar;
[V2lbar,D2lbar]=eig(full(G2lbar));
[vg2lbar,i2lbar]=sort(abs(diag(D2lbar)),'descend');
max(vg2lbar)
Ptilde=Vtilde(:,itilde(1:NC));
Rtilde=Ptilde';
G2ltilde=(eye(2*ngamma)-Ptilde*((Rtilde*Stilde*Ptilde)\(Rtilde*Stilde)))*Gtilde;
[V2ltilde,D2ltilde]=eig(full(G2ltilde));
[vg2ltilde,i2ltilde]=sort(abs(diag(D2ltilde)),'descend');
max(vg2ltilde)

