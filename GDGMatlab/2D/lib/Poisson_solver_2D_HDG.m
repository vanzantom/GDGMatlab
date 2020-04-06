%-------------------------------------------------------------------------
% Poisson_solver_2D_DG_HDG receives
% geo: structure containing a function handle describing the geometry and
% the meshsize geo.h
% data: structure containing diffusion coefficient, force term, Dir BC.
% para: structure containing several parameters such as the basis_type(type
%       of DG space, e.g linear,quadratic etc..), order of Gauss quadrature and
%       DG penalization coefficient.
% Poisson_solver_2D_DG_IP returns:
% result(array): containing the solution
% P,E,T: mesh matrices.
% Pb,Tb,Eb: matrices with data on Degrees of freedom (see
%           generate_mesh_2D for further information).
% h:         minimum length of an edge in the mesh
% A:        HDG stiffness matrix
% AII:      submatrix of A related to interior degrees of freedom
% AIGAM:    submatrix of A related to the coupling interior/substructured
%           degrees of freedom
% AGAM    subsmatrix of A for substructured DOFs, See eq (3.3) in Soheil's thesis for more details
%            available at https://archive-ouverte.unige.ch/authors/view/74979.
% b         b right hand side vector
%Nsub       vector. Nsub(j) contains last index DoF in subdomain j.
% ngamma:  total number of degrees of freedom on interface Gamma between
%            subdomains
% indexGamma: vector containing indexes DOF on GAMMA.


% author: Tommaso Vanzan
%-------------------------------------------------------------------------


function [P,E,T,Pb,Tb,Eb,h,A,AII,Aii,AIGAM,AGAM,b,Nsub,ngamma,indexGamma,result]=Poisson_solver_2D_HDG(geo,data,para)
%=== Mesh generation
[P,T,E,Pb,Tb,Eb,h,Nsub,ngamma,indexGamma]=generate_mesh_2D_HDG(geo,para.basis_type);
%=== Assembly of Stiffness matrix
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming FE spaces
matrixsize2=size(Pb,2);


A=assemble_matrix_2D(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,para.basis_type,1,0,para.basis_type,1,0,para.order_Gauss); %Assemble the matrix x gradient
A=A+assemble_matrix_2D(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,para.basis_type,0,1,para.basis_type,0,1,para.order_Gauss);  %Assemble the matrix y gradient
A=A+assemble_matrix_edges_HDG(data.c,P,T,Pb,Tb,Tb,Eb,matrixsize1,matrixsize2,para.basis_type,para.basis_type,para);
A=A+2*assemble_matrix_HDG_AGamma(data.c,P,Eb,Eb,Eb,Nsub,ngamma,para.basis_type,para.basis_type,para);
AIGAM=assemble_matrix_HDG_AIGamma(data.c,P,T,Tb,Eb,Eb,Nsub,ngamma,para.basis_type,para.basis_type,para);
A=A+AIGAM+AIGAM';
%=== Assembly RHS
b=assemble_rhs_2D(data.f,P,T,Tb,matrixsize2,para.basis_type,0,0); % Assemble the right hand side
%=== Solution
result=A\b;
% post processing to assemble
Nsubdomains=length(Nsub);
Aii{1}=A(1:Nsub(1),1:Nsub(1));
for i=2:Nsubdomains
    Aii{i}=A(Nsub(i-1)+1:Nsub(i),Nsub(i-1)+1:Nsub(i));
end
AII=A(1:Nsub(end),1:Nsub(end));
AGAM=A(end-ngamma+1:end,end-ngamma+1:end);
AIGAM=A(1:Nsub(end),end-ngamma+1:end);
end