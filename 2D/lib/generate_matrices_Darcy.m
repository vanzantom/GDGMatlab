function [ANeumann,A,ARobin,ADir,M_GG,b,e,Pb,Tb,P,T,E,boundary]=generate_matrices_Darcy(geo,data,para)
%generate_matrices_Darcy receives the data structure geo,data and para.
% It returns the matrices
%A: stiffness matrix with Neumann BC everywhere.
%ADir: stifnesss matrix with Dirichlet BC where imposed and Neumann on
%other edges (i.e. edge with label -3 are mapped to -2)
%ARobin: stifness matrix with Dirichlet BC and Neumann and Robin where imposed
%M_GG: local mass matrix on the interior DOF of the interface.
%b: force term which contains Neumann BC and Dirichlet BC.
% e: vector which contains indexes DOfs in interior edge.
%data structure matrices: Pb,Tb for FE space
% data structure for mesh; P,T,E.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
%========================================================================
geo.fun=geo.funDarcy;
[P,T,E,Pb,Tb,Eb]=generate_mesh_2D(geo,para.basis_type);
%=== Identification of Boundary
boundary=generate_boundaryedges(Tb,Eb,E,data.label,para.basis_type);% get the boundary edges and their type.
%=== Assembly of Stiffness matrix
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming FE spaces
matrixsize2=size(Pb,2);
ANeumann=assemble_matrix_2D(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,para.basis_type,1,0,para.basis_type,1,0,para.order_Gauss); %Assemble the matrix of the contribution of x gradient
ANeumann=ANeumann+assemble_matrix_2D(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,para.basis_type,0,1,para.basis_type,0,1,para.order_Gauss);  %Assemble the matrix of the contribution of  y gradient
%A now is a Neumann matrix.
b=assemble_rhs_2D(data.g1,P,T,Tb,matrixsize2,para.basis_type,1,0);%compute possible f term
b=b+assemble_rhs_2D(data.g2,P,T,Tb,matrixsize2,para.basis_type,0,1);%compute possible f term


[A,ARobin,ADir,b,M_GG,e]=impose_boundary_Darcy(ANeumann,b,boundary,Pb,Tb,Tb,para,data.Dirichlet_fun,data.Neumann_fun);% create a matrix which has Dirichlet conditions on three edge and Neumann on the top edge

end