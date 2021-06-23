%-------------------------------------------------------------------------
% Poisson_solver_2D_DG_IP receives
% geo: structure containing a function handle describing the geometry and
% the meshsize geo.h
% data: structure containing diffusion coefficient, force term, Dir BC.
% para: structure containing several parameters such as thebasis_type(type
%       of DG space, e.g linear,quadratic etc..), order of Gauss quadrature and
%       DG penalization coefficient.
% Poisson_solver_2D_DG_IP returns:
% result(array): containing the solution
% P,E,T: mesh matrices.
% Pb,Tb,Eb: matrices with data on Degrees of freedom (see
%           generate_mesh_2D for further information).
% h:         minimum length of an edge in the mesh

% author: Tommaso Vanzan
%-------------------------------------------------------------------------



function [P,E,T,Pb,Tb,Eb,h,result]=Poisson_solver_2D_DG_IP(geo,data,para)
% The function Poisson_solver_2D_DG_IP receives the mesh size 'h', the
% 'basis_type', the diffusion function 'c', the RHS 'f', the boundary data
% 'Dirichlet_fun' the order of quadrature 'order_Gauss' and the coefficient
% of the DG penalization 'alpha_coef;

%=== Mesh generation
[P,T,E,Pb,Tb,Eb,h]=generate_mesh_2D(geo,para.basis_type);
%=== Assembly of Stiffness matrix
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming FE spaces
matrixsize2=size(Pb,2);
A=assemble_matrix_2D(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,para.basis_type,1,0,para.basis_type,1,0,para.order_Gauss); %Assemble the matrix x gradient
A=A+assemble_matrix_2D(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,para.basis_type,0,1,para.basis_type,0,1,para.order_Gauss);  %Assemble the matrix y gradient
A=A+assemble_matrix_edges_IP(data.c,P,T,Tb,Tb,Eb,matrixsize1,matrixsize2,para.basis_type,para.basis_type,para);%assemble the matrix of the edge contributions for IP.
%=== Assembly RHS% 
b=assemble_rhs_2D(data.f,P,T,Tb,matrixsize2,para.basis_type,0,0); % Assemble the right hand side
%=== Solution
result=A\b;
end
