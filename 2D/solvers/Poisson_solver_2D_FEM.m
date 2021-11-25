%-------------------------------------------------------------------------
% Poisson_solver_2D receives
% geo: structure containing a function handle describing the geometry and
% the meshsize geo.h
% data: structure containing diffusion coefficient, force term, Dir BC.
% basis_type: type of FEM space (linear,quadratic etc..).
% Poisson_solver_2D returns:
% result(array): containing the solution
% P,E,T: mesh matrices.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------


function [P,E,T,Pb,Tb,Eb,result]=Poisson_solver_2D_FEM(geo,data,para)
% The function Poisson_solver_2D receives the mesh size 'h', the
% 'basis_type', the diffusion function 'c', the RHS 'f', the boundary data
% 'Dirichlet_fun'

%=== Mesh generation
[P,T,E,Pb,Tb,Eb]=generate_mesh_2D(geo,para.basis_type);
%=== Identification of Boundary
boundaryedges=generate_boundaryedges(Tb,Eb,E,data.label,para.basis_type);% get the boundary nodes and their type.
%=== Assembly of Stiffness matrix
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming FE spaces
matrixsize2=size(Pb,2);
A=assemble_matrix_2D(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,para.basis_type,1,0,para.basis_type,1,0,para.order_Gauss); %Assemble the matrix of the contribution of x gradient
A=A+assemble_matrix_2D(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,para.basis_type,0,1,para.basis_type,0,1,para.order_Gauss);  %Assemble the matrix of the contribution of  y gradient
%=== Assembly RHS
b=assemble_rhs_2D(data.f,P,T,Tb,matrixsize2,para.basis_type,0,0); % Assemble the right hand side
%=== Take Care of Boundary
[A,b,ug,int,out]=treat_boundary_symmetric(A,b,boundaryedges,Pb,Tb,Tb,para,data.Dirichlet_fun,data.Neumann_fun,data.Robin_fun);% Strong imposition of Dirichlet conditions.
%=== Solution
x=A\b;
result=zeros(size(Pb,2),1);
result(int)=x;
result(out)=ug(out);
end
