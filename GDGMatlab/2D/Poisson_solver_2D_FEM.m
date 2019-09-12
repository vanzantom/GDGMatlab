function [P,E,T,result]=Poisson_solver_2D_FEM(h,basis_type,c,f,Dirichlet_fun,order_Gauss)
% The function Poisson_solver_2D receives the mesh size 'h', the
% 'basis_type', the diffusion function 'c', the RHS 'f', the boundary data
% 'Dirichlet_fun'

%=== Mesh generation
[P,T,E,Pb,Tb,~]=generate_mesh_2D(h,basis_type);
%=== Identification of Boundary
boundarynodes=generate_boundarynodes(Pb);% get the boundary nodes and their type.
%=== Assembly of Stiffness matrix
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming FE spaces
matrixsize2=size(Pb,2);
A=assemble_matrix_2D(c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,1,0,basis_type,1,0,order_Gauss); %Assemble the matrix x gradient
A=A+assemble_matrix_2D(c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,0,1,basis_type,0,1,order_Gauss);  %Assemble the matrix y gradient
%=== Assembly RHS
b=assemble_rhs_2D(f,P,matrixsize2,T,Tb,basis_type,0,0); % Assemble the right hand side
%=== Take Care of Boundary
[A,b]=treat_boundary(A,b,boundarynodes,Pb,Dirichlet_fun);% Strong imposition of Dirichlet conditions.
%=== Solution
result=A\b;
end
