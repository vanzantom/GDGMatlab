function [P,E,T,Pb,Tb,Eb,h,result]=Poisson_solver_2D_DG_IPH(h,basis_type,c,f,Dirichlet_fun,order_Gauss,alpha_coef)
% The function Poisson_solver_2D_DG_IPH receives the mesh size 'h', the
% 'basis_type', the diffusion function 'c', the RHS 'f', the boundary data
% 'Dirichlet_fun' the order of quadrature 'order_Gauss' and the coefficient
% of the DG penalization 'alpha_coef;

%=== Mesh generation
[P,T,E,Pb,Tb,Eb,h]=generate_mesh_2D(h,basis_type);
%=== Assembly of Stiffness matrix
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming FE spaces
matrixsize2=size(Pb,2);
A=assemble_matrix_2D(c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,1,0,basis_type,1,0,order_Gauss); %Assemble the matrix x gradient
A=A+assemble_matrix_2D(c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,0,1,basis_type,0,1,order_Gauss);  %Assemble the matrix y gradient
A=A+assemble_matrix_edges_IPH(c,P,T,Tb,Tb,Eb,matrixsize1,matrixsize2,basis_type,basis_type,order_Gauss,alpha_coef);%assemble matrix of the edge contributions for IPH.

%=== Assembly RHS
b=assemble_rhs_2D(f,P,matrixsize2,T,Tb,basis_type,0,0); % Assemble the right hand side
%=== Solution
result=A\b;
end
