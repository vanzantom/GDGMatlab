function result=Poisson_solver_1D(a,b,h,basis_type,c,f,a_label,b_label,Dirichlet_fun,Neumann_fun)

[P,T,Pb,Tb]=generate_mesh_FEM(a,b,h,basis_type);% get information on the mesh and DOF according to basis_type
boundarynodes=generate_boundarynodes(Tb,a_label,b_label);% get the boundary nodes and their type.
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming FE spaces
matrixsize2=size(Pb,2);
A=assemble_matrix_1D(c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,1,basis_type,1); %Assemble the matrix
b=assemble_rhs_1D(f,P,matrixsize2,T,Tb,basis_type,0); % Assemble the right hand side
[A,b]=treat_boundary(A,b,boundarynodes,P,Pb,c,Dirichlet_fun,Neumann_fun);% Strong imposition of Dirichlet conditions.
result=A\b;%solve
end
