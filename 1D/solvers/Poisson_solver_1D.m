%-------------------------------------------------------------------------
% Poisson_solver_1D receives
% geo: structure contaning vertices domain (geo.a,geo.b) and mesh size geo.h.
% data: structure containing diffusion coefficient, force term, Dir./Neum.
% BC and labels vertices domain.
% basis_type: type of FEM space (linear,quadratic etc..).
% Poisson_solver_1D returns:
% result(array): containing the solution

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
function result=Poisson_solver_1D(geo,basis_type,data)
[P,T,Pb,Tb]=generate_mesh_FEM(geo,basis_type);% get information on the mesh and DOF according to basis_type
boundarynodes=generate_boundarynodes(Tb,data.left,data.right);% get the boundary nodes and their type.
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming FE spaces
matrixsize2=size(Pb,2);
A=assemble_matrix_1D(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,1,basis_type,1); %Assemble the matrix
b=assemble_rhs_1D(data.f,P,T,Tb,matrixsize2,basis_type,0); % Assemble the right hand side
%[A,b]=treat_boundary(A,b,boundarynodes,P,Pb,data.c,data.Dirichlet_fun,data.Neumann_fun);% Strong imposition of Dirichlet conditions.
[A,b,int]=treat_boundary_symmetric(A,b,boundarynodes,Pb,data.c,data.Dirichlet_fun,data.Neumann_fun);
result=A\b;%solve
end