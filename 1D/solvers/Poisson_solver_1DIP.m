%-------------------------------------------------------------------------
% Poisson_solver_1DIP receives
% geo: structure contaning vertices domain (geo.a,geo.b) and mesh size geo.h.
% data: structure containing diffusion coefficient, force term, Dir./Neum.
% BC and labels vertices domain.
% basis_type: type of FEM space (linear,quadratic etc..).
% penalty: penalization parameter of interior penalty method
% Poisson_solver_1DIP returns:
% result(array): containing the solution.
% Mesh matrices: P,T (indipendent on basis_type)
% Pb,Tb: matrices with information on degrees of freedom (they do depend on
% basis_type)


% author: Tommaso Vanzan
%-------------------------------------------------------------------------
function [result,P,T,Pb,Tb]=Poisson_solver_1DIP(geo,basis_type,data,penalty)
[P,T,Pb,Tb,Eb]=generate_mesh_DG(geo,basis_type);% get information on the mesh and DOF according to basis_type
boundarynodes=generate_boundarynodes(Tb,data.left,data.right);% get the boundary nodes and their type.
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming DG spaces
matrixsize2=size(Pb,2);
A=assemble_matrix_1D(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,1,basis_type,1);%Assemble the matrix looping over the element for the 'classical' terms
%A=A+assemble_matrix_edges1DIP(data.c,P,T,Tb,Tb,Eb,matrixsize1,matrixsize2,basis_type,basis_type,penalty);
A=A+assemble_matrix_1D_IP_element(data.c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,basis_type,penalty); % Assemble tipical DG terms looping over ELEMENTS. 
b=assemble_rhs_1D(data.f,P,T,Tb,matrixsize2,basis_type,0); % Assemble the right hand side
result=A\b;%solve
end