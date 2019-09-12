function [result,P,T,Pb,Tb]=Poisson_solver_1DIPH(a,b,h,basis_type,c,f,a_label,b_label,Dirichlet_fun,Neumann_fun,penalty)
[P,T,Pb,Tb,Eb]=generate_mesh_DG(a,b,h,basis_type);% get information on the mesh and DOF according to basis_type
boundarynodes=generate_boundarynodes(Tb,a_label,b_label);% get the boundary nodes and their type.
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming DG spaces
matrixsize2=size(Pb,2);
A=assemble_stiffness_matrix_1D_DG(c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,1,basis_type,1); %Assemble the matrix loop over the element for the 'classical' terms
%A=A+assemble_matrix_edges1DIP(c,P,T,Tb,Tb,Eb,matrixsize1,matrixsize2,basis_type,basis_type,penalty);%Assemble the penalization contribution looping over the edges
A=A+assemble_matrix_1D_IPH_element(c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,basis_type,penalty);%Assemble the penalization contribution looping over the element
b=assemble_rhs_1D(f,P,matrixsize2,T,Tb,basis_type,0); % Assemble the right hand side
result=A\b;%solve
end
