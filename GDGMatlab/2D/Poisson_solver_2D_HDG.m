function [P,E,T,Pb,Tb,Eb,h,A,A1,A2,A1GAMMA,A2GAMMA,AGAMMA,b,nsub1,nsub2,ngamma,result]=Poisson_solver_2D_HDG(h,basis_type,c,f,Dirichlet_fun,Neumann_fun,order_Gauss,alpha_coef)
% The function Poisson_solver_2D_HDG receives the mesh size 'h', the
% 'basis_type', the diffusion function 'c', the RHS 'f', the boundary data
% 'Dirichlet_fun' and 'Neumann_fun', the order of quadrature order Gauss
% and the penalization coefficient alpha_coef.
% The output contains the mesh information matrices P,E,T,Pb,Tb,Eb,h and the
% matrices A,A1,A2,A1GAMMA,A2GAMMA,AGAMMA. See eq (3.3) Soheil's thesis.
% available at https://archive-ouverte.unige.ch/authors/view/74979.
% It also computes the right hand side vector b and the solution stored in
% the vector result.
%=== Mesh generation
[P,T,E,Pb,Tb,Eb,h,nsub1,nsub2,ngamma]=generate_mesh_2D_HDG(h,basis_type);
%=== Assembly of Stiffness matrix
matrixsize1=size(Pb,2);% number of DOF. Supposing conforming FE spaces
matrixsize2=size(Pb,2);


A=assemble_matrix_2D(c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,1,0,basis_type,1,0,order_Gauss); %Assemble the matrix x gradient
A=A+assemble_matrix_2D(c,P,T,Tb,Tb,matrixsize1,matrixsize2,basis_type,0,1,basis_type,0,1,order_Gauss);  %Assemble the matrix y gradient
A=A+assemble_matrix_edges_HDG(c,P,T,Pb,Tb,Tb,Eb,matrixsize1,matrixsize2,basis_type,basis_type,order_Gauss,alpha_coef);
A=A+2*assemble_matrix_HDG_AGamma(c,P,Eb,Eb,Eb,nsub2,ngamma,basis_type,basis_type,order_Gauss,alpha_coef);
AIGAM=assemble_matrix_HDG_AIGamma(c,P,T,Tb,Eb,Eb,nsub2,ngamma,basis_type,basis_type,order_Gauss,alpha_coef);
A=A+AIGAM+AIGAM';
%=== Assembly RHS
b=assemble_rhs_2D(f,P,matrixsize2,T,Tb,basis_type,0,0); % Assemble the right hand side
%=== Solution
result=A\b;
% post processing
A1=A(1:nsub1,1:nsub1);
A2=A(nsub1+1:nsub2,nsub1+1:nsub2);
AGAMMA=A(end-ngamma+1:end,end-ngamma+1:end);
A1GAMMA=A(1:nsub1,end-ngamma+1:end);
A2GAMMA=A(nsub1+1:nsub2,end-ngamma+1:end);
end