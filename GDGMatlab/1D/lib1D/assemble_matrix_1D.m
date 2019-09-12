function A=assemble_matrix_1D(coe_fun,P,T,Tb_trial,Tb_test,matrixsize1,matrixsize2,basis_type_trial,der_trial,basis_type_test,der_test)
% The function assemble_matrix_1D computes the contribution of the intregrals in 1D with order 'der_trial','der_test' 
%=== It receives 
% coe_fun= Diffusion coefficient.
% P, T=> vertex matrix of the mesh
% Tb_trial, Tb_test=> DOFs of the Finite_element type on the elements.
% matrixsize1,matrixsiz2=> DOFS of the trial and test space.
% basis_type_trial/test=> FE space for trial and test
% der_(trial/test)=> order derivative in x of (test/trial) space.

A=sparse(matrixsize1,matrixsize2);%create sparse matrix
number_of_elements=size(T,2); %number of Elements
number_of_local_basis_trial=size(Tb_trial,1); %number of trial local basis function on a single element
number_of_local_basis_test=size(Tb_test,1); %number of test local basis function on a single element
for n=1:number_of_elements
    vertices=P(:,T(:,n)); % mesh information then I use P and T
    [Gauss_nodes,Gauss_weights]=generate_Gauss(vertices,2); % create Gauss nodes on the element. 1 trapezoidal 2 Simpsons rule.
    for alpha=1:number_of_local_basis_trial
        for beta=1:number_of_local_basis_test %from reference-> local -> global assembly.
            int_value=Gauss_quadrature_1D_volume_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,alpha,der_trial,basis_type_test,beta,der_test);
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+int_value;
        end
    end
end
