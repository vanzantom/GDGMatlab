%-------------------------------------------------------------------------
% assemble_matrix_1D receives 
% coe_fun: the diffusion coefficient 
% mesh matrices: T,P 
% Tb_trial,Tb_test: matrices with information on trial and test finite element space 
% matrixsize1,matrixsize2: size of the stiffness matrix 
% basis_type_trial,basis_type_test: finite element space type trial/test
% der_trial,der_test: order of derivatives in the bilinea form .
% assemble_matrix_1D returns: 
% A (matrix) whose entry A_{i,j}= \int ((\partial_x)^der_trial \phi_j)((\partial_x)^der_test \phi_i) 


% author: Tommaso Vanzan
%-------------------------------------------------------------------------
function A=assemble_matrix_1D(coe_fun,P,T,Tb_trial,Tb_test,matrixsize1,matrixsize2,basis_type_trial,der_trial,basis_type_test,der_test)
A=sparse(matrixsize1,matrixsize2);%create sparse matrix
number_of_elements=size(T,2); %number of Elements
number_of_local_basis_trial=size(Tb_trial,1); %number of trial local basis function on a single element
number_of_local_basis_test=size(Tb_test,1); %number of test local basis function on a single element
for n=1:number_of_elements %loop over number of elements
    vertices=P(:,T(:,n)); % mesh information then I use P and T
    [Gauss_nodes,Gauss_weights]=generate_Gauss(vertices,2); % create Gauss nodes on the element. 1 trapezoidal 2 Simpsons rule. It contains already the Jacobian
    for alpha=1:number_of_local_basis_trial
        for beta=1:number_of_local_basis_test %from reference-> local -> global assembly.
            int_value=Gauss_quadrature_1D_volume_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,alpha,der_trial,basis_type_test,beta,der_test,n);
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+int_value;
        end
    end
end
end