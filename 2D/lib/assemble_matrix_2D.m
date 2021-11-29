%-------------------------------------------------------------------------
% assemble_matrix_2D receives 
% coe_fun: the diffusion coefficient 
% mesh matrices: T,P 
% Tb_trial,Tb_test: matrices with information on trial and test finite element space 
% matrixsize1,matrixsize2: size of the stiffness matrix 
% basis_type_trial,basis_type_test: finite element space type trial/test
% der_x_trial,der_x_test, der_y_trial,der_y_test,: order of derivatives in the bilinea form
% order_gauss: Order of Gauss quadrature.
% assemble_matrix_2D returns: 
% A (matrix) whose entry A_{i,j}= \int (c(\partial_x)^der_x_trial (\partial_y)^der_y_trial \phi_j)((\partial_x)^der_x_test (\partial_y)^der_y_test \phi_j)


% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function A=assemble_matrix_2D(coe_fun,P,T,Tb_trial,Tb_test,matrixsize1,matrixsize2,basis_type_trial,der_x_trial,der_y_trial,basis_type_test,der_x_test,der_y_test,order_Gauss)
%=== Create sparse matrix and define parameters.
A=sparse(matrixsize1,matrixsize2);%create sparse matrix
number_of_elements=size(T,2); %number of Elements
number_of_local_basis_trial=size(Tb_trial,1); %number of trial local basis function on a single element
number_of_local_basis_test=size(Tb_test,1); %number of test local basis function on a single element

%=== Loop to assembly
for n=1:number_of_elements % Loop over the triangles.
    vertices=P(:,T((1:3),n)); % vertices of current Element. Element n, vertices T(:,n), Position P(:,T(:,n)).
    J=abs((vertices(1,2)-vertices(1,1))*(vertices(2,3)-vertices(2,1)) -(vertices(1,3)-vertices(1,1))*(vertices(2,2)-vertices(2,1)));% Compute Jacobian of the linear transformation
    
    [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,order_Gauss); % create Gauss nodes on the element. 1 baricentri 2 order accuracy.
    for alpha=1:number_of_local_basis_trial %loop over trial functions on the element for trial FE space.
        for beta=1:number_of_local_basis_test %loop over trial functions on the element for test FE space.
            int_value=J*Gauss_quadrature_2D_volume_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,alpha,der_x_trial,der_y_trial,basis_type_test,beta,der_x_test,der_y_test,n);
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+int_value;
            AA(alpha,beta)=int_value;
        end
    end
end
