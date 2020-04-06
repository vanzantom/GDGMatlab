%-------------------------------------------------------------------------
% assemble_rhs_2D receives 
% f_fun: force term
% P,T: mesh matrices ,
% Tb_test: matrices with information on trial and test FEM space
% matrixsize2: size of test FEM space
% basis_type_test: finite element space type test
% der_test_x,der_test_y: order of derivatives test function
% assemble_rhs_1D returns: 
% b (matrix) whose entry b_{i}= \int (f_fun)((\partial_x)^der_test_x (\partial_x)^der_test_y \phi_i) 


% author: Tommaso Vanzan
%-------------------------------------------------------------------------


function b=assemble_rhs_2D(f_fun,P,T,Tb_test,matrixsize2,basis_type_test,der_x_test,der_y_test)
b=zeros(matrixsize2,1);
number_of_elements=size(T,2);
number_of_local_basis_test=size(Tb_test,1);
for n=1:number_of_elements %loop over elements
    vertices=P(:,T((1:3),n));%vertices element
    [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,2);%generate Gauss quadrature
    J=abs((vertices(1,2)-vertices(1,1))*(vertices(2,3)-vertices(2,1)) -(vertices(1,3)-vertices(1,1))*(vertices(2,2)-vertices(2,1)));%Jacobian affine transformation
    for beta=1:number_of_local_basis_test %sum all contributions
        b(Tb_test(beta,n),1)=b(Tb_test(beta,n),1)+J*Gauss_quadrature_2D_volume_test(f_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_test,beta,der_x_test,der_y_test);
    end
end