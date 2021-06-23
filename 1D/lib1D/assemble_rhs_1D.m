%-------------------------------------------------------------------------
% assemble_rhs_1D receives 
% f_fun: force term
% P,T: mesh matrices ,
% Tb_test: matrices with information on trial and test FEM space
% matrixsize2: size of test FEM space
% basis_type_test: finite element space type test
% der_test: order of derivatives test function
% assemble_rhs_1D returns: 
% b (matrix) whose entry b_{i}= \int (f_fun)((\partial_x)^der_test \phi_i) 


% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function b=assemble_rhs_1D(f_fun,P,T,Tb_test,matrixsize2,basis_type_test,der_test)
b=zeros(matrixsize2,1); % create output vector
number_of_elements=size(T,2); %number of elements
number_of_local_basis_test=size(Tb_test,1);%number of local FE basis on each element
for n=1:number_of_elements % loop over elements
    vertices=P(:,T(:,n));% vertices of the element
    [Gauss_nodes,Gauss_weights]=generate_Gauss(vertices,2);
    for beta=1:number_of_local_basis_test %sum over contributions
        b(Tb_test(beta,n),1)=b(Tb_test(beta,n),1)+Gauss_quadrature_1D_volume_test(f_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_test,beta,der_test);
    end
end