function b=assemble_rhs_2D(f_fun,P,matrixsize2,T,Tb_test,basis_type_test,der_x_test,der_y_test)
% assemble_rhs_2D computes the RHS integrating f against all the test
% functions of the FE space designed by 'basis_type_test').
b=zeros(matrixsize2,1);
number_of_elements=size(T,2);
number_of_local_basis_test=size(Tb_test,1);
for n=1:number_of_elements
    vertices=P(:,T((1:3),n));
    [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,2);
    J=abs((vertices(1,2)-vertices(1,1))*(vertices(2,3)-vertices(2,1)) -(vertices(1,3)-vertices(1,1))*(vertices(2,2)-vertices(2,1)));
    for beta=1:number_of_local_basis_test
        b(Tb_test(beta,n),1)=b(Tb_test(beta,n),1)+J*Gauss_quadrature_2D_volume_test(f_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_test,beta,der_x_test,der_y_test);
    end
end