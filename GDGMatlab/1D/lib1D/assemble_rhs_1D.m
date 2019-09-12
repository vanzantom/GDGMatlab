function b=assemble_rhs_1D(f_fun,P,matrixsize2,T,Tb_test,basis_type_test,der_test)
% Function assemble_rhs_1D assembles the right hand side of the FEM
% problem.
b=zeros(matrixsize2,1);
number_of_elements=size(T,2);
number_of_local_basis_test=size(Tb_test,1);
for n=1:number_of_elements
    vertices=P(:,T(:,n));
    [Gauss_nodes,Gauss_weights]=generate_Gauss(vertices,2);
    for beta=1:number_of_local_basis_test
        b(Tb_test(beta,n),1)=b(Tb_test(beta,n),1)+Gauss_quadrature_1D_volume_test(f_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_test,beta,der_test);
    end
end