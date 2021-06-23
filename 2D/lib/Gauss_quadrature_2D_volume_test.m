function result=Gauss_quadrature_2D_volume_test(f_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_test,basis_index_test,der_x_test,der_y_test)
Gpn=size(Gauss_nodes,1);% number of Gauss points.
r=0;
for k=1:Gpn %loop on the Gaussian nodes
         %maps reference basis function on the interval
        r=r+Gauss_weights(k)*feval(f_fun,Gauss_nodes(k,1),Gauss_nodes(k,2))*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type_test,basis_index_test,der_x_test,der_y_test);
end
result=r;
end