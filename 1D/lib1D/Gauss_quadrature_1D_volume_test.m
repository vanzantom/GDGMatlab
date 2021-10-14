function result=Gauss_quadrature_1D_volume_test(f_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_test,basis_index_test,der_test)
h=vertices(end)-vertices(1);
Gpn=length(Gauss_nodes);% number of Gauss points.
r=0;
for k=1:Gpn %loop on the Gaussian nodes
        r=r+Gauss_weights(k)*feval(f_fun,Gauss_nodes(k))*FE_local_basis(Gauss_nodes(k),vertices,basis_type_test,basis_index_test,der_test);
        end
result=r;
end