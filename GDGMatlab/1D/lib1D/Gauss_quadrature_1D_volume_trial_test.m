function int_value=Gauss_quadrature_1D_volume_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,basis_index_trial,der_trial,basis_type_test,basis_index_test,der_test)
%Gauss_quadrature_1D_volume_trial_test computes the integral using Gaussian
%quadrature. It receives the gauss nodes and weigths, the vertices of the element, the basis_type, 
%the basis index and the derivative order for the FEM space.
h=vertices(end)-vertices(1);
int_value=0;%initialize first value
Gpn=length(Gauss_nodes);% number of Gauss points.
r=0;
for k=1:Gpn %loop on the Gaussian nodes
        %use explicit expression of local basis function
        %r=r+Gauss_weights(k)*feval(coe_fun,Gauss_nodes(k))*FE_local_basis(Gauss_nodes(k),vertices,basis_type_trial,basis_index_trial,der_trial)*FE_local_basis(Gauss_nodes(k),vertices,basis_type_test,basis_index_test,der_test);
        %maps reference basis function on the interval
        r=r+1/h^2*Gauss_weights(k)*feval(coe_fun,Gauss_nodes(k))*FE_reference_basis_1D((Gauss_nodes(k)-vertices(1))/h,basis_type_trial,basis_index_trial,der_trial)*FE_reference_basis_1D((Gauss_nodes(k)-vertices(1))/h,basis_type_test,basis_index_test,der_test);
end
int_value=r;
end