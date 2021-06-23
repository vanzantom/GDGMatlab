function int_value=Gauss_quadrature_1D_interface_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices1,vertices2,basis_type_trial,basis_index_trial,der_trial,basis_type_test,basis_index_test,der_test)
%Gauss_quadrature_1D_interface_trial_test computes the integral using Gaussian
%quadrature. It receives the gauss nodes and weigths, the vertices of the element, the basis_type, 
%the basis index and the derivative order for the DG space.
int_value=0;%initialize first value
Gpn=length(Gauss_nodes);% number of Gauss points.
r=0;
for k=1:Gpn %loop on the Gaussian nodes
        %use explicit expression of local basis function
        r=r+Gauss_weights(k)*feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k))*FE_local_basis_1D(Gauss_nodes(:,k),vertices1,vertices2,basis_type_trial,basis_index_trial,der_trial)*FE_local_basis_1D(Gauss_nodes(:,k),vertices1,vertices2,basis_type_test,basis_index_test,der_test);
        
end
int_value=r;
end