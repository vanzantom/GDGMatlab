%-------------------------------------------------------------------------
%Gauss_quadrature_1D_volume_trial_test computes the integral 
% int coef_fun (( \partial_x)^der_trial \phi_trial ) (( \partial_x)^der_test \phi_test )
%using Gaussian quadrature. It receives the gauss nodes and weigths, the vertices of the element, the basis_type, 
%the basis index and the derivative order for the FEM space.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------


function int_value=Gauss_quadrature_1D_volume_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,basis_index_trial,der_trial,basis_type_test,basis_index_test,der_test)
h=vertices(end)-vertices(1);
int_value=0;%initialize first value
Gpn=length(Gauss_nodes);% number of Gauss points.
r=0;
for k=1:Gpn %loop on the Gaussian nodes
        r=r+Gauss_weights(k)*feval(coe_fun,Gauss_nodes(k))*FE_local_basis(Gauss_nodes(k),vertices,basis_type_trial,basis_index_trial,der_trial)*FE_local_basis(Gauss_nodes(k),vertices,basis_type_test,basis_index_test,der_test);
end
int_value=r;
end