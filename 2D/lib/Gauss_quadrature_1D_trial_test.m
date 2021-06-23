function int_value=Gauss_quadrature_1D_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices1,vertices2,basis_type_trial,basis_index_trial,der_x_trial,der_y_trial,basis_type_test,basis_index_test,der_x_test,der_y_test)
%Gauss_quadrature_1D_trial_test computes the integral using Gaussian
%quadrature. It receives the gauss nodes and weigths, the vertices of the element, the basis_type, 
%the basis index and the derivative order for the FEM space both for trial
%and test FE spaces.
% Vertices 1 refers to the vertices of the trial basis, Vertices 2
% correspond to the vertices of the test triangle.

int_value=0;%initialize first value
Gpn=size(Gauss_nodes,2);% number of Gauss points.
r=0;
for k=1:Gpn %loop on the Gaussian nodes
        %maps reference basis function on the reference element
        r=r+Gauss_weights(k)*feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k))*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices1,basis_type_trial,basis_index_trial,der_x_trial,der_y_trial)*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices2,basis_type_test,basis_index_test,der_x_test,der_y_test);
end
int_value=r;
end