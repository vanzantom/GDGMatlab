function int_value=Gauss_quadrature_2D_volume_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,basis_index_trial,der_x_trial,der_y_trial,basis_type_test,basis_index_test,der_x_test,der_y_test,n)
%Gauss_quadrature_2D_volume_trial_test computes the integral using Gaussian
%quadrature. It receives the gauss nodes and weigths, the vertices of the element, the basis_type, 
%the basis index and the derivative order for the FEM space both for trial
%and test FE spaces.

 Gpn=size(Gauss_nodes,1);% number of Gauss points.


 int_value=0;%initialize first value
% r=0;
% for k=1:Gpn %loop on the Gaussian nodes
%         %maps reference basis function on the reference element
%         r=r+Gauss_weights(k)*feval(coe_fun,Gauss_nodes(k,1),Gauss_nodes(k,2),n)*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type_trial,basis_index_trial,der_x_trial,der_y_trial)*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type_test,basis_index_test,der_x_test,der_y_test);
% end
% int_value=r;
c=zeros(Gpn,1);
for k=1:Gpn
    c(k)=feval(coe_fun,Gauss_nodes(k,1),Gauss_nodes(k,2),n);
end

int_value=Gauss_weights*(c.*FE_local_basis_2D(Gauss_nodes(:,1),Gauss_nodes(:,2),vertices,basis_type_trial,basis_index_trial,der_x_trial,der_y_trial).*FE_local_basis_2D(Gauss_nodes(:,1),Gauss_nodes(:,2),vertices,basis_type_test,basis_index_test,der_x_test,der_y_test));

end