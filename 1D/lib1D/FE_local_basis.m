%-------------------------------------------------------------------------
% FE_local_basis receives a point in x in which evaluete the basis function
% on the local element.
% basis_type: specifies  the kind of FE Space.
% basis_index: is a variable containing the index of the local basis
% function inside the local element.(Ex: 1D linear I have two basis functions.)
% der: specifies if I compute the value of the function or of its
% derivative.

%author: Tommaso Vanzan
%-------------------------------------------------------------------------
function result=FE_local_basis(x,vertices,basis_type,basis_index,der)
h=vertices(2)-vertices(1);
x_hat=(x-vertices(1))/h;
if der==0
    result=FE_reference_basis_1D(x_hat,basis_type,basis_index,der);
elseif der==1
    result=(1/h)*FE_reference_basis_1D(x_hat,basis_type,basis_index,der);
else
    fprintf('Order derivative higher than expected');
end