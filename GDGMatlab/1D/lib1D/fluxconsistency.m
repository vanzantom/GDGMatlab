function result=fluxconsistency(coe_fun,vertices1,vertices2,vertex,basis_type_trial,basis_index_trial,basis_type_test,basis_index_test,der_trial,der_test)
% the function fluxconsistency() helps to compute the terms \int_e {u}[v] and \int_e [u]{v}. See blue
% notebook for the calculations. At each interface I have 4 terms like
% -1/2 nu(x)p'(x-)v(x-) -1/2 nu(x)p(x-)v'(x-)
% It receveis verticesj to move from a refence to a local framework. The
% variable vertex contains the position of the interface.
h=vertices1(end)-vertices1(1);
result=1/h*coe_fun(vertex)*FE_reference_basis_1D((vertex-vertices1(1))/h,basis_type_trial,basis_index_trial,der_trial)*FE_reference_basis_1D((vertex-vertices2(1))/h,basis_type_test,basis_index_test,der_test);
end

