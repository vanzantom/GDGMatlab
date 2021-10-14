%--------------------------------------------------------------------------
% the function stabilization() helps to compute the terms \mu\int_e [u][v] present in the IP formulation (See blue
% notebook for the calculations) At each interface I have 4 terms like
% u(x+)*v(x+) or u(x+)*v(x-) 
% It receveis verticesj to move from local to a r framework. The
% variable vertex contains the position of the interface.
% CAREFUL: vertices1 refers to the vertices of the element for the trial
% basis. Vertices 2 intestead refers to vertices for the test's element

% author: Tommaso Vanzan
%--------------------------------------------------------------------------
function result=stabilization(penalty,vertices1,vertices2,vertex,basis_type_trial,basis_index_trial,basis_type_test,basis_index_test)

h=vertices1(end)-vertices1(1);
result=(penalty)*FE_local_basis(vertex,vertices1,basis_type_trial,basis_index_trial,0)*FE_local_basis(vertex,vertices2,basis_type_test,basis_index_test,0);
end

