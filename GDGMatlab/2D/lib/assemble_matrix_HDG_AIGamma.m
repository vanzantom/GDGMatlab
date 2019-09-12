function A=assemble_matrix_HDG_AIGamma(coe_fun,P,T,Tb_trial,Eb_test,Eb,nsub2,ngamma,basis_type_trial,basis_type_test,order_Gauss,alpha_coef)

ii=find(Eb(5,:)==1); %find edges on the interface.
number_of_local_basis_trial=size(Tb_trial,1);
number_of_local_basis_test=size(Eb(6:end,1),1);
A=sparse(nsub2+ngamma,nsub2+ngamma);
for n=1:length(ii)
    vertex1=Eb(1,ii(n));
    vertex2=Eb(2,ii(n));
    vertices1=P(:,vertex1);
    vertices2=P(:,vertex2);
    edge_length=norm(P(:,vertex1)-P(:,vertex2),2);
    mu=alpha_coef/edge_length;
    T_1=Eb(3,ii(n));
    T_2=Eb(4,ii(n)); %don't need to check if T_j are -1 since the interface is in the interior.
    vertices_T_1=P(:,T((1:3),T_1));
    vertices_T_2=P(:,T((1:3),T_2));
    [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(vertices1,vertices2,order_Gauss); % create Gauss nodes on the edge. 1 trapezoidal 2 Simpsons rule.
    n1=find_normal(vertex1,vertex2,T((1:3),T_1),P);
    n2=find_normal(vertex1,vertex2,T((1:3),T_2),P);
    for alpha=1:number_of_local_basis_trial
        for beta=1:number_of_local_basis_test 
            %contribution of -mu\int v_i\phi
            int_value=Gauss_quadrature_1D_trial_test_coupling(coe_fun,Gauss_nodes,Gauss_weights,vertices_T_1,vertices1,vertices2,basis_type_trial,alpha,0,0,basis_type_test,beta,0);
            A(Eb_test(5+beta,ii(n)),Tb_trial(alpha,T_1))=A(Eb_test(5+beta,ii(n)),Tb_trial(alpha,T_1))-mu*int_value;
            int_value=Gauss_quadrature_1D_trial_test_coupling(coe_fun,Gauss_nodes,Gauss_weights,vertices_T_2,vertices1,vertices2,basis_type_trial,alpha,0,0,basis_type_test,beta,0);
            A(Eb_test(5+beta,ii(n)),Tb_trial(alpha,T_2))=A(Eb_test(5+beta,ii(n)),Tb_trial(alpha,T_2))-mu*int_value;
            %contribution of \int \partial_n v_i\phi
            int_value=n1(1)*Gauss_quadrature_1D_trial_test_coupling(coe_fun,Gauss_nodes,Gauss_weights,vertices_T_1,vertices1,vertices2,basis_type_trial,alpha,1,0,basis_type_test,beta,0);
            int_value=int_value+n1(2)*Gauss_quadrature_1D_trial_test_coupling(coe_fun,Gauss_nodes,Gauss_weights,vertices_T_1,vertices1,vertices2,basis_type_trial,alpha,0,1,basis_type_test,beta,0);
            A(Eb_test(5+beta,ii(n)),Tb_trial(alpha,T_1))=A(Eb_test(5+beta,ii(n)),Tb_trial(alpha,T_1))+int_value;
            int_value=n2(1)*Gauss_quadrature_1D_trial_test_coupling(coe_fun,Gauss_nodes,Gauss_weights,vertices_T_2,vertices1,vertices2,basis_type_trial,alpha,1,0,basis_type_test,beta,0);
            int_value=int_value+n2(2)*Gauss_quadrature_1D_trial_test_coupling(coe_fun,Gauss_nodes,Gauss_weights,vertices_T_2,vertices1,vertices2,basis_type_trial,alpha,0,1,basis_type_test,beta,0);
            A(Eb_test(5+beta,ii(n)),Tb_trial(alpha,T_2))=A(Eb_test(5+beta,ii(n)),Tb_trial(alpha,T_2))+int_value;
        end
    end
  
 end

 end