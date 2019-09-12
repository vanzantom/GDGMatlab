function A=assemble_matrix_HDG_AGamma(coe_fun,P,Eb_trial,Eb_test,Eb,nsub2,ngamma,basis_type_trial,basis_type_test,order_Gauss,alpha_coef)

ii=find(Eb(5,:)==1); %find edges on the interface.
number_of_local_basis_trial=size(Eb(6:end,1),1);
number_of_local_basis_test=size(Eb(6:end,1),1);
A=sparse(nsub2+ngamma,nsub2+ngamma);
for n=1:length(ii)
    vertex1=Eb(1,ii(n));
    vertex2=Eb(2,ii(n));
    vertices1=P(:,vertex1);
    vertices2=P(:,vertex2);
    edge_length=norm(P(:,vertex1)-P(:,vertex2),2);
    mu=alpha_coef/edge_length;
    [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),order_Gauss); % create Gauss nodes on the edge. 1 trapezoidal 2 Simpsons rule.
    for alpha=1:number_of_local_basis_trial
        for beta=1:number_of_local_basis_test %from reference-> local -> global assembly.
            int_value=Gauss_quadrature_1D_interface_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices1,vertices2,basis_type_trial,alpha,0,basis_type_test,beta,0);
            A(Eb_test(5+beta,ii(n)),Eb_trial(5+alpha,ii(n)))=A(Eb_test(5+beta,ii(n)),Eb_trial(5+alpha,ii(n)))+mu*int_value;
        end
    end
  
 end

 end