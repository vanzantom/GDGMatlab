%-------------------------------------------------------------------------
% assemble_matrix_HDG_AGamma receives 
% coe_fun: the diffusion coefficient 
% mesh matrix:  P
% Eb_trial,Eb_test: matrices with information on trial and test finite
% element space on Gamma
% Nsub: is a vector. Nsub(jj) last index of DOF in subdomain jj.
% ngamma: number DOFs on interface Gamma.
% basis_type_trial,basis_type_test: finite element space type trial/test
% der_x_trial,der_x_test, der_y_trial,der_y_test,: order of derivatives in the bilinea form
% order_gauss: Order of Gauss quadrature.
% alpha_coef: DG penalization parameter
% assemble_matrix_HDG_AGamma returns: 
% A_Gamma (matrix) whose entry A_{i,j}= \alpha/h \int_\Gamma  \phi_j \phi_i)


% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function A=assemble_matrix_HDG_AGamma(coe_fun,P,Eb_trial,Eb_test,Eb,Nsub,ngamma,basis_type_trial,basis_type_test,para)

order_Gauss=para.order_Gauss;
alpha_coef=para.alpha_coef;
ii=find(Eb(5,:)==1); %find edges on the interface.
number_of_local_basis_trial=size(Eb(6:end,1),1);
number_of_local_basis_test=size(Eb(6:end,1),1);
A=sparse(Nsub(end)+ngamma,Nsub(end)+ngamma);
for n=1:length(ii) %for each edge on interface Gamma, compute contributions.
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