%-------------------------------------------------------------------------
% assemble_matrix_1D_IP_element assembles the extra `DG' terms of an
% interior penalty discretization looping over elements, see Prof. Antonietti slides
% assemble_matrix_1D_IP_element receives
% coe_fun: the diffusion coefficient 
% mesh matrices: T,P 
% Tb_trial,Tb_test: matrices with information on trial and test finite element space 
% matrixsize1,matrixsize2: size of the stiffness matrix 
% basis_type_trial,basis_type_test: finite element space type trial/test
% der_trial,der_test: order of derivatives in the bilinea form .
% assemble_matrix_1D returns: 
% A (matrix) whose entry A_{i,j}= \int ((\partial_x)^der_trial \phi_j)((\partial_x)^der_test \phi_i) 


% author: Tommaso Vanzan
%-------------------------------------------------------------------------
function A=assemble_matrix_1D_IP_element(coe_fun,P,T,Tb_trial,Tb_test,matrixsize1,matrixsize2,basis_type_trial,basis_type_test,penalty)
A=sparse(matrixsize1,matrixsize2);%create sparse matrix
number_of_elements=size(T,2); %number of Elements
number_of_local_basis_trial=size(Tb_trial,1); %number of trial local basis function on a single element
number_of_local_basis_test=size(Tb_test,1); %number of test local basis function on a single element
for n=1:number_of_elements %loop over elements
    vertices=P(:,T(:,n)); % mesh information then I use P and T
    for alpha=1:number_of_local_basis_trial
        for beta=1:number_of_local_basis_test %from reference-> local -> global assembly.
            if n>1 && n<number_of_elements
            vertices0=P(:,T(:,n-1)); %vertices element on the left
            vertices2=P(:,T(:,n+1)); %vertices element on the right
            vertex1=vertices(1); %vertex interface element n on the left
            vertex2=vertices(end); %vertex interface element n on the right
            %stabilization terms
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+stabilization(penalty,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha)+stabilization(penalty,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))-stabilization(penalty,vertices2,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))-stabilization(penalty,vertices0,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha);
            %flux consistency terms int_{E1 U E2} {{\nabla\psi_i}}[\psi_j]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1)+1/2*fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,0,1);
            % second integral int_{E1 U E2} {{\nabla\psi_j}}[\psi_i]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0)+1/2*fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,1,0);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,1,0);
            elseif n==1 % I treat separately first and last subdomain
            vertices2=P(:,T(:,n+1)); %vertices element on the right
            vertex1=vertices(1); %vertex interface element n on the left
            vertex2=vertices(end); %vertex interface element n on the right.  I first compute the terms involing the 'interior' edge located at vertex2
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+stabilization(penalty,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))-stabilization(penalty,vertices2,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha);
            %flux consistency terms int_{E1 U E2} {{\nabla\psi_i}}[\psi_j]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1);
            % second integral int_{E1 U E2} {{\nabla\psi_j}}[\psi_i]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0);
            % Boundary conditions. Now i add the term located at vertex1.
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+stabilization(penalty,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha); %stabilization term
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,0,1)+fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,1,0);
            elseif n==number_of_elements
            vertices0=P(:,T(:,n-1)); %vertices element on the left
            vertex1=vertices(1); %vertex interface element n on the left
            vertex2=vertices(end); %vertex interface element n on the right I first compute the terms involing the 'interior' edge located at vertex1
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+stabilization(penalty,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))-stabilization(penalty,vertices0,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha);
            %flux consistency terms int_{E1 U E2} {{\nabla\psi_i}}[\psi_j]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,0,1);
            % second integral int_{E1 U E2} {{\nabla\psi_j}}[\psi_i]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,1,0);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,1,0);
            %Boundary condition 
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+stabilization(penalty,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha);
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1)-fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0);
            end
        end
    end
end
