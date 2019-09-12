function A=assemble_matrix_1D_DG_element(coe_fun,P,T,Tb_trial,Tb_test,matrixsize1,matrixsize2,basis_type_trial,basis_type_test,penalty)
% The function assemble_matrix_1D_DG_element loops over the elements to
% build the extra term of the SIP method. See Antonietti slides and
% calculations on the blue notebook. I do not loop over the edges but over
% the elements!
A=sparse(matrixsize1,matrixsize2);%create sparse matrix
number_of_elements=size(T,2); %number of Elements
number_of_local_basis_trial=size(Tb_trial,1); %number of trial local basis function on a single element
number_of_local_basis_test=size(Tb_test,1); %number of test local basis function on a single element
for n=1:number_of_elements
    vertices=P(:,T(:,n)); % mesh information then I use P and T
    for alpha=1:number_of_local_basis_trial
        for beta=1:number_of_local_basis_test %from reference-> local -> global assembly.
            if n>1 && n<number_of_elements
            vertices0=P(:,T(:,n-1)); %vertices element on the left
            vertices2=P(:,T(:,n+1)); %vertices element on the right
            vertex1=vertices(1); %vertex interface element n on the left
            vertex2=vertices(end); %vertex interface element n on the right
            %stabilization terms
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+penalty/2*stabilizationIPH(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,0,basis_type_test,alpha,0)+penalty/2*stabilizationIPH(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,0,basis_type_test,alpha,0);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))-penalty/2*stabilizationIPH(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,0,basis_type_test,alpha,0);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))-penalty/2*stabilizationIPH(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,0,basis_type_test,alpha,0);
            %flux consistency terms int_{E1 U E2} {{\nabla\psi_i}}[\psi_j]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1)+1/2*fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,0,1);
            % second integral int_{E1 U E2} {{\nabla\psi_j}}[\psi_i]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0)+1/2*fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,1,0);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,1,0);
            % Additional term IPH
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/(2*penalty)*stabilizationIPH(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,1,basis_type_test,alpha,1)-1/(2*penalty)*stabilizationIPH(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,1,basis_type_test,alpha,1);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))+1/(2*penalty)*stabilizationIPH(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,1,basis_type_test,alpha,1);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))+1/(2*penalty)*stabilizationIPH(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,1,basis_type_test,alpha,1);
            elseif n==1 % I treat separately first and last subdomain
            vertices2=P(:,T(:,n+1)); %vertices element on the right
            vertex1=vertices(1); %vertex interface element n on the left
            vertex2=vertices(end); %vertex interface element n on the right.  I first compute the terms involing the 'interior' edge located at vertex2
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+penalty/2*stabilizationIPH(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,0,basis_type_test,alpha,0);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))-penalty/2*stabilizationIPH(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,0,basis_type_test,alpha,0);
            %flux consistency terms int_{E1 U E2} {{\nabla\psi_i}}[\psi_j]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1);
            % second integral int_{E1 U E2} {{\nabla\psi_j}}[\psi_i]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0);
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0);
            % Additional IPH term
            A(Tb_test(beta,n+1),Tb_trial(alpha,n))=A(Tb_test(beta,n+1),Tb_trial(alpha,n))+1/(2*penalty)*stabilizationIPH(coe_fun,vertices2,vertices,vertex2,basis_type_trial,beta,1,basis_type_test,alpha,1);
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/(2*penalty)*stabilizationIPH(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,1,basis_type_test,alpha,1);
            % Boundary conditions. Now I add the term located at vertex1.
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+penalty/2*stabilizationIPH(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,0,basis_type_test,alpha,0); %stabilization term
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,0,1)+fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,1,0);
            elseif n==number_of_elements
            vertices0=P(:,T(:,n-1)); %vertices element on the left
            vertex1=vertices(1); %vertex interface element n on the left
            vertex2=vertices(end); %vertex interface element n on the right I first compute the terms involing the 'interior' edge located at vertex1
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+penalty/2*stabilizationIPH(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,0,basis_type_test,alpha,0);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))-penalty/2*stabilizationIPH(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,0,basis_type_test,alpha,0);
            %flux consistency terms int_{E1 U E2} {{\nabla\psi_i}}[\psi_j]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))-1/2*fluxconsistency(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,0,1);
            % second integral int_{E1 U E2} {{\nabla\psi_j}}[\psi_i]
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,1,0);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))+1/2*fluxconsistency(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,basis_type_test,alpha,1,0);
            %additional IPH
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-1/(2*penalty)*stabilizationIPH(coe_fun,vertices,vertices,vertex1,basis_type_trial,beta,1,basis_type_test,alpha,1);
            A(Tb_test(beta,n-1),Tb_trial(alpha,n))=A(Tb_test(beta,n-1),Tb_trial(alpha,n))+1/(2*penalty)*stabilizationIPH(coe_fun,vertices0,vertices,vertex1,basis_type_trial,beta,1,basis_type_test,alpha,1);
            %Boundary condition 
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+penalty/2*stabilizationIPH(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,0,basis_type_test,alpha,0);
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))-fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1)-fluxconsistency(coe_fun,vertices,vertices,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0);
            end
        end
    end
end
