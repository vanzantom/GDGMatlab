function A=assemble_matrix_edges1DIP(coe_fun,P,T,Tb_trial,Tb_test,Eb,matrixsize1,matrixsize2,basis_type_trial,basis_type_test,penalty)
% The function assemble_matrix_edgeds1D loops over the edges and computes
% the edge terms of the SIP. See blue notebook.

A=sparse(matrixsize1,matrixsize2);%create sparse matrix
number_of_elements=size(T,2); %number of Elements
number_of_local_basis_trial=size(Tb_trial,1); %number of trial local basis function on a single element
number_of_local_basis_test=size(Tb_test,1); %number of test local basis function on a single element
number_of_edges=size(Eb,2); %number of Elements
  for n=1:number_of_edges   %cycle over the interior edges.
      T1=Eb(1,n);% get index triangles which share the edge n
      T2=Eb(2,n);
      vertices1=P(:,T(:,T1));
      vertices2=P(:,T(:,T2));
      vertex=vertices1(end);% the edge is located at the end of T1
    for alpha=1:number_of_local_basis_trial % Contribution from term involving basis on same element.
        for beta=1:number_of_local_basis_test %from reference-> local -> global assembly.
            %penalization \mu block
            A(Tb_test(beta,T1),Tb_trial(alpha,T1))=A(Tb_test(beta,T1),Tb_trial(alpha,T1))+stabilization(penalty,vertices1,vertices1,vertex,basis_type_trial,beta,basis_type_test,alpha);
            A(Tb_test(beta,T2),Tb_trial(alpha,T1))=A(Tb_test(beta,T2),Tb_trial(alpha,T1))-stabilization(penalty,vertices2,vertices1,vertex,basis_type_trial,beta,basis_type_test,alpha);
            A(Tb_test(beta,T1),Tb_trial(alpha,T2))=A(Tb_test(beta,T1),Tb_trial(alpha,T2))-stabilization(penalty,vertices1,vertices2,vertex,basis_type_trial,beta,basis_type_test,alpha);
            A(Tb_test(beta,T2),Tb_trial(alpha,T2))=A(Tb_test(beta,T2),Tb_trial(alpha,T2))+stabilization(penalty,vertices2,vertices2,vertex,basis_type_trial,beta,basis_type_test,alpha);
            % symmetric blocks
            A(Tb_test(beta,T1),Tb_trial(alpha,T1))=A(Tb_test(beta,T1),Tb_trial(alpha,T1))-1/2*fluxconsistency(coe_fun,vertices1,vertices1,vertex,basis_type_trial,beta,basis_type_test,alpha,1,0)-1/2*fluxconsistency(coe_fun,vertices1,vertices1,vertex,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,T2),Tb_trial(alpha,T1))=A(Tb_test(beta,T2),Tb_trial(alpha,T1))-1/2*fluxconsistency(coe_fun,vertices2,vertices1,vertex,basis_type_trial,beta,basis_type_test,alpha,1,0)+1/2*fluxconsistency(coe_fun,vertices2,vertices1,vertex,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,T1),Tb_trial(alpha,T2))=A(Tb_test(beta,T1),Tb_trial(alpha,T2))+1/2*fluxconsistency(coe_fun,vertices1,vertices2,vertex,basis_type_trial,beta,basis_type_test,alpha,1,0)-1/2*fluxconsistency(coe_fun,vertices1,vertices2,vertex,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,T2),Tb_trial(alpha,T2))=A(Tb_test(beta,T2),Tb_trial(alpha,T2))+1/2*fluxconsistency(coe_fun,vertices2,vertices2,vertex,basis_type_trial,beta,basis_type_test,alpha,1,0)+1/2*fluxconsistency(coe_fun,vertices2,vertices2,vertex,basis_type_trial,beta,basis_type_test,alpha,0,1);
        end
    end
  end
% I treat the boundary terms in penalization way.
    T1=1; % First triangle.
    vertices1=P(:,T(:,T1));
    vertex=vertices1(1);% left boundary
    T2=number_of_elements;
    vertices2=P(:,T(:,T2));
    vertex2=vertices2(end);%right boundary
    for alpha=1:number_of_local_basis_trial % Contribution from term involving basis on same element.
        for beta=1:number_of_local_basis_test %from reference-> local -> global assembly.
            %penalization \mu block
            A(Tb_test(beta,T1),Tb_trial(alpha,T1))=A(Tb_test(beta,T1),Tb_trial(alpha,T1))+stabilization(penalty,vertices1,vertices1,vertex,basis_type_trial,beta,basis_type_test,alpha);
            A(Tb_test(beta,T2),Tb_trial(alpha,T2))=A(Tb_test(beta,T2),Tb_trial(alpha,T2))+stabilization(penalty,vertices2,vertices2,vertex2,basis_type_trial,beta,basis_type_test,alpha);
            %symmetry block
            A(Tb_test(beta,T1),Tb_trial(alpha,T1))=A(Tb_test(beta,T1),Tb_trial(alpha,T1))+fluxconsistency(coe_fun,vertices1,vertices1,vertex,basis_type_trial,beta,basis_type_test,alpha,1,0)+fluxconsistency(coe_fun,vertices1,vertices1,vertex,basis_type_trial,beta,basis_type_test,alpha,0,1);
            A(Tb_test(beta,T2),Tb_trial(alpha,T2))=A(Tb_test(beta,T2),Tb_trial(alpha,T2))-fluxconsistency(coe_fun,vertices2,vertices2,vertex2,basis_type_trial,beta,basis_type_test,alpha,1,0)-fluxconsistency(coe_fun,vertices2,vertices2,vertex2,basis_type_trial,beta,basis_type_test,alpha,0,1);
        end
    end
end
