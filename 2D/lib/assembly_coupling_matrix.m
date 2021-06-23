%The function assembly_coupling_matrix computes the integral_edge_label
%phi_iphi_j where phi_i and phi_j are test and trial function possible
%living on different meshes.
%It receives
%Pbtest,Tb_test,Pbtrial,Tb_trial:data matrices for FE space trial and test
%boundarytest/trial: data structure for the boundary of the two meshes
%label_edge: label of the edge over which we compute the integral
%oaratrial/paratest: parameters for FE space trial and test
%x0,x1,y0: extrema of the segment edge_label. Needed if we aim to impose Dirichlet BC at the two vertex (x0,y0) and (x1,y1). 
%It returns:
%C : coupling matrix
%C1,C2: modified coupling matrix which do not include the extrema of the
%points as DOF for the test space.


function [C,C1,C2]=assembly_coupling_matrix(Pbtest,Tb_test,Pbtrial,Tb_trial,boundarytest,boundarytrial,label_edge,paratest,paratrial,x1,x2,y0)
    
    ndofve=size(Pbtest,2);
    ndofq=size(Pbtrial,2);
    C=sparse(ndofve,ndofq);
   
    
    indexedgeS=find(boundarytest(1,:)==label_edge);
    indexedgeD=find(boundarytrial(1,:)==label_edge);
    number_of_local_basis_test=size(Tb_test,1);
    number_of_local_basis_trial=size(Tb_trial,1);
    for l=1:length(indexedgeS)
        k=indexedgeS(l);
        El=boundarytest(2,k);% Element of Stokes which has this boundary
         %number of local basis function on element El
        Vertex_El=Pbtest(:,Tb_test(1:3,El));%vertices of the element
        vertices=Pbtest(:,[boundarytest(3,k);boundarytest(4,k)]);%vertices of the edge
        j=1;
        flag=0;
        while(j<=length(indexedgeD) && flag==0)%find index of the current edge k=indexedgeS(l) in the data structure boundarytrial
            ii=indexedgeD(j);
            vertices_Darcy=Pbtrial(:,[boundarytrial(3,ii);boundarytrial(4,ii)]);
            if( norm(vertices-vertices_Darcy)<10^-8 || norm(fliplr(vertices)-vertices_Darcy)<10^-8)
                 flag=1;
                 j=j-1;
                 end
            j=j+1;
        end
        El_Darcy=boundarytrial(2,ii);% Element of Darcy which has this boundary  
        Vertex_ElDarcy=Pbtrial(:,Tb_trial(1:3,El_Darcy));
        [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(vertices(:,1),vertices(:,2),paratest.order_Gauss);%already multiplied by length edge
                for alpha=1:number_of_local_basis_trial %loop over trial functions on the element for trial FE space.
                    for beta=1:number_of_local_basis_test %loop over trial functions on the element for test FE space.
                        int_value=Gauss_quadrature_1D_trial_test(@(x,y) 1,Gauss_nodes,Gauss_weights,Vertex_ElDarcy,Vertex_El,paratrial.basis_type,alpha,0,0,paratest.basis_type_ve,beta,0,0);
                        C(Tb_test(beta,El),Tb_trial(alpha,El_Darcy))=C(Tb_test(beta,El),Tb_trial(alpha,El_Darcy))+int_value;
                    end
                end
    end
    %Attention! I built  C correctly, but it counts also the degree of
    %freedom which is on the angle where I imposed Dirichlet BC!
    l=1;
    for ii=1:size(Pbtest,2)
         if( (abs(Pbtest(1,ii)-x1)<10^-8  ||abs(Pbtest(1,ii)-x2)<10^-8) &&  abs(Pbtest(2,ii)-y0)<10^-8)
         vertexS(l)=ii;
         l=l+1;
        end
    
    end
    l=1;
    for ii=1:size(Pbtrial,2)
         if( (abs(Pbtrial(1,ii)-x1)<10^-8  ||abs(Pbtrial(1,ii)-x2)<10^-8) &&  abs(Pbtrial(2,ii)-y0)<10^-8)
         vertexD(l)=ii;
         l=l+1;
        end
    
    end
    C1=-C;
    C2=C';
    C1(vertexS,:)=0;
    C2(vertexD,:)=0;
    
end