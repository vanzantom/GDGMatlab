function [error_L2,error_H1,error_DG]=compute_errors(u,uex,P,T,Tb,Eb,basis_type_trial,order_Gauss,alpha_coef,method)
% The function compute_errors computes the L_2 error and H1 error of the finite element solution
% It receives the solution u, the exact vector uexvec and the mesh info in
% P and T, Tb, Eb, the basis type and order of gauss quadrature.
% If method=='DG' then it computes also the contribution on edges.
% BE CAREFUL!: Error_Edges works only if uex=0 on boundary. Because
% [[uex]]=0 on the interior edges but not on the boundary unless uex=0.
error_L2=0;
error_H1=0;
number_of_elements=size(T,2); %number of Elements
number_of_edges=size(Eb,2);
number_of_local_basis_trial=size(Tb,1);
for n=1:number_of_elements % Loop over the triangles.
    vertices=P(:,T((1:3),n)); % vertices of current Element. Element n, vertices T(:,n), Position P(:,T(:,n)).
    J=abs((vertices(1,2)-vertices(1,1))*(vertices(2,3)-vertices(2,1)) -(vertices(1,3)-vertices(1,1))*(vertices(2,2)-vertices(2,1)));% Compute Jacobian of the linear transformation
    [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,order_Gauss); % create Gauss nodes on the element. 1 baricentri 2 order accuracy.
    u_app_n=u(Tb((1:3),n)); % get indeces of DOF inside the element.
    u_exact_n=uex(Tb((1:3),n)); % get indeces of DOF inside the element.
    for k=1:size(Gauss_nodes,1) % for each gauss node
        err_element=0;% construction error inside each element. \sum FE basis(uapp-uex)
        err_element_x=0; % construction derivative_x error inside each element
        err_element_y=0;% construction derivative_y error inside each element
        for alpha=1:number_of_local_basis_trial %loop over trial functions on the element for trial FE space.
            err_element=err_element +(u_app_n(alpha)-u_exact_n(alpha))*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type_trial,alpha,0,0);
            err_element_x=err_element_x + (u_app_n(alpha)-u_exact_n(alpha))*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type_trial,alpha,1,0);
            err_element_y=err_element_y + (u_app_n(alpha)-u_exact_n(alpha))*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type_trial,alpha,0,1);
           
        end
        error_L2=error_L2+J*Gauss_weights(k)*(err_element)^2;
        error_H1=error_H1+J*Gauss_weights(k)*( err_element_x^2+err_element_y^2);
    end
        
end
if(strcmp(method,'DG')==1)
error_Edges=0;
for n=1:number_of_edges % Loop over the edges
    vertex1=Eb(1,n);%First vertex of the edge
    vertex2=Eb(2,n);%Second vertex of the edge
    T_left=Eb(3,n);% triangle on the left
    T_right=Eb(4,n);% triangle on the right
    [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),order_Gauss);
    vertices_plus=P(:,T((1:3),T_left));
    vertex_left_index=Tb((1:3),T_left);%index vertices triangle on the left
    u_app_plus=u(vertex_left_index); % get indeces of DOF inside the element.
    edge_length=norm(P(:,vertex1)-P(:,vertex2),2);
    mu=alpha_coef/edge_length;
    if(T_right~=-1)
        vertex_right_index=Tb((1:3),T_right);% index vertices triangle on the right
        vertices_minus=P(:,T((1:3),T_right)); %position of the geometric vertices of the triangles
        u_app_minus=u(vertex_right_index); % get indeces of DOF inside the element.
        for k=1:size(Gauss_nodes,2) % for each gauss node
            int_square_plus=0;% term \int_e (uapp+)^2
            int_square_minus=0; % term \int_e (uapp-)^2
            int_cross=0;% construction derivative_y error inside each element
            uapp_plus=0;
            uapp_minus=0;
            for alpha=1:number_of_local_basis_trial %loop over trial functions on the element for trial FE space.
                uapp_plus=uapp_plus +(u_app_plus(alpha))*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices_plus,basis_type_trial,alpha,0,0);
                uapp_minus=uapp_minus +(u_app_minus(alpha))*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices_minus,basis_type_trial,alpha,0,0);
            end
            int_square_plus=int_square_plus + Gauss_weights(k)*uapp_plus^2;
            int_square_minus=int_square_minus + Gauss_weights(k)*uapp_minus^2;
            int_cross=int_cross+ Gauss_weights(k)*uapp_plus*uapp_minus;
        end
        error_Edges=error_Edges + mu*(int_square_plus+int_square_minus-2*int_cross);  
    elseif(T_right==-1)
        for k=1:size(Gauss_nodes,2) % for each gauss node
            int_square_plus=0;% term \int_e (uapp+)^2
            uapp_plus=0;
            for alpha=1:number_of_local_basis_trial %loop over trial functions on the element for trial FE space.
                uapp_plus=uapp_plus +(u_app_plus(alpha))*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices_plus,basis_type_trial,alpha,0,0);
            end
            int_square_plus=int_square_plus + Gauss_weights(k)*uapp_plus^2;
         end
        error_Edges=error_Edges + mu*(int_square_plus);  
    end

end
error_DG=error_H1+error_Edges;
error_DG=sqrt(error_DG);
end
error_L2=sqrt(error_L2);
error_H1=sqrt(error_H1);
end