%-------------------------------------------------------------------------
% The function compute_errors computes the L_2 error and H1 error of the finite element solution
% It receives:
% u: solution of the discrete problem
% uexvec: the exact solution vector
% P and T, Tb, Eb:   the mesh info and FE info.
%Order gauss quadrature
%basis_type: type of Finite element space.


% Attention: DG norm is computed correctly only if uex=0 on boundary.
% Because, up to now, it is assumed [[uex]]=0 on the interior edges but not on the boundary unless uex=0.


% author: Tommaso Vanzan
%-------------------------------------------------------------------------


function [error_L2,error_H1,error_DG]=compute_errors_Stokes(u,uex,P,T,Tb,Eb,order_Gauss,basis_type)

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
            err_element=err_element +(u_app_n(alpha)-u_exact_n(alpha))*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type,alpha,0,0);
            err_element_x=err_element_x + (u_app_n(alpha)-u_exact_n(alpha))*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type,alpha,1,0);
            err_element_y=err_element_y + (u_app_n(alpha)-u_exact_n(alpha))*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type,alpha,0,1);
           
        end
        error_L2=error_L2+J*Gauss_weights(k)*(err_element)^2;
        error_H1=error_H1+J*Gauss_weights(k)*( err_element_x^2+err_element_y^2);
    end
        
end
error_L2=sqrt(error_L2);
error_H1=sqrt(error_H1);
end