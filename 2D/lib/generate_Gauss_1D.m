function [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(vertices1,vertices2,type)
%the function generate_Gauss_1D generates the quadrature nodes on the
%line defined by vertices1 and  vertices2. It contains already the
%Jacobian which in 1D is equalt to edge_length.
edge_length=norm(vertices1-vertices2,2);
if type==1
Gauss_nodes(1)=vertices1;
Gauss_nodes(2)=vertices2;
Gauss_weights(1)=1/2;
Gauss_weights(2)=1/2;
Gauss_weights=Gauss_weights*edge_length;
elseif type>=2
Gauss_nodes(:,1)=vertices1;
Gauss_nodes(:,2)=(vertices1+vertices2)/2;
Gauss_nodes(:,3)=vertices2;
Gauss_weights(1)=1/3;
Gauss_weights(2)=4/3;
Gauss_weights(3)=1/3;
Gauss_weights=Gauss_weights*edge_length/2;
end