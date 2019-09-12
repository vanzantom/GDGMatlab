function [Gauss_nodes,Gauss_weights]=generate_Gauss(vertices,type)
%generate_Gauss generates Netwon-Cotes formulae on the interval described
%by vertices.
%type=1=> trapeizodal rule
%type=2=>Simpson rule
h=(vertices(2)-vertices(1));
if type==1
Gauss_nodes(1)=vertices(1);
Gauss_nodes(2)=vertices(end);
Gauss_weights(1)=1/2;
Gauss_weights(2)=1/2;
Gauss_weights=Gauss_weights*h;
elseif type==2
Gauss_nodes(1)=vertices(1);
Gauss_nodes(2)=(vertices(end)+vertices(1))/2;
Gauss_nodes(3)=vertices(end);
Gauss_weights(1)=1/3;
Gauss_weights(2)=4/3;
Gauss_weights(3)=1/3;
Gauss_weights=Gauss_weights*h/2;

end