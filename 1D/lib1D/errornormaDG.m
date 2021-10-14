%-------------------------------------------------------------------------
% errornormaDG receives 
% P,T: mesh matrices ,
% Tb: matrices with information on trial FEM space
% uex_der: function handle for the derivative of the exact solution
% u: array containing discrete solution
% basis_type_trial: finite element space type trial
% penalty: penalization parameter
% assemble_rhs_1D returns: 
% result: error between discrete and exact solution in norm DG.


% author: Tommaso Vanzan
%-------------------------------------------------------------------------


function result=errornormaDG(P,T,Tb,uex_der,u,basis_type_trial,penalty)
number_elements=size(Tb,2);
result=0;
%cycle over the element to compute the 1 norm. \|uex_der-u_der\|_L2
for k=1:number_elements
    vertices=P(:,T(:,k));
    h=vertices(end)-vertices(1);
    [Gauss_nodes,Gauss_weights]=generate_Gauss(vertices,2);
    uelement=u(Tb(:,k));
    result=result+Gauss_quadrature_1D_volume_err(uex_der,uelement,Tb,Gauss_nodes,Gauss_weights,vertices,basis_type_trial);
end
% Now compute edge term |[e]|^2/h= (1/h) *[uex(x+)-uapp(x+)-uex(x-)+uapp(x-)]^2
for k=1:number_elements-1
    uelement1=u(Tb(:,k));
    uelement2=u(Tb(:,k+1));
    result=result+penalty*(uelement1(end)-uelement2(1))^2;
end
result=sqrt(result);
end