function result=errornormaDG(P,T,Tb,uex_der,u,basis_type_trial,penalty)
number_elements=size(Tb,2);
result=0;
%cycle over the element to compute the 1 norm.
for k=1:number_elements
    vertices=P(:,T(:,k));
    [Gauss_nodes,Gauss_weights]=generate_Gauss(vertices,2);
    uelement=u(Tb(:,k));
    result=result+Gauss_quadrature_1D_volume_err(uex_der,uelement,Tb,Gauss_nodes,Gauss_weights,vertices,basis_type_trial);
end
% Ora calcolo |[e]|^2/h= (1/h) *[uex(x+)-uapp(x+)-uex(x-)+uapp(x-)]^2
for k=1:number_elements-1
    uelement1=u(Tb(:,k));
    uelement2=u(Tb(:,k+1));
    result=result+penalty*(uelement1(end)-uelement2(1))^2;
end
result=sqrt(result);
end