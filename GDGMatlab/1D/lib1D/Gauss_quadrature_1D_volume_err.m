function result=Gauss_quadrature_1D_volume_err(f_exact,uapp,Tb,Gauss_nodes,Gauss_weights,vertices,basis_index_trial)
h=vertices(end)-vertices(1);
number_of_local_basis_trial=size(Tb,1);
Gpn=length(Gauss_nodes);% number of Gauss points.
r=0;
for k=1:Gpn %loop on the Gaussian nodes
         upp=0;
         for j=1:number_of_local_basis_trial
             upp=upp+uapp(j)*1/h*FE_reference_basis_1D((Gauss_nodes(k)-vertices(1))/h,basis_index_trial,j,1);
         end
        r=r+Gauss_weights(k)*(feval(f_exact,Gauss_nodes(k))-upp)^2;
end
result=r;
end