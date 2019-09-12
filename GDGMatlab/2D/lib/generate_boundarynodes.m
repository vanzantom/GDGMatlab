function boundary=generate_boundarynodes(Pb)
%generate_boundarynodes generates a vector with the indexes of the
%boundary nodes and specify their Type
%-1 Dirichlet; -2 Neumann, -3 Robin.
j=1;
for k=1:size(Pb,2)
    if( Pb(1,k)==0 || Pb(1,k)==1 || Pb(2,k)==1 || Pb(2,k)==0)
        boundary(1,k)=-1;
        boundary(2,k)=k;
    end
end