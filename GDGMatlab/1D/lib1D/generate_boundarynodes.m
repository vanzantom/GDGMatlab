function boundary=generate_boundarynodes(Tb,a_label,b_label)
%generate_boundarynodes generates a vector with the indexes of the
%boundary nodes and specify their Type
%-1 Dirichlet; -2 Neumann, -3 Robin. (Robin still to implement

boundary=zeros(2,2);
number_of_elements=size(Tb,2);
boundary(1,1)=a_label;
boundary(1,2)=b_label;
boundary(2,1)=Tb(1,1);
boundary(2,2)=Tb(end,number_of_elements);
end