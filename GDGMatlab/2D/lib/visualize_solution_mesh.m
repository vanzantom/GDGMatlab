%-------------------------------------------------------------------------
% The function visualize_solution_mesh receives:
% Pb,Tb and T: matrices with information on FE space(Pb,Tb) and mesh(T). 
% u: solution
% Nsub: vector. Nsub(jj) last index of DOF in subdomain jj
% flag:  a scalar variable. If flag==1, it plots the solution
% if flag==2, it shows the mesh elements. 
%If flag==3, it shows the subdomain  decomposition

% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function []=visualize_solution_mesh(T,Pb,Tb,u,Nsub,flag)

n=length(u);
Tb=[Tb;T(4,:)];
Pb=Pb(:,1:n);
Nsubdomains=length(Nsub);
if(flag==1)
    pdeplot(Pb,[],Tb,'xydata',u,'xystyle','flat',...
           'zdata',u,'zstyle','continuous','colorbar','on');
elseif(flag==2)
u=zeros(n,1);
for i=1:size(Tb,2)
   u(Tb((1:3),i))=i;
end
elseif(flag==3)
u=zeros(n,1);

for i=1:size(Tb,2)
    for jj=1:Nsubdomains
        if(Tb(4,i)==jj)
            u(Tb((1:3),i))=jj;
        end
    end
end
end

pdeplot(Pb,[],Tb,'xydata',u,'xystyle','flat',...
           'zdata',u,'zstyle','continuous','colorbar','on');
end