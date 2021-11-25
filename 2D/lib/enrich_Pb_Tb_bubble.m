function [Pb,Tb]=enrich_Pb_Tb_bubble(Pb,Tb,basis_type)

n=size(Pb,2);
for k=1:size(Tb,2)
    vertices=Pb(:,Tb(1:3,k));
    G=mean(Pb(:,Tb(1:3,k)),2);
    Pb(1,n+k)=G(1);
    Pb(2,n+k)=G(2);
    if(basis_type==2011)
        Tb(4,k)=n+k;
    elseif basis_type==2021
        Tb(7,k)=n+k;
    end
end


end