function [R]=Prolongatortilde(NgammaC)
    Numbercoarseelements=NgammaC/2;
    B=zeros(4,2);
    B(1,1)=1;
    B(2,1)=0;
    B(2,2)=0;
    B(3,1)=0;
    B(3,2)=0;
    B(4,2)=1;
    R=[];
    for k=1:Numbercoarseelements
        R=blkdiag(R,B);
    end
end