function [R]=Prolongator(NgammaC)
    Numbercoarseelements=NgammaC/2;
    B=zeros(4,2);
    B(1,1)=1;
    B(2,1)=1/2;
    B(2,2)=1/2;
    B(3,1)=1/2;
    B(3,2)=1/2;
    B(4,2)=1;
    R=[];
    for k=1:Numbercoarseelements
        R=blkdiag(R,B);
    end
end