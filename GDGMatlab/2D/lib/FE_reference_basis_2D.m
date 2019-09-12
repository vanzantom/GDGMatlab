function result=FE_reference_basis_2D(x_hat,y_hat,basis_type,basis_index,der_x,der_y)
%FE_reference_basis_1D creates, according to the basis type the FE
%functions on a reference elemennt [0,1]



if basis_type==201 || basis_type==2010 % Continous P1(201) and DG P1 basis_type==2010 (DOF in the vertices)
    if basis_index==2
        if der_x==0 && der_y==0
           result=x_hat;
        elseif der_x==1 && der_y==0
                result=1;
        else
                result=0;
        end
    
    elseif basis_index==3
        if der_x==0 && der_y==0
           result=y_hat;
        elseif der_x==0 && der_y==1
                result=1;
        else
                result=0;
        end
    
    elseif basis_index==1
        if der_x==0 && der_y==0
           result=1-x_hat-y_hat;
        elseif der_x==0 && der_y==1
                result=-1;
        elseif der_x==1 && der_y==0
                result=-1;
        else
                result=0;
        end
    end
% elseif basis_type==2010 % Discontinous P1
%     if basis_index==2
%         if der_x==0 && der_y==0
%            result=x_hat;
%         elseif der_x==1 && der_y==0
%                 result=1;
%         else
%                 result=0;
%         end
%     
%     elseif basis_index==3
%         if der_x==0 && der_y==0
%            result=y_hat;
%         elseif der_x==0 && der_y==1
%                 result=1;
%         else
%                 result=0;
%         end
%     
%     elseif basis_index==1
%         if der_x==0 && der_y==0
%            result=1-x_hat-y_hat;
%         elseif der_x==0 && der_y==1
%                 result=-1;
%         elseif der_x==1 && der_y==0
%                 result=-1;
%         else
%                 result=0;
%         end
%     end
else
    fprintf('no such type stored');
end