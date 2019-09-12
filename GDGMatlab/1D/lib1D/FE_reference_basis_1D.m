function result=FE_reference_basis_1D(x_hat,basis_type,basis_index,der)
%FE_reference_basis_1D creates, according to the basis type the FE
%functions on a reference elemennt [0,1]
if basis_type==101
    if basis_index==1
        if der==0
           result=1-x_hat;
        elseif der==1
                result=-1;
        elseif der>1==0
                der=0;
        else
         fprintf('error in order derivative');
        end
    end
    if basis_index==2
        if der==0
           result=x_hat;
        elseif der==1
                result=1;
        elseif der>1==0
                result=0;
        else
         fprintf('error in order derivative');
        end
    end

elseif basis_type==102
    if basis_index==1 % function 1 on the left boundary.
        if der==0
           result=2*x_hat.^2-3*x_hat+1;
        elseif der==1
                result=4*x_hat-3;
        elseif der==2
                result=4;
        elseif der>2
        else
         fprintf('error in order derivative');
        end
    end
    if basis_index==3 % function 1 on the right boundary.
        if der==0
           result=2*x_hat.^2-x_hat;
        elseif der==1
                result=4*x_hat-1;
        elseif der==2
                result=4;
        elseif der>2
        else
         fprintf('error in order derivative');
        end
    end
    if basis_index==2  % function 1 in the middle.
        if der==0
           result=-4*x_hat.^2+4*x_hat;
        elseif der==1
                result=-8*x_hat+4;
        elseif der==2
                result=-8;
        elseif der>2
        else
         fprintf('error in order derivative');
        end
    end
    
else
    fprintf('no such type stored');
end