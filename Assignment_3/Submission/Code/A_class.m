classdef A_class   
    properties        
        theta_list; 
        num_rows;
        num_cols;
    end    
    methods      

       function object = A_class(m,n,theta_list)            
            object.num_rows=m;
            object.num_cols=n;
            object.theta_list=theta_list;            
       end

       function product = mtimes(object,matrix)
%            As discussed in class, we use the radon(idct2()) function to
%            express multiplication by A. The final reshape is done as we
%            know that the input matrix will always be in the vectorised
%            format.
            product=reshape(radon(idct2(reshape(matrix,sqrt(object.num_cols),[])),object.theta_list), [],1);
       end  
    end    
end



