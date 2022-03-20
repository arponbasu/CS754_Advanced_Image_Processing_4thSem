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
            product=reshape(radon(idct2(reshape(matrix,sqrt(object.num_cols),sqrt(object.num_cols))),object.theta_list), [],1);
       end  
    end    
end



