classdef At_class   
    properties        
        theta_list; 
        num_rows;
        num_cols;
    end    
    methods      

       function object = At_class(m,n,theta_list)            
            object.num_rows=m;
            object.num_cols=n;
            object.theta_list=theta_list;            
       end

       function product = mtimes(object,matrix)
            product=reshape(dct2(iradon(reshape(matrix,object.num_cols,[]),object.theta_list,'nearest','Ram-Lak',1,sqrt(object.num_rows))),[],1);
       end  
    end    
end



